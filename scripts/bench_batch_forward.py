#!/usr/bin/env python3
"""Benchmark batched vs sequential model.forward() for MetatomicPotential.

Run on cosmolab with the eongpu pixi env:
    pixi run -e eongpu python scripts/bench_batch_forward.py

Measures wall time and peak GPU memory for:
1. N sequential forward+backward passes (one system each)
2. 1 batched forward+backward pass (N systems)
"""
import torch, time, argparse
torch._C._jit_set_profiling_mode(False)
import metatomic.torch, metatensor.torch, vesin.torch
import numpy as np

def build_system(model, pos, types, cell, dtype, device):
    t = torch.tensor(types, dtype=torch.int32, device=device)
    p = torch.tensor(pos, dtype=torch.float64).to(dtype).to(device).requires_grad_(True)
    c = torch.tensor(cell, dtype=torch.float64).to(dtype).to(device)
    pbc = torch.tensor([False, False, False], device=device)
    sys = metatomic.torch.System(t, p, c, pbc)
    for req in model.requested_neighbor_lists():
        cutoff = req.engine_cutoff("angstrom")
        nl = vesin.torch.NeighborList(cutoff=cutoff, full_list=True)
        i_idx, j_idx, S, d, D = nl.compute(
            points=torch.tensor(pos, dtype=torch.float64),
            box=torch.tensor(cell, dtype=torch.float64),
            periodic=torch.tensor([False, False, False]),
            quantities="ijSdD", copy=True,
        )
        n = i_idx.shape[0]
        pairs = torch.stack([i_idx.to(torch.int32), j_idx.to(torch.int32)], dim=1)
        samples_vals = torch.cat([pairs, S.to(torch.int32)], dim=1).to(device)
        samples = metatensor.torch.Labels(
            ["first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"],
            samples_vals,
        )
        comp = metatensor.torch.Labels("xyz", torch.arange(3, dtype=torch.int32, device=device).reshape(-1, 1))
        prop = metatensor.torch.Labels("distance", torch.zeros(1, 1, dtype=torch.int32, device=device))
        block = metatensor.torch.TensorBlock(D.reshape(n, 3, 1).to(dtype).to(device), samples, [comp], prop)
        metatomic.torch.register_autograd_neighbors(sys, block, False)
        sys.add_neighbor_list(req, block)
    return sys, p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_path", help="Path to .pt model file")
    parser.add_argument("--nimages", type=int, default=8)
    parser.add_argument("--nreps", type=int, default=20)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    model = metatomic.torch.load_atomistic_model(args.model_path)
    device = torch.device(args.device)
    model.to(device)
    dtype = torch.float32

    types = np.array([7, 6, 8, 1, 1, 1, 1], dtype=np.int32)
    cell = np.zeros((3, 3))
    base_pos = np.array(
        [[0, 0, 0], [1.2, 0.3, 0], [2.4, -0.1, 0], [0, -1, 0.5],
         [0, 1, -0.5], [-1, 0, 0], [3.5, 0, 0]], dtype=np.float64
    )
    all_pos = []
    for k in range(args.nimages):
        p = base_pos.copy()
        p[6, 0] += k * 0.3
        p[0, 1] += k * 0.05
        all_pos.append(p)

    eval_opts = metatomic.torch.ModelEvaluationOptions()
    eval_opts.length_unit = "angstrom"
    req_out = metatomic.torch.ModelOutput()
    req_out.per_atom = False
    req_out.quantity = "energy"
    req_out.unit = "eV"
    eval_opts.outputs = {"energy": req_out}

    # Warmup
    for pos in all_pos[:2]:
        sys, p = build_system(model, pos, types, cell, dtype, device)
        out = model([sys], eval_opts, False)
        out["energy"].block(0).values.sum().backward()
    torch.cuda.synchronize()

    # Sequential
    t0 = time.perf_counter()
    for _ in range(args.nreps):
        for pos in all_pos:
            sys, p = build_system(model, pos, types, cell, dtype, device)
            out = model([sys], eval_opts, False)
            out["energy"].block(0).values.sum().backward()
        torch.cuda.synchronize()
    t_seq = (time.perf_counter() - t0) / args.nreps

    # Batched
    t0 = time.perf_counter()
    for _ in range(args.nreps):
        systems, ptensors = [], []
        for pos in all_pos:
            sys, p = build_system(model, pos, types, cell, dtype, device)
            systems.append(sys)
            ptensors.append(p)
        out = model(systems, eval_opts, False)
        out["energy"].block(0).values.sum().backward()
        torch.cuda.synchronize()
    t_batch = (time.perf_counter() - t0) / args.nreps

    print(f"Sequential ({args.nimages} images): {t_seq * 1000:.1f} ms")
    print(f"Batched    ({args.nimages} images): {t_batch * 1000:.1f} ms")
    print(f"Speedup: {t_seq / t_batch:.2f}x")

    # Memory
    if args.device == "cuda":
        torch.cuda.reset_peak_memory_stats()
        sys, p = build_system(model, all_pos[0], types, cell, dtype, device)
        model([sys], eval_opts, False)["energy"].block(0).values.sum().backward()
        torch.cuda.synchronize()
        mem1 = torch.cuda.max_memory_allocated() / 1e6

        torch.cuda.reset_peak_memory_stats()
        systems, ptensors = [], []
        for pos in all_pos:
            sys, p = build_system(model, pos, types, cell, dtype, device)
            systems.append(sys)
            ptensors.append(p)
        model(systems, eval_opts, False)["energy"].block(0).values.sum().backward()
        torch.cuda.synchronize()
        memN = torch.cuda.max_memory_allocated() / 1e6

        print(f"\nPeak GPU (1 sys):  {mem1:.0f} MB")
        print(f"Peak GPU ({args.nimages} sys): {memN:.0f} MB")
        print(f"Memory ratio: {memN / mem1:.1f}x")

if __name__ == "__main__":
    main()
