"""Write reactant.con (bcc) and product.con (fcc) for the Bain-path band.

Run: python make_endpoints.py

Both endpoints are the two-atom body-centred tetragonal cell diag(a, a, c),
atoms at (0, 0, 0) and (1/2, 1/2, 1/2), repeated 5 x 5 x 4 (200 atoms).
The lattice constants are the relaxed minima of LJ (u0 1 eV, psi 1 A) with
the C^2 switch over [2.0, 2.2] A, from validation/bain_lj_reference.py in
rgsaddle: bcc a = c = 1.228441371 A, fcc (as bct) a = 1.100691058 A,
c = 1.556612222 A = sqrt(2) a. Each cell width stays above 4.4 A, twice
the cutoff, so rgpot's minimum-image pair list sees every pair.
"""

from pathlib import Path

REPS = (5, 5, 4)
BCC = (1.228441371, 1.228441371)
FCC = (1.100691058, 1.556612222)


def write(path, a, c):
    nx, ny, nz = REPS
    rows = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for b in (0.0, 0.5):
                    rows.append(((i + b) * a, (j + b) * a, (k + b) * c))
    lines = [
        "Bain path endpoint, LJ cut at 2.2 A",
        "",
        f"{nx * a:.9f}\t{ny * a:.9f}\t{nz * c:.9f}",
        "90.000000\t90.000000\t90.000000",
        "",
        "",
        "1",
        f"{len(rows)}",
        "39.948000",
        "Ar",
        "Coordinates of Component 1",
    ]
    for n, (x, y, z) in enumerate(rows):
        lines.append(f"{x:22.15f} {y:22.15f} {z:22.15f} 0 {n:5d}")
    Path(path).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    write("reactant.con", *BCC)
    write("product.con", *FCC)
