import json
import sys
from pathlib import Path

# This helps find leaks, since potentials with unexpected force call counts
# indicate unpaired new/delete or missing cleanup.


def analyze_potcalls(json_path):
    data = json.loads(Path(json_path).read_text())
    total_force_calls = sum(r["force_calls"] for r in data)
    total_instances = len(data)

    # Group by type
    by_type = {}
    for r in data:
        t = r["type"]
        if t not in by_type:
            by_type[t] = {"instances": 0, "force_calls": 0}
        by_type[t]["instances"] += 1
        by_type[t]["force_calls"] += r["force_calls"]

    print(f"Total instances destroyed: {total_instances}")
    print(f"Total force calls: {total_force_calls}")
    print()
    for t, stats in sorted(by_type.items()):
        print(f"  {t}: {stats['instances']} instances, {stats['force_calls']} force calls")

    # Flag instances with suspiciously low force calls (potential diagnostic)
    low_call = [r for r in data if r["force_calls"] <= 1]
    if low_call:
        print(f"\nInstances with <= 1 force call ({len(low_call)}):")
        for r in low_call:
            print(f"  id={r['id']} type={r['type']} force_calls={r['force_calls']}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "_potcalls.json"
    analyze_potcalls(path)
