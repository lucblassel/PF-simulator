"""Export the trusted legacy alpha pickle, preserving order and multiplicity."""

import argparse
import hashlib
import math
import pickle
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Trusted hogenom_alphas.txt pickle")
    parser.add_argument("--output", type=Path, default=Path("data/empirical_gamma.tsv"))
    args = parser.parse_args()
    source = args.source.read_bytes()
    values = [float(value) for value in pickle.loads(source)]
    if not values or not all(math.isfinite(value) and value > 0 for value in values):
        raise ValueError("Expected finite positive empirical alpha values")
    payload = "".join(f"{value!r}\n" for value in values)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "# hogenom_alphas.txt; original order and multiplicity, no filtering\n"
        f"# source_sha256={hashlib.sha256(source).hexdigest()}\n"
        f"# values_sha256={hashlib.sha256(payload.encode()).hexdigest()}\n"
        f"# count={len(values)}\n" + payload
    )


if __name__ == "__main__":
    main()
