"""Export the legacy pickle lists to Rust-readable weighted text (trusted inputs)."""

import argparse
import hashlib
import pickle
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source", type=Path, help="Directory containing the two legacy pickle files"
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/empirical_diameters.tsv")
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as stream:
        stream.write("# diameter\tweight; strict filter 0.02 < diameter < 15\n")
        stream.write(
            "# HOGENOM weight 1; RAxML weight 10 (Python list replication, not diameter multiplication)\n"
        )
        for name, weight in [("hogenom_diams.txt", 1), ("raxml_diams.txt", 10)]:
            data = (args.source / name).read_bytes()
            stream.write(f"# {name} sha256={hashlib.sha256(data).hexdigest()}\n")
            for value in pickle.loads(data):
                if 0.02 < float(value) < 15:
                    stream.write(f"{float(value)!r}\t{weight}\n")


if __name__ == "__main__":
    main()
