"""Compare Rust alpha draws with actual legacy alisim.py sample_scale draws."""
import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import pickle
import random
import subprocess
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler", type=Path, required=True, help="Built sample_gamma example")
    parser.add_argument("--python-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("validation/empirical_gamma"))
    parser.add_argument("--count", type=int, default=20000)
    args = parser.parse_args()
    if args.count < 458:
        parser.error("at least 458 draws required for the two-sample DKWM bound")
    bundle = Path(__file__).resolve().parents[1] / "data/empirical_gamma.tsv"
    text = bundle.read_text()
    headers = dict(line[2:].split("=", 1) for line in text.splitlines() if line.startswith("# ") and "=" in line)
    payload = "".join(line + "\n" for line in text.splitlines() if not line.startswith("#"))
    source = (args.python_script.parent / "hogenom_alphas.txt").read_bytes()
    assert hashlib.sha256(source).hexdigest() == headers["source_sha256"]
    assert hashlib.sha256(payload.encode()).hexdigest() == headers["values_sha256"]
    alphas = pickle.loads(source)
    assert np.array_equal(np.array(alphas, dtype=np.float64), np.array([float(x) for x in payload.splitlines()]))
    assert len(alphas) == int(headers["count"])
    spec = importlib.util.spec_from_file_location("legacy_alisim", args.python_script)
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    random.seed(20261102)
    np.random.seed(20261103)
    python = np.array([legacy.sample_scale(alphas) for _ in range(args.count)])
    command = [str(args.sampler.resolve()), "--seed", "42", "--count", str(args.count)]
    rust_text = subprocess.check_output(command, text=True)
    rust = np.array([float(x) for x in rust_text.splitlines()])
    assert len(rust) == args.count and np.all(rust >= 0.05)
    points = np.sort(np.concatenate([python, rust]))
    d = float(np.max(np.abs(np.searchsorted(np.sort(python), points, side="right") / len(python)
                               - np.searchsorted(np.sort(rust), points, side="right") / len(rust))))
    threshold = math.sqrt(math.log(2 / 0.01) / args.count)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, draws in [("python", python), ("rust", rust)]:
        with gzip.open(args.output / f"{name}.csv.gz", "wt") as stream:
            stream.write("alpha\n" + "".join(f"{x!r}\n" for x in draws.tolist()))
    result = {
        "count_per_backend": args.count, "ks_d": d, "threshold_at_0.01": threshold,
        "p_upper": min(1.0, 2 * math.exp(-args.count * d * d)), "detected_difference": d > threshold,
        "empirical_source_sha256": headers["source_sha256"], "values_sha256": headers["values_sha256"],
        "python_script_sha256": hashlib.sha256(args.python_script.read_bytes()).hexdigest(),
        "sampler_sha256": hashlib.sha256(args.sampler.read_bytes()).hexdigest(),
        "seeds": {"rust": 42, "python_random": 20261102, "numpy": 20261103},
        "numpy_version": np.__version__,
        "floor_fraction": {"python": float(np.mean(python == 0.05)), "rust": float(np.mean(rust == 0.05))},
        "sampling": "Uniform observation including multiplicity, Normal(a,a/10), floor at 0.05",
        "bundle_matches_source_exactly": True,
    }
    (args.output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["detected_difference"]:
        raise SystemExit("Gamma prior comparison detected a difference")


if __name__ == "__main__":
    main()
