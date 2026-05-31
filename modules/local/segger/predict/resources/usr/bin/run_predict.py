#!/usr/bin/env python3
"""
Run segger predict with spatialxe-specific preprocessing.

Wraps segger's predict_fast.py with:
  - GPU enumeration (replaces inline python3 -c torch check)
  - WORKAROUND: patch predict_parquet.py at runtime to add torch.no_grad() for ~30-50% VRAM savings
  - WORKAROUND: seed random.choice for deterministic GPU assignment (avoids stochastic OOM)

Both WORKAROUNDs should be removable once the patches are upstreamed to segger.
"""

import argparse
import os
import subprocess
import sys


SEGGER_CLI = "/workspace/segger_dev/src/segger/cli/predict_fast.py"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--models-dir", required=True)
    p.add_argument("--segger-data-dir", required=True)
    p.add_argument("--transcripts-file", required=True)
    p.add_argument("--benchmarks-dir", required=True)
    p.add_argument("--batch-size", type=int, required=True)
    p.add_argument("--use-cc", required=True)
    p.add_argument("--knn-method", required=True)
    p.add_argument("--num-workers", type=int, required=True)
    args, extra = p.parse_known_args()
    return args, extra


def detect_gpus():
    """Return comma-separated list of available CUDA device ids (or "0" if none)."""
    import torch

    print("=== GPU Detection (SEGGER_PREDICT) ===")
    print(f"PyTorch CUDA available: {torch.cuda.is_available()}")
    n = torch.cuda.device_count()
    print(f"CUDA device count: {n}")
    print("======================================")
    if n > 0:
        return ",".join(str(i) for i in range(n))
    return "0"


def patch_predict_parquet():
    """
    WORKAROUND: patch segger.prediction.predict_parquet at runtime.

    Avoids rebuilding the segger Docker image. Two patches:
      1. Add torch.no_grad() to disable gradient graphs during inference (~30-50% VRAM savings).
      2. Seed random for deterministic GPU assignment (avoids stochastic OOM).

    Remove this function once the patches are upstreamed to segger.
    """
    import segger.prediction.predict_parquet as m

    pred_py = m.__file__
    print(f"Patching {pred_py}: torch.no_grad() + round-robin GPU assignment")
    # Use sed via subprocess for in-place edit (matches the original behavior exactly)
    subprocess.run(
        [
            "sed",
            "-i",
            "s/with cp.cuda.Device(gpu_id):/with cp.cuda.Device(gpu_id), torch.no_grad():/",
            pred_py,
        ],
        check=True,
    )
    subprocess.run(
        [
            "sed",
            "-i",
            "s/gpu_id = random.choice(gpu_ids)/random.seed(0); gpu_id = random.choice(gpu_ids)/",
            pred_py,
        ],
        check=True,
    )


def run_segger_cli(args, extra, gpu_ids):
    cmd = [
        "python3",
        SEGGER_CLI,
        "--models_dir",
        args.models_dir,
        "--segger_data_dir",
        args.segger_data_dir,
        "--transcripts_file",
        args.transcripts_file,
        "--benchmarks_dir",
        args.benchmarks_dir,
        "--batch_size",
        str(args.batch_size),
        "--use_cc",
        str(args.use_cc),
        "--knn_method",
        args.knn_method,
        "--num_workers",
        str(args.num_workers),
        "--gpu_ids",
        gpu_ids,
        *extra,
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main():
    args, extra = parse_args()

    # Limit cupy GPU memory to 80% so PyTorch has headroom for graph attention ops
    os.environ.setdefault("CUPY_GPU_MEMORY_LIMIT", "80%")
    # Belt-and-suspenders: ensure PyTorch uses expandable segments
    os.environ.setdefault(
        "PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True,max_split_size_mb:512"
    )
    # Numba cache directory
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.getcwd(), ".numba_cache"))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

    gpu_ids = detect_gpus()
    print(f"Using GPUs: {gpu_ids}")

    patch_predict_parquet()

    run_segger_cli(args, extra, gpu_ids)


if __name__ == "__main__":
    main()
