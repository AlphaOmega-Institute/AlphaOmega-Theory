#!/usr/bin/env python3
"""Run every verification check and report pass/fail.

    python checks/run_all.py            # standard run
    python checks/run_all.py --full     # also run the slow checks
    python checks/run_all.py -k orbital # only checks whose name contains "orbital"

Exit status is 0 when every selected check passes and 1 otherwise.
"""
import argparse
import importlib.util
import os
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(path):
    name = ".".join(path.relative_to(HERE).with_suffix("").parts)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true", help="also run the slow checks")
    parser.add_argument("-k", default="", help="only run checks whose name contains this text")
    args = parser.parse_args()
    if args.full:
        os.environ["CHECKS_FULL"] = "1"
    sys.path.insert(0, str(HERE))
    from lib import SkipCheck

    passed = failed = skipped = 0
    start = time.time()
    for path in sorted(HERE.rglob("test_*.py")):
        rel = path.relative_to(HERE).as_posix()
        try:
            module = load(path)
        except Exception:
            print(f"FAIL  {rel} (import)")
            traceback.print_exc()
            failed += 1
            continue
        for name, fn in vars(module).items():
            if not (name.startswith("test_") and callable(fn)):
                continue
            label = f"{rel}::{name}"
            if args.k and args.k not in label:
                continue
            t0 = time.time()
            try:
                fn()
            except SkipCheck as exc:
                print(f"SKIP  {label}  ({exc})")
                skipped += 1
            except Exception:
                print(f"FAIL  {label}")
                traceback.print_exc()
                failed += 1
            else:
                print(f"PASS  {label}  ({time.time() - t0:.1f}s)")
                passed += 1
    print(f"\n{passed} passed, {failed} failed, {skipped} skipped in {time.time() - start:.0f}s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
