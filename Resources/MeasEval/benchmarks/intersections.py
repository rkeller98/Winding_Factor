"""Reproducible geometry/intersection timings; run from the project venv.

Before editing: python .../intersections.py --output results/measeval/before.json
After editing:  python .../intersections.py --baseline results/measeval/baseline_edge_intersector.py
Use --output to save raw measurements; no benchmark dependency is required.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time
import tracemalloc

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "Resources/MeasEval/python"))
from vice_measeval.backend.numeric.geometry.levelset import edge_intersector as current
from vice_measeval.backend.numeric.model.point_set import PointSet
from vice_measeval.backend.numeric.model.sampled_field import SampledField


def load_baseline(path):
    spec = importlib.util.spec_from_file_location("edge_baseline", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def measure(function, repeats):
    function()  # warm up
    samples = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        function()
        samples.append((time.perf_counter_ns() - start) / 1e6)
    return samples


def paired(before, after, repeats):
    before(); after()
    a, b = [], []
    for index in range(repeats):
        # Alternate order to reduce bias from CPU drift.
        order = [(before, a), (after, b)] if index % 2 else [(after, b), (before, a)]
        for function, samples in order:
            start = time.perf_counter_ns()
            function()
            samples.append((time.perf_counter_ns() - start) / 1e6)
    return {"before_ms": a, "after_ms": b,
            "ratio": float(np.median(b) / np.median(a))}


def peak(function):
    tracemalloc.start()
    result = function()
    _, highest = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return {"peak_mib": highest / 2**20,
            "pieces": len(result.pieces) if hasattr(result, "pieces") else len(result.coordinates)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=21)
    args = parser.parse_args()
    baseline = load_baseline(args.baseline) if args.baseline else None
    report = {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
              "platform": platform.platform(), "repeats": args.repeats, "cases": {}}
    source = Path(current.__file__).read_bytes()
    report["source_sha256"] = hashlib.sha256(source).hexdigest()
    if args.baseline:
        report["baseline_sha256"] = hashlib.sha256(args.baseline.read_bytes()).hexdigest()
    for dim, count in [(2, 2000), (3, 1000), (4, 250)]:
        points = np.random.default_rng(710 + dim).random((count, dim))
        basis = PointSet(points)
        field = SampledField(basis, np.column_stack((points.sum(axis=1), points[:, 0] - points[:, 1])))
        geometry = current.EdgeIntersector(basis)
        levels = np.linspace(.2 * dim, .8 * dim, 17)
        functions = {
            "prepare": lambda: current.EdgeIntersector(basis),
            "single": lambda: geometry.intersect(field, 0, levels),
            "repeated": lambda: [geometry.intersect(field, 0, value) for value in levels[:8]],
        }
        case = {"points": count, "simplices": len(geometry._simplices), "edges": len(geometry._edges)}
        if baseline:
            old = baseline.EdgeIntersector(basis)
            old_functions = {
                "prepare": lambda: baseline.EdgeIntersector(basis),
                "single": lambda: old.intersect(field, 0, levels),
                "repeated": lambda: [old.intersect(field, 0, value) for value in levels[:8]],
            }
            for name, function in functions.items():
                case[name] = paired(old_functions[name], function, args.repeats)
        else:
            for name, function in functions.items():
                case[name] = {"ms": measure(function, args.repeats)}
        if hasattr(geometry, "intersect_constraints"):
            conditions = [current.EqualityConstraint("field", 0), current.EqualityConstraint("coordinate", 0)]
            targets = [[.5 * dim, .4], [.6 * dim, .6]]
            operation = lambda: geometry.intersect_constraints(field, conditions, targets)
            start = time.perf_counter()
            first = operation()
            case["constraints_first_ms"] = 1000 * (time.perf_counter() - start)
            case["constraints_pieces"] = len(first.pieces)
            case["constraints_warm_ms"] = measure(operation, 3)
            del first
        report["cases"][str(dim)] = case
        print("D", dim, {name: round(value["ratio"], 4) if "ratio" in value else
                         round(float(np.median(value["ms"])), 3)
                         for name, value in case.items() if isinstance(value, dict)}, flush=True)
    points = np.random.default_rng(718).random((150, 3))
    geometry = current.EdgeIntersector(points)
    field = points.sum(axis=1)
    report["memory"] = {}
    for queries in (32, 256, 1024):
        report["memory"][f"edges_{queries}"] = peak(lambda: geometry.intersect(field, 0, np.linspace(.1, 2.9, queries)))
        if hasattr(geometry, "intersect_constraints"):
            constraints = [current.EqualityConstraint("coordinate", 0), current.EqualityConstraint("coordinate", 1)]
            # Empty requests expose working memory without unavoidable output growth.
            targets = np.full((queries, 2), 2.)
            report["memory"][f"constraints_empty_{queries}"] = peak(
                lambda: geometry.intersect_constraints(field, constraints, targets, chunk_size=128))
    if hasattr(geometry, "intersect_constraints"):
        for queries in (8, 32, 128):
            targets = np.column_stack((np.linspace(.2, .8, queries), np.full(queries, .5)))
            report["memory"][f"constraints_output_{queries}"] = peak(
                lambda: geometry.intersect_constraints(field, constraints, targets, chunk_size=128))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["memory"], indent=2))


if __name__ == "__main__":
    main()
