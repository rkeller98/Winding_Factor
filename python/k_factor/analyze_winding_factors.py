"""Reproduce MATLAB's conventional contour estimator, without fitting.

Reference: matlab/common/estimate_contour_factor.m and paper section 3.2.
Only psi_d is intersected; iq is not constrained. Currents are in A, flux in Vs,
speed references in mechanical rpm, and k in A/A. Importing runs no analysis.
"""

from pathlib import Path
import sys

import numpy as np
from scipy.spatial import Delaunay, QhullError

ROOT = Path(__file__).resolve().parents[2]
if __package__ in (None, ""):
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "Resources/MeasEval/python"))

from python.k_factor.prepare_reduction_data import ReductionDataPreparation
from vice_measeval.backend.numeric.geometry.levelset.edge_intersector import (
    batched_edge_intersections,
)


def estimate_contour_factor(data):
    """Estimate k for one speed using prepared, named measurement arrays.

    Groups are ordered by nominal ie_ref; geometry and interpolated excitation
    use measured currents. Return (group, level) counts/means and (pair, level)
    ratios, validity and exclusion reasons. edges contains the filtered topology
    in each group's original point indexing, for inspecting the geometry.

    Shared vertices count once per intersected edge, without deduplication.
    An absent contour leaves NaN means and never bridges over a missing group.
    MATLAB excludes nonfinite ratios and abs(ratio)<eps; it does NOT impose a
    denominator tolerance or discard negative denominators. The final arithmetic
    mean weights each surviving pair/level ratio equally.
    """
    speeds = np.unique(data["rpm_ref"])
    if len(speeds) != 1 or not np.isfinite(speeds[0]):
        raise ValueError("estimate_contour_factor requires exactly one known rpm_ref")
    groups = np.unique(data["ie_ref"][np.isfinite(data["ie_ref"])])
    psi_min, psi_max = np.min(data["psi_d"]), np.max(data["psi_d"])
    span = psi_max - psi_min
    # Value-span fractions, not quantiles; both endpoints are included.
    levels = np.linspace(psi_min + 0.25 * span, psi_min + 0.75 * span, 100)
    shape = (len(groups), len(levels))
    pair_shape = (max(len(groups) - 1, 0), len(levels))
    result = {
        "rpm_ref": float(speeds[0]),
        "k": np.nan,
        "status": "not_determinable",
        "reason": "",
        "ie_ref": groups,
        "levels": levels,
        "intersection_counts": np.zeros(shape, dtype=int),
        "mean_id": np.full(shape, np.nan),
        "mean_iq": np.full(shape, np.nan),
        "mean_ie": np.full(shape, np.nan),
        "delta_ie": np.full(pair_shape, np.nan),
        "ratios": np.full(pair_shape, np.nan),
        "valid": np.zeros(pair_shape, dtype=bool),
        "exclusion_reason": np.full(pair_shape, "missing_contour", dtype="U32"),
        "group_errors": [""] * len(groups),
        "edges": [np.empty((0, 2), dtype=int) for _ in groups],
    }
    if not np.isfinite(data["ie_ref"]).all():
        result["reason"] = (
            "Missing ie_ref: excitation groups cannot be assigned reliably for this speed."
        )
        result["exclusion_reason"][:] = "missing_ie_ref"
        return result

    for g, excitation in enumerate(groups):
        use = data["ie_ref"] == excitation
        samples = np.column_stack(
            [data[name][use] for name in ("id_meas", "iq_meas", "ie_meas", "psi_d")]
        )
        basis = samples[:, :2]
        try:
            # Original 2-D coordinates, unique undirected edges, strict length
            # test on those unique edges: exactly the MATLAB estimator's rule.
            triangles = Delaunay(basis).simplices
            edges = np.vstack([triangles[:, pair] for pair in ((0, 1), (0, 2), (1, 2))])
            edges = np.unique(np.sort(edges, axis=1), axis=0)
            length_squared = np.sum(
                (basis[edges[:, 0]] - basis[edges[:, 1]]) ** 2, axis=1
            )
            edges = edges[length_squared < 9 * np.mean(length_squared)]
            result["edges"][g] = edges
            # Reuse the existing MATLAB-arithmetic kernel for ALL channels.
            # Unlike EdgeIntersector, this does not canonically reorder points,
            # rescale the length test or change interpolation arithmetic.
            hits = batched_edge_intersections(
                samples,
                search_column=3,
                levels=levels,
                edges=edges,
                index_base=0,
                numeric_policy="matlab-ieee",
            )
        except (QhullError, ValueError) as error:
            result["group_errors"][g] = (
                f"ie_ref={excitation:g} A: {str(error).splitlines()[0]}"
            )
            continue
        offsets = hits["level_offsets"]
        result["intersection_counts"][g] = np.diff(offsets)
        for l, (start, end) in enumerate(zip(offsets[:-1], offsets[1:])):
            if start == end:
                continue
            # Equal weight per edge hit, including shared-vertex repetitions.
            # Measured ie was interpolated with the SAME lambda as id and iq.
            mean = hits["samples"][start:end, :3].mean(axis=0)
            (
                result["mean_id"][g, l],
                result["mean_iq"][g, l],
                result["mean_ie"][g, l],
            ) = mean

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        delta_id = np.diff(result["mean_id"], axis=0)
        delta_ie = np.diff(result["mean_ie"], axis=0)
        ratios = np.abs(delta_id) / delta_ie
    result["delta_ie"] = delta_ie
    # MATLAB eps at 1 is applied to the ratio, NOT the denominator.
    valid = np.isfinite(ratios) & (np.abs(ratios) >= np.finfo(float).eps)
    reason = result["exclusion_reason"]
    reason[:] = ""
    reason[~np.isfinite(ratios)] = "nonfinite_ratio"
    reason[delta_ie == 0] = "zero_denominator"
    reason[np.isfinite(ratios) & ~valid] = "below_eps"
    counts = result["intersection_counts"]
    reason[(counts[:-1] == 0) | (counts[1:] == 0)] = "missing_contour"
    result["valid"] = valid
    result["ratios"] = np.where(valid, ratios, np.nan)
    if any(result["group_errors"]):
        # MATLAB aborts failed triangulations. Do not aggregate a different
        # estimator using only the remaining geometry groups.
        result["reason"] = "Geometry failed: " + "; ".join(
            filter(None, result["group_errors"])
        )
    elif len(groups) < 2:
        result["reason"] = "At least two nominal excitation groups are required."
    elif not valid.any():
        result["reason"] = "No valid adjacent-group contour factors."
    else:
        result["k"] = float(np.mean(ratios[valid]))
        result["status"] = "ok"
    return result


def analyze_winding_factors(data):
    """Return results keyed by rpm_ref; never combine different speeds.

    Input is the dictionary from ReductionDataPreparation.run(). Missing rpm_ref
    makes dataset grouping impossible and raises ValueError. Missing ie_ref
    makes the affected speed explicitly not_determinable.
    """
    if not len(data["rpm_ref"]):
        raise ValueError("No prepared measurement points are available.")
    unknown = np.count_nonzero(~np.isfinite(data["rpm_ref"]))
    if unknown:
        raise ValueError(
            f"Missing rpm_ref for {unknown} points: reliable speed grouping is impossible."
        )
    results = {}
    for rpm in dict.fromkeys(data["rpm_ref"]):
        use = data["rpm_ref"] == rpm
        results[float(rpm)] = estimate_contour_factor(
            {name: values[use] for name, values in data.items()}
        )
    return results


def main(argv=None):
    """Read supplied Meas/raw JSON files (default: research exports of Set1–3)."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "files", nargs="*", type=Path, help="Measurement .json.gz files"
    )
    args = parser.parse_args(argv)
    files = args.files or [
        ROOT
        / f"Resources/DataSets/Set{s}/extracted/evaluation/research_flux_measurements.json.gz"
        for s in (1, 2, 3)
    ]
    results = {}
    for filename in files:
        data = ReductionDataPreparation(filename).run()
        try:
            results[str(filename)] = analyze_winding_factors(data)
        except ValueError as error:
            print(f"{filename}: not determinable: {error}")
            results[str(filename)] = {"error": str(error)}
            continue
        for rpm, result in results[str(filename)].items():
            n_valid = np.count_nonzero(result["valid"])
            print(
                f"{filename}: {rpm:g} rpm, {len(result['ie_ref'])} groups, "
                f"{n_valid}/{result['valid'].size} valid factors, k={result['k']:.12g} A/A"
            )
            if result["status"] != "ok":
                print(f"  not determinable: {result['reason']}")
    return results


if __name__ == "__main__":
    main()
