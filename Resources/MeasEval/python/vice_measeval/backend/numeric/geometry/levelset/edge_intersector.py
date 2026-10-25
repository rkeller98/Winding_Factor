"""Delaunay-edge level-set intersection for sampled numeric fields.

The module provides a reusable ``EdgeIntersector`` for fields sampled on a
fixed point basis. The basis is normalized to a validated ``PointSet``,
triangulated once, and reduced to unique geometrically valid edges.

Field input is normalized to ``SampledField`` on the same ordered basis.
Repeated intersections with one or more scalar levels therefore reuse the
prepared geometry without rebuilding the triangulation.

The implementation is independent of motor-specific MeasEval semantics.
"""

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import scipy as sp
from numpy.typing import ArrayLike, NDArray

from vice_measeval.backend.numeric.model.point_set import PointSet
from vice_measeval.backend.numeric.model.sampled_field import SampledField

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def batched_edge_intersections(
    samples: ArrayLike,
    search_column: int,
    levels: ArrayLike,
    edges: ArrayLike,
    *,
    index_base: int,
    numeric_policy: str,
) -> dict[str, NDArray]:
    """Shadow kernel for MATLAB ``edge_intersections`` with supplied topology.

    ``samples`` is an N-by-D real matrix containing coordinates and any extra
    channels. ``edges`` is E-by-2; orientation, repetitions and row order are
    preserved. No triangulation, length filtering or geometric deduplication
    occurs. ``search_column`` may select any sample column.

    ``index_base`` must be 0 or 1 and applies to the search column, input/output
    vertices, output edge rows and output level positions. ``numeric_policy``
    must be ``"matlab-ieee"``: binary64 arithmetic, inclusive lambda in [0, 1],
    start + lambda * (end - start), search column assigned the exact level.
    NaN/Inf, overflow and zero denominators follow this arithmetic unchanged.

    Levels may be scalar, a 1-D sequence or a MATLAB row/column vector (including
    empty). Results are level-major, then input-edge-major. Repeated levels are
    distinct requests. The dictionary contains only arrays: ``samples`` (M,D),
    ``edges`` (M,2), and 1-D ``levels``, ``lambdas``, ``edge_indices`` and
    ``level_indices``. ``level_offsets`` (L+1,) always contains zero-based row
    offsets, even for index_base=1, so empty levels remain recoverable by slicing.

    This function does not change the separate SciPy/PointSet-based class API.
    A single input edge with nonempty levels is rejected: the original MATLAB
    kernel fails on its vector-indexing shape in this case as well.
    """
    def real_array(value: ArrayLike, name: str) -> FloatArray:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf":
            raise TypeError(f"{name} must contain real numeric values.")
        return np.asarray(raw, dtype=np.float64)

    if isinstance(index_base, (bool, np.bool_)) or not isinstance(
        index_base, (int, np.integer)
    ) or index_base not in (0, 1):
        raise ValueError("index_base must be the integer 0 or 1.")
    if numeric_policy != "matlab-ieee":
        raise ValueError("numeric_policy must be 'matlab-ieee'.")
    x = real_array(samples, "samples")
    if x.ndim != 2 or x.shape[1] == 0:
        raise ValueError("samples must have shape (N, D) with D >= 1.")
    if isinstance(search_column, (bool, np.bool_)) or not isinstance(
        search_column, (int, np.integer)
    ):
        raise TypeError("search_column must be an integer.")
    column = int(search_column) - index_base
    if not 0 <= column < x.shape[1]:
        raise IndexError("search_column is outside samples.")
    level = real_array(levels, "levels")
    if level.ndim > 2 or (level.ndim == 2 and min(level.shape) > 1):
        raise ValueError("levels must be a scalar or vector.")
    level = level.reshape(-1)
    raw_edges = real_array(edges, "edges")
    if raw_edges.ndim != 2 or raw_edges.shape[1] != 2:
        raise ValueError("edges must have shape (E, 2).")
    if raw_edges.shape[0] == 1 and level.size:
        raise ValueError("MATLAB edge_intersections does not support a single input edge.")
    if not np.all(np.isfinite(raw_edges) & (raw_edges == np.floor(raw_edges))):
        raise ValueError("edges must contain finite integer indices.")
    if np.any(raw_edges < index_base) or np.any(
        raw_edges >= x.shape[0] + index_base
    ):
        raise IndexError("edge vertex index is outside samples.")
    vertices = raw_edges.astype(np.int64) - index_base
    start = x[vertices[:, 0], column]
    end = x[vertices[:, 1], column]
    # Keep operation order identical to MATLAB; do not rescale overflowing
    # differences or prefilter nonfinite endpoints (Inf can yield lambda=0).
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        factors = (level[:, None] - start) / (end - start)
        valid = (factors >= 0.0) & (factors <= 1.0)
        level_rows, edge_rows = np.nonzero(valid)
        lambdas = factors[level_rows, edge_rows]
        selected = vertices[edge_rows]
        first = x[selected[:, 0]]
        output = first + lambdas[:, None] * (x[selected[:, 1]] - first)
    output[:, column] = level[level_rows]
    offsets = np.empty(level.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(np.count_nonzero(valid, axis=1), out=offsets[1:])
    return {
        "samples": output,
        "edges": selected + index_base,
        "levels": level[level_rows].copy(),
        "lambdas": lambdas,
        "edge_indices": edge_rows + index_base,
        "level_indices": level_rows + index_base,
        "level_offsets": offsets,
    }


@dataclass(frozen=True)
class EdgeIntersections:
    """Immutable results of edge-level intersections.

    Each row represents one intersection between one prepared triangulation
    edge and one requested scalar level.

    Parameters
    ----------
    coordinates : FloatArray
        Interpolated coordinates with shape ``(M, D)``, where ``M`` is the
        number of detected intersections and ``D`` is the dimension of the
        underlying ``PointSet``.
    values : FloatArray
        Interpolated field values with shape ``(M, K)``, where ``K`` is the
        number of components in the intersected ``SampledField``.
    levels : FloatArray
        Requested level associated with each intersection, shape ``(M,)``.
    edges : IntArray
        Point indices of the intersected prepared edges, shape ``(M, 2)``.
    lambdas : FloatArray
        Linear interpolation factors along the intersected edges, shape
        ``(M,)``. Values lie in the closed interval ``[0, 1]``.

    Notes
    -----
    The result preserves edge-level intersections rather than unique geometric
    points. Consequently, if a requested level passes exactly through a vertex
    shared by multiple edges, the same geometric coordinate may appear more
    than once with different edge provenance.

    All result arrays are copied and exposed as read-only views after
    construction. Mutating constructor inputs cannot change an existing result.
    """

    coordinates: FloatArray
    values: FloatArray
    levels: FloatArray
    edges: IntArray
    lambdas: FloatArray

    def __post_init__(self) -> None:
        """Own read-only copies of all result arrays."""
        for name in ("coordinates", "values", "levels", "edges", "lambdas"):
            dtype = np.int64 if name == "edges" else np.float64
            owned = np.array(getattr(self, name), dtype=dtype, copy=True)
            owned.flags.writeable = False

            view = owned.view()
            view.flags.writeable = False

            object.__setattr__(self, name, view)


class EdgeIntersector:
    """Intersect sampled field levels with edges of a fixed Delaunay basis.

    An ``EdgeIntersector`` owns the prepared geometric state required for
    repeated level-set intersections on one immutable ``PointSet`` basis.

    Array-like basis input is normalized to a validated ``PointSet`` at the
    public API boundary. The numeric core therefore operates only on validated
    and immutable geometry.

    Construction performs three steps:

    1. Canonically order the basis for deterministic triangulation.
    2. Build a Delaunay triangulation.
    3. Extract unique undirected edges and remove excessively long edges.

    Once constructed, the prepared edges can be reused for arbitrary fields
    defined on the same ordered point basis.

    Parameters
    ----------
    basis : PointSet or ArrayLike
        Point coordinates defining the triangulation basis. Array-like input is
        normalized to a validated ``PointSet`` before triangulation.
    validity_range : float, default=3.0
        Maximum accepted edge length relative to the RMS edge length of the
        triangulation. The comparison is performed using squared lengths, so
        no square root is required.

    Raises
    ------
    TypeError
        If ``validity_range`` is not a real numeric value or the supplied
        basis cannot be converted by ``PointSet`` because of an invalid type.
    ValueError
        If the supplied basis cannot be normalized to a valid ``PointSet``,
        ``validity_range`` is not finite and strictly positive, the point basis
        cannot be Delaunay triangulated, or the validity filter removes all
        triangulation edges.

    Notes
    -----
    The normalized basis is immutable through ``PointSet``. Therefore
    triangulation and edge preparation are performed only once per
    ``EdgeIntersector`` instance.

    Triangulation is performed after canonical lexicographic ordering of the
    coordinates. The resulting simplex indices are mapped back to the original
    ``PointSet`` row indices so that prepared edges remain aligned with field
    rows.
    """

    def __init__(
        self,
        basis: PointSet | ArrayLike,
        validity_range: float = 3.0,
    ) -> None:
        try:
            self._basis = basis if isinstance(basis, PointSet) else PointSet(basis)
        except (TypeError, ValueError) as exc:
            exc.add_note(
                "Error occurred while converting the EdgeIntersector basis to PointSet."
            )
            raise

        if isinstance(validity_range, (bool, np.bool_)) or not isinstance(
            validity_range,
            (int, float, np.integer, np.floating),
        ):
            raise TypeError(
                f"validity_range must be a real number; "
                f"got {type(validity_range).__name__}."
            )

        if not np.isfinite(validity_range) or validity_range <= 0:
            raise ValueError(
                f"validity_range must be finite and greater than zero; "
                f"got {validity_range}."
            )

        self._validity_range = float(validity_range)

        self._simplices = self._triangulate()
        self._edges = self._prepare_edges()

    def _triangulate(self) -> IntArray:
        """Triangulate the fixed point basis into deterministic simplices."""
        coordinates = self._basis.coordinates

        sort_keys = tuple(
            coordinates[:, dimension]
            for dimension in range(coordinates.shape[1] - 1, -1, -1)
        )
        order = np.lexsort(sort_keys)

        try:
            triangulation = sp.spatial.Delaunay(coordinates[order])
        except sp.spatial.QhullError as exc:
            raise ValueError(
                "EdgeIntersector requires a non-degenerate point basis that "
                "can be Delaunay triangulated."
            ) from exc
        except ValueError as exc:
            raise ValueError(
                "EdgeIntersector requires a non-degenerate point basis that "
                "can be Delaunay triangulated with at least two spatial "
                "dimensions."
            ) from exc

        return np.asarray(
            order[triangulation.simplices],
            dtype=np.int64,
        )

    def _prepare_edges(self) -> IntArray:
        """Extract unique edges and remove edges exceeding the validity range."""
        vertex_pairs = np.array(
            list(combinations(range(self._simplices.shape[1]), 2)),
            dtype=np.int64,
        )

        # Extract every edge from every simplex.
        edges = self._simplices[:, vertex_pairs].reshape(-1, 2)

        # Treat edges as undirected and remove duplicates shared by simplices.
        edges = np.sort(edges, axis=1)
        edges = np.unique(edges, axis=0)

        # Compute squared edge lengths in the PointSet coordinate space.
        edge_vectors = (
            self._basis.coordinates[edges[:, 0]] - self._basis.coordinates[edges[:, 1]]
        )

        # A common scale cancels from the RMS comparison and keeps squares
        # representable for very small or large edge lengths.
        length_scale = np.max(np.abs(edge_vectors))
        edge_lengths_squared = np.sum(
            (edge_vectors / length_scale) ** 2,
            axis=1,
        )

        # Compare squared lengths directly:
        #
        #     l_e^2 < validity_range^2 * mean(l^2)
        #
        # This is equivalent to limiting the edge length relative to the
        # RMS edge length without evaluating square roots.
        with np.errstate(over="ignore", under="ignore"):
            # Original point indices change edge order under permutations.
            # Sort lengths before reduction to keep a threshold tie
            # deterministic.
            mean_length_squared = np.mean(np.sort(edge_lengths_squared))

            threshold = (
                np.square(np.float64(self._validity_range)) * mean_length_squared
            )

        valid_mask = edge_lengths_squared < threshold
        edges = edges[valid_mask, :]

        if edges.shape[0] == 0:
            raise ValueError(
                "EdgeIntersector preparation removed all triangulation edges; "
                f"validity_range={self._validity_range} is too restrictive "
                "for this basis."
            )

        return edges

    def intersect(
        self,
        field: SampledField | ArrayLike,
        search_component: int,
        levels: ArrayLike,
    ) -> EdgeIntersections:
        """Intersect scalar field levels with the prepared triangulation edges.

        Array-like field input is normalized to a ``SampledField`` on the
        ``EdgeIntersector`` basis. An existing ``SampledField`` must use the
        same ordered ``PointSet`` basis.

        The selected field component is linearly interpolated along every
        prepared edge. For each requested level, an intersection is returned
        whenever the corresponding interpolation factor lies in ``[0, 1]``.

        The same interpolation factor is applied to the basis coordinates and
        to all field components.

        Parameters
        ----------
        field : SampledField or ArrayLike
            Field values to intersect. Array-like input is normalized to a
            ``SampledField`` on the ``EdgeIntersector`` basis.
        search_component : int
            Zero-based index of the scalar field component whose levels are
            intersected.
        levels : ArrayLike
            One scalar level or a one-dimensional collection of levels.

        Returns
        -------
        EdgeIntersections
            Coordinates, interpolated field values, requested levels, source
            edges, and interpolation factors for every valid edge-level
            intersection.

        Raises
        ------
        TypeError
            If ``search_component`` is not an integer or the supplied field
            cannot be converted by ``SampledField`` because of an invalid
            input type.
        ValueError
            If the supplied field cannot be normalized to a valid
            ``SampledField``, an existing ``SampledField`` uses a different
            ordered point basis, or ``levels`` is invalid, empty, or contains
            non-finite values.
        IndexError
            If ``search_component`` is outside the available field-component
            range.

        Notes
        -----
        Field values and prepared edges are aligned by their common ordered
        ``PointSet`` basis. Existing ``SampledField`` instances are checked
        against the basis fingerprint.

        Edges whose search-component endpoints contain ``NaN`` are ignored.

        Edges with equal search-component values are ignored because they do
        not define a unique linear intersection factor. If such an edge lies
        exactly on a requested level, its intersection is the complete edge
        rather than one unique point and is therefore not represented by this
        result type.

        Intersections at shared vertices are preserved per edge. No geometric
        deduplication is performed.

        The interpolated search-component value is not overwritten with the
        requested level. Tests can therefore verify the interpolation
        invariant without masking implementation errors.
        """
        try:
            sampled_field = (
                field
                if isinstance(field, SampledField)
                else SampledField(self._basis, field)
            )
        except (TypeError, ValueError) as exc:
            exc.add_note(
                "Error occurred while converting EdgeIntersector field input "
                "to SampledField."
            )
            raise

        if sampled_field.points.fingerprint != self._basis.fingerprint:
            raise ValueError(
                "Field rows must be defined on the same ordered support "
                "points referenced by the prepared EdgeIntersector edges; "
                f"field basis fingerprint="
                f"{sampled_field.points.fingerprint}, "
                f"basis fingerprint={self._basis.fingerprint}."
            )

        if isinstance(search_component, (bool, np.bool_)) or not isinstance(
            search_component,
            (int, np.integer),
        ):
            raise TypeError(
                f"search_component must be an integer; "
                f"got {type(search_component).__name__}."
            )

        if not (0 <= search_component < sampled_field.n_components):
            raise IndexError(
                f"search_component {search_component} is out of range; "
                f"expected 0 <= search_component < "
                f"{sampled_field.n_components}."
            )

        levels = np.asarray(levels, dtype=np.float64)

        if levels.ndim > 1:
            raise ValueError(
                "levels must be a scalar or one-dimensional array; "
                f"got shape={levels.shape}."
            )

        levels = np.atleast_1d(levels)

        if levels.size == 0:
            raise ValueError("levels must contain at least one level.")

        if not np.isfinite(levels).all():
            raise ValueError("levels must contain only finite values.")

        search_values = sampled_field.values[:, search_component]

        # Evaluate the selected field component at both vertices of every edge.
        start_values = search_values[self._edges[:, 0]]
        end_values = search_values[self._edges[:, 1]]

        with np.errstate(over="ignore", invalid="ignore"):
            denominator = end_values - start_values

        # A unique point intersection requires finite endpoint values and a
        # non-constant search component along the edge.
        finite_edges = np.isfinite(start_values) & np.isfinite(end_values)
        non_constant_edges = denominator != 0.0
        usable_edges = finite_edges & non_constant_edges

        # One interpolation factor is calculated for every edge-level pair:
        #
        #     lambda = (level - f_start) / (f_end - f_start)
        #
        # Shape: (E, L), with E prepared edges and L requested levels.
        lambdas = np.full(
            (self._edges.shape[0], levels.size),
            np.nan,
            dtype=np.float64,
        )

        with np.errstate(over="ignore", invalid="ignore"):
            np.divide(
                levels[None, :] - start_values[:, None],
                denominator[:, None],
                out=lambdas,
                where=usable_edges[:, None],
            )

        # Finite endpoints of opposite sign can have an overflowing difference.
        # Rescale only these rows; the interpolation factor is scale invariant.
        overflow_edges = finite_edges & ~np.isfinite(denominator)

        if np.any(overflow_edges):
            scale = np.maximum(
                np.abs(start_values[overflow_edges]),
                np.abs(end_values[overflow_edges]),
            )

            start_scaled = start_values[overflow_edges] / scale
            end_scaled = end_values[overflow_edges] / scale

            lambdas[overflow_edges] = (
                levels[None, :] / scale[:, None] - start_scaled[:, None]
            ) / (end_scaled - start_scaled)[:, None]

        # A level intersects the finite line segment exactly when lambda lies
        # within the closed interval [0, 1].
        valid_mask = (lambdas >= 0.0) & (lambdas <= 1.0)

        edge_indices, level_indices = np.nonzero(valid_mask)

        intersection_lambdas = lambdas[
            edge_indices,
            level_indices,
        ]
        intersection_edges = self._edges[edge_indices]
        intersection_levels = levels[level_indices]

        # Interpolate the geometric coordinates.
        start_coordinates = self._basis.coordinates[intersection_edges[:, 0]]
        end_coordinates = self._basis.coordinates[intersection_edges[:, 1]]

        weights = intersection_lambdas[:, None]

        coordinates = (1.0 - weights) * start_coordinates + weights * end_coordinates

        # Interpolate every field component using the same geometric factor.
        start_field_values = sampled_field.values[
            intersection_edges[:, 0],
            :,
        ]
        end_field_values = sampled_field.values[
            intersection_edges[:, 1],
            :,
        ]

        values = (1.0 - weights) * start_field_values + weights * end_field_values

        return EdgeIntersections(
            coordinates=coordinates,
            values=values,
            levels=intersection_levels,
            edges=intersection_edges,
            lambdas=intersection_lambdas,
        )


if __name__ == "__main__":
    """Visual smoke test using Himmelblau's function."""
    import matplotlib.pyplot as plt

    from vice_measeval.backend.numeric.benchmark.test_functions import (
        himmelblau,
    )

    x = np.linspace(-6, 6, 101)
    x_grid, y_grid = np.meshgrid(x, x)

    coordinates = np.column_stack(
        (
            x_grid.ravel(),
            y_grid.ravel(),
        )
    )

    field_values = himmelblau(
        x_grid.ravel(),
        y_grid.ravel(),
    )

    levels = np.array(
        [
            25.0,
            50.0,
            100.0,
        ]
    )

    intersector = EdgeIntersector(coordinates)

    result = intersector.intersect(
        field=field_values,
        search_component=0,
        levels=levels,
    )

    plt.figure()

    plt.contour(
        x_grid,
        y_grid,
        field_values.reshape(x_grid.shape),
        levels=levels,
    )

    for level in levels:
        mask = result.levels == level

        plt.scatter(
            result.coordinates[mask, 0],
            result.coordinates[mask, 1],
            s=10,
            label=f"Edge intersections: {level:g}",
        )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.legend()
    plt.show()
