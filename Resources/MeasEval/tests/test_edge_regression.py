"""Contracts of the two pre-existing APIs, frozen before adding constraints."""
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from vice_measeval.backend.numeric.geometry.levelset.edge_intersector import (
    EdgeIntersector, EdgeIntersections, batched_edge_intersections,
)
from vice_measeval.backend.numeric.model.point_set import PointSet
from vice_measeval.backend.numeric.model.sampled_field import SampledField

TRIANGLE = np.array([[0., 0.], [1., 0.], [0., 1.]])


def test_edge_order_duplicates_and_constant_edges():
    geometry = EdgeIntersector(TRIANGLE)
    result = geometry.intersect(TRIANGLE[:, 0], 0, [0, .25, 1, .25])
    assert_array_equal(result.edges, [[0, 1]] * 4 + [[1, 2]] * 4)
    assert_array_equal(result.levels, [0, .25, 1, .25] * 2)
    assert_allclose(result.lambdas, [0, .25, 1, .25, 1, .75, 0, .75])
    assert_allclose(result.coordinates, [[0, 0], [.25, 0], [1, 0], [.25, 0],
                                         [0, 1], [.25, .75], [1, 0], [.25, .75]])
    assert_allclose(result.values[:, 0], result.coordinates[:, 0])


def test_interpolation_uses_same_weights_and_keeps_extra_nan():
    result = EdgeIntersector(TRIANGLE).intersect([[0, np.nan], [1, 7], [0, 9]], 0, [0, .5])
    assert_array_equal(result.edges, [[0, 1], [0, 1], [1, 2], [1, 2]])
    assert np.isnan(result.values[:2, 1]).all()  # Existing 0*NaN endpoint behavior.
    assert_allclose(result.values[2:, 1], [9, 8])


def test_overflow_rescaling_in_class_only():
    result = EdgeIntersector(TRIANGLE).intersect([-1e308, 1e308, -1e308], 0, [0])
    assert_allclose(result.lambdas, [.5, .5])
    assert_allclose(result.values, 0)


def test_nan_search_edges_and_empty_result():
    result = EdgeIntersector(TRIANGLE).intersect([np.nan, 1, 0], 0, [.5])
    assert_array_equal(result.edges, [[1, 2]])
    empty = EdgeIntersector(TRIANGLE).intersect([0, 1, 0], 0, [5])
    assert empty.coordinates.shape == (0, 2)
    assert empty.values.shape == (0, 1)


def test_ordered_basis_and_copy_contract():
    original = TRIANGLE.copy()
    basis = PointSet(original)
    field = SampledField(basis, original[:, 0])
    original[:] = 99
    assert_array_equal(basis.coordinates, TRIANGLE)
    with pytest.raises(ValueError, match="same ordered"):
        EdgeIntersector(TRIANGLE[::-1]).intersect(field, 0, [.5])
    with pytest.raises(ValueError, match="unique"):
        EdgeIntersector([[0, 0], [1, 0], [0, 0]])
    with pytest.raises(ValueError, match="Inf"):
        SampledField(basis, [0, np.inf, 1])
    result = EdgeIntersector(basis).intersect(field, 0, [.5])
    for array in [result.coordinates, result.values, result.levels, result.edges, result.lambdas]:
        with pytest.raises(ValueError):
            array.setflags(write=True)
    values = np.ones((1, 2))
    owned = EdgeIntersections(values, values, [0], [[0, 1]], [.5])
    values[:] = 8
    assert_array_equal(owned.coordinates, [[1, 1]])


@pytest.mark.parametrize("levels", [[], [[1]], [np.nan], [np.inf]])
def test_invalid_levels(levels):
    with pytest.raises(ValueError):
        EdgeIntersector(TRIANGLE).intersect([0, 1, 0], 0, levels)


def test_geometric_filter_determinism_under_permutation():
    points = np.array([[0, 0], [.1, 0], [0, .1], [5, 0], [5.1, 0], [5, .1]])
    first = EdgeIntersector(points, validity_range=1)
    permutation = [4, 2, 0, 5, 1, 3]
    second = EdgeIntersector(points[permutation], validity_range=1)
    def edges(geometry):
        return {tuple(sorted(map(tuple, geometry._basis.coordinates[e]))) for e in geometry._edges}
    assert edges(first) == edges(second)
    assert all(abs(points[a, 0] - points[b, 0]) < 1 for a, b in first._edges)


def test_matlab_level_major_supplied_order_duplicates_and_indices():
    samples = np.column_stack((TRIANGLE, TRIANGLE[:, 0]))
    edges = [[2, 1], [0, 1], [0, 1], [0, 2]]
    result = batched_edge_intersections(samples, 2, [0, .25, 1, .25], edges,
                                       index_base=0, numeric_policy="matlab-ieee")
    assert_array_equal(result["edge_indices"], [0, 1, 2] * 4)
    assert_array_equal(result["level_indices"], np.repeat(np.arange(4), 3))
    assert_array_equal(result["level_offsets"], [0, 3, 6, 9, 12])
    assert_allclose(result["lambdas"], [0, 0, 0, .25, .25, .25, 1, 1, 1, .25, .25, .25])
    assert_array_equal(result["samples"][:, 2], np.repeat([0, .25, 1, .25], 3))
    one = batched_edge_intersections(samples, 3, [[.25]], np.array(edges) + 1,
                                    index_base=1, numeric_policy="matlab-ieee")
    assert_array_equal(one["edge_indices"], [1, 2, 3])
    assert_array_equal(one["level_indices"], [1, 1, 1])
    assert_array_equal(one["level_offsets"], [0, 3])


def test_matlab_ieee_and_single_edge_quirk():
    samples = [[0, 0, 0], [1, 0, np.inf], [0, 1, 1]]
    result = batched_edge_intersections(samples, 2, [0], [[0, 1], [0, 2]],
                                       index_base=0, numeric_policy="matlab-ieee")
    assert_array_equal(result["lambdas"], [0, 0])
    assert_array_equal(result["samples"], [[0, 0, 0], [0, 0, 0]])
    with pytest.raises(ValueError, match="single input edge"):
        batched_edge_intersections(samples, 2, [0], [[0, 1]], index_base=0, numeric_policy="matlab-ieee")
    empty = batched_edge_intersections(samples, 2, [], [[0, 1]], index_base=0, numeric_policy="matlab-ieee")
    assert empty["samples"].shape == (0, 3)
    assert_array_equal(empty["level_offsets"], [0])
