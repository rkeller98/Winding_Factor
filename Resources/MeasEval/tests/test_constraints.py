"""Analytic cuts and independent linear-programming references."""
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize import linprog

from vice_measeval.backend.numeric.geometry.levelset.edge_intersector import (
    EdgeIntersector, EdgeIntersections, EqualityConstraint as C, ConstraintIntersections,
)
from vice_measeval.backend.numeric.model.point_set import PointSet
from vice_measeval.backend.numeric.model.sampled_field import SampledField


def simplex(dim):
    return np.vstack((np.zeros(dim), np.eye(dim)))


def ordered(rows):
    return np.array(sorted(map(tuple, np.round(rows, 10))))


def test_triangle_interior_and_all_field_components():
    points = simplex(2)
    values = np.column_stack((points.sum(axis=1), 7 + 2*points[:, 0] - 3*points[:, 1]))
    g = EdgeIntersector(points)
    result = g.intersect_constraints(values, [C('field', 0), C('coordinate', 1)], [[.6, .2]])
    assert isinstance(result, ConstraintIntersections)
    assert result.query_status == ('intersected',)
    piece, = result.pieces
    assert piece.dimension == 0 and piece.rank == 2 and piece.affine_dimension == 0
    assert_allclose(piece.coordinates, [[.4, .2]], atol=1e-12)
    assert_allclose(piece.values, [[.6, 7.2]])
    assert_allclose(piece.weights @ points[piece.vertices], piece.coordinates)
    assert_allclose(piece.weights @ values[piece.vertices], piece.values)
    assert_allclose(piece.weights.sum(axis=1), 1)
    assert (piece.weights >= 0).all()
    assert piece.simplex_index == 0 and piece.query_index == 0


def test_tetrahedron_segment():
    result = EdgeIntersector(simplex(3)).intersect_constraints(
        None, [C('coordinate', 0), C('coordinate', 1)], [.2, .3])
    piece, = result.pieces
    assert piece.dimension == 1 and piece.affine_dimension == 1
    assert_allclose(ordered(piece.coordinates), [[.2, .3, 0], [.2, .3, .5]], atol=1e-12)
    assert piece.values.shape == (2, 0)


@pytest.mark.parametrize('m,dimension,count', [(2, 2, 3), (3, 1, 2), (4, 0, 1)])
def test_four_dimensional_slices(m, dimension, count):
    target = [.1, .2, .1, .2][:m]
    result = EdgeIntersector(simplex(4)).intersect_constraints(
        None, [C('coordinate', j) for j in range(m)], target)
    piece, = result.pieces
    assert piece.dimension == dimension and len(piece.coordinates) == count
    assert_allclose(piece.coordinates[:, :m], np.tile(target, (count, 1)), atol=1e-12)
    assert np.all(piece.coordinates.sum(axis=1) <= 1 + 1e-12)


def test_more_conditions_than_dimension_and_contradictions():
    points = simplex(4)
    conditions = [C('coordinate', j) for j in range(4)] + [C('field', 0)]
    g = EdgeIntersector(points)
    result = g.intersect_constraints(points.sum(axis=1), conditions,
                                    [[.1, .2, .1, .2, .6], [.1, .2, .1, .2, .8]])
    piece, = result.pieces
    assert piece.rank == 4 and piece.redundant_conditions == 1
    assert result.query_status == ('intersected', 'empty')
    assert result.inconsistent_counts[1] == 1


def test_explicit_target_tuples_not_cartesian_product():
    g = EdgeIntersector(simplex(2))
    result = g.intersect_constraints(None, [C('coordinate', 0), C('coordinate', 1)],
                                    [[.2, .3], [.4, .1], [.2, .3]])
    assert len(result.pieces) == 3
    assert [p.query_index for p in result.pieces] == [0, 1, 2]
    assert_allclose([p.coordinates[0] for p in result.pieces], [[.2, .3], [.4, .1], [.2, .3]])


def test_redundant_whole_edge_and_whole_simplex():
    points = simplex(2)
    g = EdgeIntersector(points)
    result = g.intersect_constraints(points[:, 0], [C('field', 0), C('coordinate', 0)], [0, 0])
    piece, = result.pieces
    assert piece.rank == 1 and piece.redundant_conditions == 1 and piece.dimension == 1
    assert_allclose(ordered(piece.coordinates), [[0, 0], [0, 1]], atol=1e-12)
    whole = g.intersect_constraints(np.ones((3, 2)), [C('field', 0), C('field', 1)], [1, 1])
    piece, = whole.pieces
    assert piece.rank == 0 and piece.dimension == 2
    assert_allclose(ordered(piece.coordinates), ordered(points))
    empty = g.intersect_constraints(points[:, 0], [C('field', 0), C('coordinate', 0)], [.2, .3])
    assert not empty.pieces and empty.inconsistent_counts[0] == 1


def test_affine_dimension_can_exceed_actual_boundary_dimension():
    points = simplex(2)
    g = EdgeIntersector(points)
    result = g.intersect_constraints(points.sum(axis=1), [C('field', 0), C('field', 0)], [0, 0])
    piece, = result.pieces
    assert piece.rank == 1 and piece.affine_dimension == 1 and piece.dimension == 0
    assert_allclose(piece.coordinates, [[0, 0]], atol=1e-12)


def test_nearly_dependent_and_numerically_redundant():
    points = simplex(2)
    g = EdgeIntersector(points)
    conditions = [C('field', 0), C('field', 1)]
    values = np.column_stack((points[:, 0], points[:, 0] + 1e-8 * points[:, 1]))
    result = g.intersect_constraints(values, conditions, [.2, .2 + 3e-9])
    piece, = result.pieces
    assert piece.rank == 2 and piece.near_singular
    assert_allclose(piece.coordinates, [[.2, .3]], atol=2e-8)
    values[:, 1] = points[:, 0] + 1e-14 * points[:, 1]
    result = g.intersect_constraints(values, conditions, [.2, .2])
    piece, = result.pieces
    assert piece.rank == 1 and piece.near_singular
    assert piece.scaled_residual <= 1e-10


def test_different_units_large_offsets_and_overflow():
    points = simplex(2)
    g = EdgeIntersector(points)
    values = np.column_stack((points[:, 0] * 1e-200, points[:, 1] * 1e200))
    result = g.intersect_constraints(values, [C('field', 0), C('field', 1)], [2e-201, 3e199])
    assert_allclose(result.pieces[0].coordinates, [[.2, .3]], atol=1e-12)
    values = np.column_stack((1e12 + 4*points[:, 0], 1e-12 * points[:, 1]))
    result = g.intersect_constraints(values, [C('field', 0), C('field', 1)], [1e12 + 1, .5e-12])
    assert_allclose(result.pieces[0].coordinates, [[.25, .5]], atol=1e-12)
    result = g.intersect_constraints([-1e308, 1e308, -1e308],
                                    [C('field', 0), C('coordinate', 1)], [0, .2])
    assert_allclose(result.pieces[0].coordinates, [[.5, .2]], atol=1e-12)


def test_shared_boundary_provenance_is_not_globally_deduplicated():
    points = np.array([[0, 0], [1, 0], [0, 1], [1, 1.]])
    g = EdgeIntersector(points)
    shared = sorted(set(g._simplices[0]) & set(g._simplices[1]))
    target = points[shared].mean(axis=0)
    result = g.intersect_constraints(None, [C('coordinate', 0), C('coordinate', 1)], target)
    assert len(result.pieces) == 2
    assert len({p.simplex_index for p in result.pieces}) == 2
    for piece in result.pieces:
        assert piece.dimension == 0
        assert_allclose(piece.coordinates, [target], atol=1e-12)


def test_nan_constraints_skip_support_but_extra_components_are_local():
    points = simplex(2)
    g = EdgeIntersector(points)
    result = g.intersect_constraints([np.nan, 1, 0], [C('field', 0), C('coordinate', 1)], [.2, .3])
    assert not result.pieces and result.query_status == ('unsupported',)
    assert_array_equal(result.nonfinite_simplex_indices, [0])
    values = np.column_stack((points[:, 0], [np.nan, 7, 9]))
    result = g.intersect_constraints(values, [C('field', 0), C('coordinate', 1)], [[.2, .3], [1, 0]])
    assert np.isnan(result.pieces[0].values[0, 1])
    assert_allclose(result.pieces[1].values[0, 1], 7)


def test_filtered_simplices_do_not_bridge_holes():
    points = np.array([[0, 0], [.1, 0], [0, .1], [5, 0], [5.1, 0], [5, .1]])
    g = EdgeIntersector(points, validity_range=1)
    result = g.intersect_constraints(None, [C('coordinate', 0), C('coordinate', 1)], [2.5, .025])
    assert not result.pieces
    assert len(result.rejected_simplex_indices) > 0


def test_single_field_delegates_without_topology_preparation(monkeypatch):
    g = EdgeIntersector(simplex(2))
    monkeypatch.setattr(np.linalg, 'svd', lambda *a, **k: pytest.fail('SVD in edge path'))
    ordinary = g.intersect([0, 1, 0], 0, [.25, 0, .25])
    result = g.intersect_constraints([0, 1, 0], [C('field', 0)], [[.25], [0], [.25]])
    assert isinstance(result, EdgeIntersections)
    for name in ['coordinates', 'values', 'levels', 'edges', 'lambdas']:
        assert_array_equal(getattr(result, name), getattr(ordinary, name))
    assert not hasattr(g, '_constraint_topology')


def test_single_field_simplex_opt_out_represents_constant_pieces():
    g = EdgeIntersector(simplex(2))
    result = g.intersect_constraints([1, 1, 1], [C('field', 0)], [1], edge_fast_path=False)
    assert result.pieces[0].dimension == 2


def test_reuse_geometry_without_stale_field_values(monkeypatch):
    points = simplex(2)
    g = EdgeIntersector(points)
    monkeypatch.setattr(g, '_triangulate', lambda: pytest.fail('retriangulated'))
    c = [C('field', 0), C('coordinate', 1)]
    a = g.intersect_constraints(points[:, 0], c, [.2, .3])
    topology = g._constraint_topology
    b = g.intersect_constraints(2*points[:, 0], c, [.2, .3])
    assert topology is g._constraint_topology
    assert_allclose(a.pieces[0].coordinates, [[.2, .3]])
    assert_allclose(b.pieces[0].coordinates, [[.1, .3]])


def test_immutable_results_and_ordered_field_basis():
    points = simplex(2)
    g = EdgeIntersector(points)
    targets = np.array([[.2, .3]])
    result = g.intersect_constraints(None, [C('coordinate', 0), C('coordinate', 1)], targets)
    targets[:] = 99
    assert_allclose(result.targets, [[.2, .3]])
    piece = result.pieces[0]
    for array in [result.targets, piece.weights, piece.coordinates, piece.values, piece.vertices,
                  result.rejected_simplex_indices, result.numerical_failure_counts]:
        with pytest.raises(ValueError):
            array.setflags(write=True)
    other = SampledField(PointSet(points[::-1]), [0, 1, 2])
    with pytest.raises(ValueError, match='ordered'):
        g.intersect_constraints(other, [C('field', 0), C('coordinate', 1)], [.2, .3])


@pytest.mark.parametrize('dim,m', [(2, 2), (3, 2), (4, 2), (4, 3)])
def test_vertices_match_independent_linear_programming_extrema(dim, m):
    rng = np.random.default_rng(381 + dim + m)
    points = simplex(dim)
    coefficients = rng.normal(size=(m, dim))
    values = points @ coefficients.T
    target = np.full(dim, .1) @ coefficients.T
    result = EdgeIntersector(points).intersect_constraints(values, [C('field', j) for j in range(m)], target)
    piece, = result.pieces
    equations = np.vstack((np.ones(len(points)), values.T))
    for direction in rng.normal(size=(12, dim)):
        reference = linprog(points @ direction, A_eq=equations, b_eq=np.r_[1, target],
                            bounds=(0, None), method='highs')
        assert reference.success
        assert_allclose(np.min(piece.coordinates @ direction), reference.fun, atol=1e-9)


def test_empty_requests_and_outside_ranges():
    g = EdgeIntersector(simplex(2))
    conditions = [C('coordinate', 0), C('coordinate', 1)]
    empty = g.intersect_constraints(None, conditions, np.empty((0, 2)))
    assert empty.pieces == () and empty.query_status == ()
    missing = g.intersect_constraints(None, conditions, [[2, .2], [.8, .8]])
    assert missing.query_status == ('empty', 'empty')
    assert missing.range_rejected_counts[0] == 1


@pytest.mark.parametrize('kwargs', [dict(chunk_size=0), dict(chunk_size=True), dict(rtol=-1),
                                  dict(atol=[0]), dict(rank_tol=0), dict(barycentric_tol=-1)])
def test_invalid_options(kwargs):
    with pytest.raises((TypeError, ValueError)):
        EdgeIntersector(simplex(2)).intersect_constraints(None, [C('coordinate', 0), C('coordinate', 1)], [.2, .3], **kwargs)


@pytest.mark.parametrize('constraints,targets', [([], []), ([C('coordinate', 8)], [.2]),
                                               ([C('coordinate', 0), C('coordinate', 1)], [.2]),
                                               ([C('field', 0)], [np.nan])])
def test_invalid_requests(constraints, targets):
    with pytest.raises((TypeError, ValueError, IndexError)):
        EdgeIntersector(simplex(2)).intersect_constraints(None, constraints, targets)
