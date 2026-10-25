# Reusable MeasEval geometry

Unmodified source files from [GUI_VICE-MeasurementEvalKit](https://github.com/Weg-Weiser/GUI_VICE-MeasurementEvalKit),
branch `modernization-baseline`, commit
[`8cf4dc2629c12e36d71b3784cf1ac1f5d4979fed`](https://github.com/Weg-Weiser/GUI_VICE-MeasurementEvalKit/commit/8cf4dc2629c12e36d71b3784cf1ac1f5d4979fed).
Paths below `python/` match upstream `SourceCode/python/`; paths below `matlab/`
match upstream `SourceCode/`. No GUI installation is needed.

## Python

Python 3.11 or newer, NumPy and SciPy:

```sh
python -m pip install numpy scipy
```

Run from the repository root:

```python
import sys
import numpy as np

sys.path.insert(0, "Resources/MeasEval/python")
from vice_measeval.backend.numeric.geometry.levelset.edge_intersector import (
    EdgeIntersector, batched_edge_intersections,
)

points = np.array([[0., 0.], [2., 0.], [0., 1.], [2., 2.]])
values = points[:, 0] + points[:, 1]
geometry = EdgeIntersector(points)
hits = geometry.intersect(values, search_component=0, levels=[0.5, 1.5])
print(hits.coordinates)

# Reuse the triangulation with another field or other levels.
other_hits = geometry.intersect(2 * values, search_component=0, levels=[1.])

# Use supplied edges and MATLAB arithmetic instead of rebuilding topology.
samples = np.column_stack((points, values))
result = batched_edge_intersections(
    samples, search_column=2, levels=[0.5, 1.5],
    edges=np.array([[0, 1], [0, 2], [1, 3], [2, 3], [1, 2]]),
    index_base=0, numeric_policy="matlab-ieee",
)
print(result["samples"])
```

`EdgeIntersector` prepares SciPy Delaunay topology once, removes duplicate edges
and filters long edges. `PointSet` and `SampledField` are included dependencies;
the original package imports are retained. Points must be finite and unique.
The module's optional visual demo also needs Matplotlib; its existing
`himmelblau` function is included under the original `benchmark/test_functions.py`
path. This is a mathematical example function, not a test suite.

The class and `batched_edge_intersections` have different contracts: the latter
preserves supplied edge order and repetitions, performs no triangulation or
length filtering, and fixes the searched column to the requested level. It can
also use one-based MATLAB indices (`index_base=1`). The upstream restriction
against a single supplied edge with nonempty levels is retained.

## MATLAB

```matlab
addpath('Resources/MeasEval/matlab');
points = [0 0; 2 0; 0 1; 2 2];
values = sum(points, 2);
triangles = delaunayn(points);
edges = backend.numeric.geometry.unique_edges(triangles, points);
hits = backend.numeric.geometry.edge_intersections([points values], 3, 0.5, edges);
```

Also included: `valid_edges` (simplex edge validity) and
`triangle_intersections` (plane intersections with interpolated extra columns).
`unique_edges` filters lengths over the supplied coordinate columns and unique
edges; `valid_edges` uses simplex edge occurrences. These are distinct upstream
policies. The research calculations and saved study results are unchanged.
