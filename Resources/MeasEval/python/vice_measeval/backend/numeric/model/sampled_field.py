from vice_measeval.backend.numeric.model.point_set import PointSet
import numpy as np
from numpy.typing import NDArray, ArrayLike

FloatArray = NDArray[np.float64]


class SampledField:
    """Immutable values of one sampled field defined on a PointSet.

    A ``SampledField`` associates every point in a ``PointSet`` with one or
    more value components.

    Values are stored row-wise in a two-dimensional array with shape
    ``(N, K)``, where ``N`` is the number of sampled points and ``K`` is the
    number of value components per point.

    A one-dimensional input with shape ``(N,)`` is interpreted as a scalar-valued
    field and normalized internally to ``(N, 1)``. A field with
    multiple value components uses shape ``(N, K)`` with ``K > 1``.

    Parameters
    ----------
    points : PointSet
        Points at which the field values are sampled.
    values : ArrayLike
        Field values convertible to ``float64`` with shape ``(N,)`` or ``(N, K)``.
        The first dimension must match ``points.n_points``.

        ``NaN`` values are allowed to represent missing or invalid field
        values. Positive and negative infinity are not allowed.

    Raises
    ------
    TypeError
        If ``points`` is not a ``PointSet``.
    ValueError
        If ``values`` is neither one- nor two-dimensional, has no value components,
        contains a different number of samples than ``points``, or contains
        positive or negative infinity.

    Notes
    -----
    The supplied values are copied during construction, converted to
    ``float64``, and made read-only.

    The associated ``PointSet`` is referenced directly because ``PointSet``
    already guarantees immutable geometry.
    """

    def __init__(self, points: PointSet, values: ArrayLike) -> None:
        self._initialize(points, values)

    def _initialize(self, points: PointSet, values: ArrayLike) -> None:
        """Normalize, validate, and store the sampled field data."""
        if not isinstance(points, PointSet):
            raise TypeError(f"points must be a PointSet; got {type(points).__name__}.")
        values = np.array(values, dtype=np.float64, copy=True)

        if values.ndim == 1:
            values = values.reshape(-1, 1).copy()

        if values.ndim != 2:
            raise ValueError(
                f"SampledField requires values with shape (N,) or (N, K), "
                f"where N is the number of sampled points and K is the number of "
                f"value components per point; got ndim={values.ndim} with "
                f"shape={values.shape}. "
                f"Shape (N,) means a scalar-valued field."
            )

        if values.shape[1] < 1:
            raise ValueError(
                "SampledField requires at least one value component per point; got K=0."
            )

        if points.n_points != values.shape[0]:
            raise ValueError(
                f"PointSet and values must contain the same number of samples N; "
                f"got {points.n_points} points and {values.shape[0]} value rows."
            )

        if np.isinf(values).any():
            raise ValueError(
                "SampledField values must not contain +Inf or -Inf. "
                "NaN is allowed to represent missing or invalid field values."
            )

        self._points = points
        self._values = values
        self._values.flags.writeable = False

    @property
    def points(self) -> PointSet:
        """Return the PointSet on which the field is sampled."""
        return self._points

    @property
    def values(self) -> FloatArray:
        """Return the read-only field values with shape ``(N, K)``."""
        values = self._values.view()
        values.flags.writeable = False
        return values

    @property
    def n_points(self) -> int:
        """Return the number of sampled points ``N``."""
        return self._points.n_points

    @property
    def n_components(self) -> int:
        """Return the number of value components ``K`` per sampled point."""
        return self._values.shape[1]
