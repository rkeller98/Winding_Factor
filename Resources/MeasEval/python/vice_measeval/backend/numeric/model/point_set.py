import numpy as np
from numpy.typing import NDArray, ArrayLike
import hashlib

FloatArray = NDArray[np.float64]


class PointSet:
    """Immutable collection of unique points in a D-dimensional coordinate space.

    A ``PointSet`` stores ``N`` points row-wise in a two-dimensional coordinate
    array with shape ``(N, D)``, where ``N`` is the number of points and ``D``
    is the spatial dimension.

    The class represents geometry only. Associated field values, measurement
    validity, units, grid semantics, and algorithm-specific requirements are
    outside its responsibility.

    Parameters
    ----------
    coordinates : ArrayLike
        Point coordinates convertible to ``float64`` with shape ``(N,)`` or
        ``(N, D)``.
        Each row represents one point and each column one coordinate component.

        A one-dimensional input with shape ``(N,)`` is interpreted as N points
        in one spatial dimension and normalized internally to ``(N, 1)``.
        A single point in ``D`` dimensions must use shape ``(1, D)``.

    Raises
    ------
    ValueError
        If ``coordinates`` is neither one- nor two-dimensional, contains no points, has no
        spatial dimensions, contains non-finite values, or contains duplicate
        points.

    Notes
    -----
    Coordinates are copied during construction, converted to ``float64``, and
    made read-only. Consequently, a successfully constructed ``PointSet`` owns
    immutable coordinate data.

    Point uniqueness is exact. Numerically close but non-identical points are
    considered distinct; tolerance-based duplicate handling belongs to
    preprocessing or the consuming numerical algorithm.

    Signed zero is normalized before fingerprinting. The unchanged BLAKE2b
    fingerprint covers canonical shape, dtype, and ordered coordinate bytes;
    point order is part of basis identity.
    """

    def __init__(self, coordinates: ArrayLike) -> None:
        self._initialize_coordinates(coordinates)
        self._fingerprint = self._calculate_fingerprint()

    def _initialize_coordinates(self, coordinates: ArrayLike) -> None:
        """Normalize, validate, and store the point coordinates."""
        coordinates = np.array(coordinates, dtype=np.float64, copy=True)

        if coordinates.ndim == 1:
            coordinates = coordinates.reshape(-1, 1).copy()

        if coordinates.ndim != 2:
            raise ValueError(
                f"PointSet requires shape (N,) or (N, D), where N is the number "
                f"of points and D is the spatial dimension; got ndim={coordinates.ndim} "
                f"with shape={coordinates.shape}. "
                f"Shape (N,) means N points in one dimension. "
                f"For one point in D dimensions, use shape (1, D)."
            )

        if coordinates.shape[0] < 1:
            raise ValueError("PointSet requires at least one point; got N=0.")

        if coordinates.shape[1] < 1:
            raise ValueError(
                "PointSet requires at least one spatial dimension; got D=0."
            )

        if not np.isfinite(coordinates).all():
            raise ValueError(
                "PointSet coordinates must be finite; NaN, +Inf, and -Inf are not allowed."
            )

        if np.unique(coordinates, axis=0).shape[0] != coordinates.shape[0]:
            raise ValueError(
                "PointSet requires unique points; duplicate coordinates were found."
            )

        coordinates[coordinates == 0.0] = 0.0
        self._coordinates = coordinates
        self._coordinates.flags.writeable = False

    def _calculate_fingerprint(self) -> str:
        """Calculate a stable fingerprint of the normalized coordinate data."""
        hasher = hashlib.blake2b(digest_size=16)

        hasher.update(np.asarray(self._coordinates.shape, dtype=np.int64).tobytes())
        hasher.update(self._coordinates.dtype.str.encode())
        hasher.update(self._coordinates.tobytes(order="C"))

        return hasher.hexdigest()

    @staticmethod
    def _validate_index(index, count: int) -> None:
        """Validate a zero-based index against a fixed element count."""
        if isinstance(index, (bool, np.bool_)) or not isinstance(
            index, (int, np.integer)
        ):
            raise TypeError(f"Index must be an integer; got {type(index).__name__}.")

        if not (0 <= index < count):
            raise IndexError(
                f"Index {index} is out of range; expected 0 <= index < {count}."
            )

    def get_point(self, index) -> FloatArray:
        """Return one point from the point set.

        Parameters
        ----------
        index : int
            Zero-based index of the requested point. NumPy integer types are
            accepted as well.

        Returns
        -------
        FloatArray
            Read-only point coordinates with shape ``(D,)``.

        Raises
        ------
        TypeError
            If ``index`` is not an integer or is a boolean.
        IndexError
            If ``index`` is outside ``[0, N)``.
        """
        self._validate_index(index, self.n_points)
        return self._coordinates[index, :]

    def get_component(self, index) -> FloatArray:
        """Return one coordinate component for all points.

        Parameters
        ----------
        index : int
            Zero-based index of the requested coordinate component. NumPy
            integer types are accepted as well.

        Returns
        -------
        FloatArray
            Read-only coordinate values with shape ``(N,)``.

        Raises
        ------
        TypeError
            If ``index`` is not an integer or is a boolean.
        IndexError
            If ``index`` is outside ``[0, D)``.
        """
        self._validate_index(index, self.dimension)
        return self._coordinates[:, index]

    @property
    def coordinates(self) -> FloatArray:
        """Return the complete read-only coordinate matrix with shape ``(N, D)``."""
        coordinates = self._coordinates.view()
        coordinates.flags.writeable = False
        return coordinates

    @property
    def shape(self) -> tuple[int, int]:
        """Return the coordinate matrix shape ``(N, D)``."""
        return self._coordinates.shape[0], self._coordinates.shape[1]

    @property
    def n_points(self) -> int:
        """Return the number of points ``N``."""
        return self.shape[0]

    @property
    def dimension(self) -> int:
        """Return the spatial dimension ``D`` of each point."""
        return self.shape[1]

    @property
    def fingerprint(self) -> str:
        """Return the stable fingerprint of this point set."""
        return self._fingerprint
