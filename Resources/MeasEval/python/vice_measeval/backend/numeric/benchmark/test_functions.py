import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def himmelblau(x: ArrayLike, y: ArrayLike) -> FloatArray:
    """Evaluate f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2.

    Inputs are converted to float64 and follow NumPy broadcasting. Array
    inputs produce the broadcast shape; scalar inputs produce a NumPy scalar.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    return (x**2 + y - 11) ** 2 + (x + y**2 - 7) ** 2
