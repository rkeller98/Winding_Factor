"""Run the implemented conventional contour estimator; other studies stay idle."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from python.k_factor.analyze_winding_factors import main as analyze_contours


def main(argv=None):
    """Forward optional input paths to the contour analysis."""
    return analyze_contours(argv)


if __name__ == "__main__":
    main()
