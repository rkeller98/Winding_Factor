"""Run the implemented Python k-factor study without starting placeholder studies."""

from python.k_factor.run_k_factor_analysis import main as run_k_factor_analysis


def main(argv=None):
    """Run conventional contours only; accept the same optional input paths."""
    return run_k_factor_analysis(argv)


if __name__ == "__main__":
    main()
