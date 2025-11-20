
import numpy as np
from earthbox.utilities.mathtools import linear_interp_bisection


def run_test():
    """ Test the linear interpolation with bisection method.
    """
    gridx = np.linspace(0, 500_000.0, 100, dtype=np.float64)
    nx = len(gridx)
    gridy = np.linspace(3, 10.0, nx, dtype=np.float64)

    # print min and max y
    print(f'min y = {gridy.min()}')
    print(f'max y = {gridy.max()}')

    x = 800_000.0
    y = linear_interp_bisection(gridx, gridy, x)
    print(f'x {x} -> y = {y}')


if __name__ == '__main__':
    run_test()