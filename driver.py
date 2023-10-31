import numpy as np
# from system import System as sys

def random_spin(s0, alpha=0.1):
    """Generate a new random spin based on the original one.

    Parameters
    ----------
    s0: np.ndarray

        The original spin that needs to be changed.

    alpha: float

        Larger alpha, larger the modification of the spin. Defaults to 0.1.

    Returns
    -------
    np.ndarray

        New updated spin, normalised to 1.

    """
    delta_s = (2 * np.random.random(3) - 1) * alpha
    s1 = s0 + delta_s
    return s1 / np.linalg.norm(s1)


class Driver:
    """Driver class.

    Driver class does not take any input parameters at initialisation.

    """

    def __init__(self):
        pass

    def drive(self, system, n, alpha=0.1):
        x=0
        while x < n:
            # step 1
            E0 = system.energy()

            # step 2
            i = np.random.randint(0, system.s.array.shape[0])
            j = np.random.randint(0, system.s.array.shape[1])

            # step 3
            finalSpin = None
            s0 = system.s.array[i,j,...]
            s1 = random_spin(s0, alpha=alpha)

            # step 4
            system.s.array[i,j] = s1
            E1= system.energy()

            # step 5
            deltaE = E1-E0
            if (deltaE < 0):
                # finalSpin= s1
                x+=1
            elif deltaE:
                finalSpin=s0
                system.s.array[i,j]= s0
            return 1





