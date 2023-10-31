import numbers
import matplotlib.pyplot as plt
import numpy as np


class Spins:
    """Field of spins on a two-dimensional lattice.

    Each spin is a three-dimensional vector s = (sx, sy, sz). Underlying data
    stucture (``self.array``) is a numpy array (``np.ndarray``) with shape
    ``(nx, ny, 3)``, where ``nx`` and ``ny`` are the number of spins in the x
    and y directions, respectively, and 3 to hold all three vector components of
    the spin.

    Parameters
    ----------
    n: Iterable

        Dimensions of a two-dimensional lattice ``n = (nx, ny)``, where ``nx``
        and ``ny`` are the number of atoms in x and y directions, respectively.
        Values of ``nx`` and ``ny`` must be positive integers.

    value: Iterable

        The value ``(sx, sy, sz)`` that is used to initialise all spins in the
        lattice. All elements of ``value`` must be real numbers. Defaults to
        ``(0, 0, 1)``.

    """

    def __init__(self, n, value=(0, 0, 1)):
        # Checks on input parameters.
        if len(n) != 2:
            raise ValueError(f"Length of iterable n must be 2, not {len(n)=}.")
        if any(i <= 0 or not isinstance(i, int) for i in n):
            raise ValueError("Elements of n must be positive integers.")

        if len(value) != 3:
            raise ValueError(f"Length of iterable value must be 3, not {len(n)=}.")
        if any(not isinstance(i, numbers.Real) for i in n):
            raise ValueError("Elements of value must be real numbers.")

        self.n = n
        self.array = np.empty((*self.n, 3), dtype=np.float64)
        self.array[..., :] = value

        if not np.isclose(value[0] ** 2 + value[1] ** 2 + value[2] ** 2, 1) :
            # we ensure all spins' magnitudes are normalised to 1.
            self.normalise()

    @property
    def mean(self):
        total_x = np.sum(self.array[..., 0])
        total_y = np.sum(self.array[..., 1])
        total_z = np.sum(self.array[..., 2])

        nx, ny = self.n
        mx = total_x / (nx*ny)
        my = total_y / (nx * ny)
        mz = total_z / (nx * ny)
        return np.array([mx,my,mz])

    def __abs__(self):
        sum = np.sum(self.array**2,axis=-1)
        norm = np.sqrt(sum)
        return norm[..., np.newaxis]

    def normalise(self):
        """Normalise the magnitude of all spins to 1."""
        self.array = self.array / abs(self)

    def randomise(self):
        """Initialise the lattice with random spins.

        Components of each spin are between -1 and 1: -1 <= si <= 1, and all
        spins are normalised to 1.

        """
        self.array = 2 * np.random.random((*self.n, 3)) - 1
        self.normalise()

    def plot(self, **kwargs):
        figureSize = kwargs.get('figsize', (7,7))
        nx, ny, nz = self.array.shape
        x,y = np.meshgrid(range(nx), range(ny))

        u , v = self.array[:,:,0], self.array[:,:,1]
        plt.figure(figsize=figureSize)
        # After tinkering around, i discovered that at a scale of 100, we can see the uniform lattice arrangment of the
        # spins, but i have to scale it down in order to see the direction of the arrows.
        # plt.quiver(x, y, xx, yy, scale=100, cmap=cmap)

        plt.quiver(x, y, u, v, scale=10, cmap='viridis')

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Spin Lattice Visualisation')
        plt.show()




# x= np.arange(8,5,3)
# array = np.random.rand(4, 3, 3)
# print(array)
