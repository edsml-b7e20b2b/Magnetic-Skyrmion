import numpy as np


class System:
    """System object with the spin configuration and necessary parameters.

    Parameters
    ----------
    s: mcsim.Spins

        Two-dimensional spin field.

    B: Iterable

        External magnetic field, length 3.

    K: numbers.Real

        Uniaxial anisotropy constant.

    u: Iterable(float)

        Uniaxial anisotropy axis, length 3. If ``u`` is not normalised to 1, it
        will be normalised before the calculation of uniaxial anisotropy energy.

    J: numbers.Real

        Exchange energy constant.

    D: numbers.Real

        Dzyaloshinskii-Moriya energy constant.

    """

    def __init__(self, s, B, K, u, J, D):
        self.s = s
        self.J = J
        self.D = D
        self.B = B
        self.K = K
        self.u = u

    def energy(self):
        """Total energy of the system.

        The total energy of the system is computed as the sum of all individual
        energy terms.

        Returns
        -------
        float

            Total energy of the system.

        """
        return self.zeeman() + self.anisotropy() + self.exchange() + self.dmi()

    def zeeman(self):
        totalZeemanEnergy= -np.sum(self.s.array * self.B)
        # raise NotImplementedError
        return totalZeemanEnergy

    def anisotropy(self):
        # normalising u
        normalised = self.u / np.linalg.norm(self.u)
        sumOfNormalised = np.sum(self.s.array * normalised, axis=-1)

        totalAnisotropyEnergy = -self.K * np.sum(sumOfNormalised**2)
        # raise NotImplementedError
        return totalAnisotropyEnergy

    # Eex = -J[Σ(i=0, to i = nx - 1)(Σ(j=0, to j=ny-2) Si, j * Si, j + 1) + Σ(j=0, to j = ny - 1)(Σ(i=0, to i=nx-2) Si, j * Si + 1, j)]
    def exchange(self):

        Eex = 0.0
        nx, ny, _ = self.s.array.shape


        for i in range(nx):
            for j in range(ny):
                if i > nx-1:
                    Eex -= np.dot(self.s.array[i, j], self.s.array[i - 1, j])

                # calculating the contribution from the spin below (if exists)
                # if i is not the last spin,we're sure there are spins remaining, so we can proceed to increment the index
                if i < nx - 1:
                    Eex -= np.dot(self.s.array[i, j], self.s.array[i + 1, j])

                # calculating the contribution from the spin to the left if its there
                if j > ny-1:
                    Eex -= np.dot(self.s.array[i, j], self.s.array[i, j - 1])

                # calculating the contribution from the spin to the right if thers a spin there
                if j < ny - 1:
                    Eex -= np.dot(self.s.array[i, j], self.s.array[i, j + 1])

        # Here, i multiply by the exchange energy constant J, which is the final part of the formula
        final= (self.J * Eex)
        return final


    # please refer tp readme file for additional information
    def dmi(self):
        # assigning the x and y coordinates/dimensions from the spin tuple to nx&ny variables to use in the comditional
        # statement below
        nx, ny, _ = self.s.array.shape
        totalEdmi = 0
        rX, rY = 0, 0

        # Adding up the
        Edmi1, Edmi2 = 0, 0
        for i in range(nx):
            for j in range(ny):
                if i < nx - 1:
                    rX += 1
                if j < ny - 1:
                    rY += 1

                if (i > nx-1) & (j > ny -2):
                    Edmi1 += np.sum((rX *self.s.array[i,j, 0] *self.s.array[i+rX, j,1]))

                if (j > ny-1) & (i > nx-2):
                    Edmi2 += np.sum((rY * self.s.array[i,j, 0] * self.s.array[i+rY, j, 1] ))

                totalEdmi += np.dot(Edmi1,Edmi2)
        final= -self.D * totalEdmi

        return final
