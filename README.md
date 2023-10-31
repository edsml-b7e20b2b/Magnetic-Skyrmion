# Search for a Magnetic Skyrmion

## Overview
Computational science is emerging as the third pillar of research and development in academia and industry across all science and engineering disciplines. Computational studies complement experimental and theoretical approaches and are sometimes the only feasible way to address research challenges, effective industrial design, and engineering of various products and systems
This project addresses the research challenge by showcasing a simulation environment for studying the behavior of magnetic skyrmions in a 2D spin-lattice. The simulation is based on a comprehensive implementation of energy terms, such as Zeeman energy, uniaxial anisotropy energy, exchange energy, and Dzyaloshinskii-Moriya energy. The Metropolis Monte Carlo algorithm is employed to minimise the system's energy and identify all equilibrium states.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Classes and Modules](#classes-and-modules)
  - [Spins Class](#spins-class)
    - [Spins Test Results](#spin-test-results)
  - [System Class](#system-class)
    - [System Test Results](#system-test-results)
  - [Driver Class](#driver-class)
    - [Performance Considerations](#performance-considerations)
- [Technical Details](#technical-details)
  - [Metropolis Monte Carlo Algorithm](#metropolis-monte-carlo-algorithm)
  - [Energy Computations](#energy-computations)
  - [Random Spin Generation](#random-spin-generation)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Prerequisites

Before running this simulation, please make sure you have the following prerequisites installed:

- Python (>= 3.6)
- NumPy (>= 1.18)
- Matplotlib (>= 3.2)

## Installation


- git clone https://github.com/edsml-b7e20b2b/magnetic-skyrmion-simulation.git
- cd magnetic-skyrmion-simulation
- pip install -r requirements.txt


## Usage

To perform a simulation, follow these steps:

1. Initialise a spin lattice using the Spins class.
2. Randomise the spins or set an initial spin configuration.
3. Visualise the spin lattice.
4. Initialise the Driver class.
5. Drive the system using the Metropolis Monte Carlo algorithm.


# Classes and Modules

## Spins Class

The Spins class represents a field of spins on a two-dimensional lattice. Each spin is a three-dimensional vector (sx, sy, sz). The class provides methods for initialisation, normalisation, and visualiation of the spin lattice.

The spins class consists of 4 methods; mean, normalise, randomise and plot.
## `mean(self)`

- **Description:** This method calculates the mean direction of all the spins in the lattice.
- **Details:**
    - I compute the sum of the x, y, and z components of the spins in the lattice.
    - Then, it calculates the mean by dividing the sum by the product of the lattice dimensions (nx * ny).
    - The result is returned as an array with three components (mx, my, mz).
```
    def mean(self):
        total_x = np.sum(self.array[..., 0])
        total_y = np.sum(self.array[..., 1])
        total_z = np.sum(self.array[..., 2])

        nx, ny = self.n
        mx = total_x / (nx*ny)
        my = total_y / (nx * ny)
        mz = total_z / (nx * ny)
        return np.array([mx,my,mz])
```


## `__abs__(self)`

- **Description:** This special method computes the absolute value (magnitude) of each spin in the lattice.
- **Details:**
    - I first square each component of the spins and sum them along the last axis ([-1]).
    - Then, I take the square root of the resulting sum, which gives the magnitude of each spin.
    - The result is returned as an array of magnitudes.
```
    def __abs__(self):
        sum = np.sum(self.array**2,axis=-1)
        norm = np.sqrt(sum)
        return norm[..., np.newaxis]
```

## `normalise(self)`

- **Description:** This method normalises the magnitude of all spins in the lattice to 1.
- **Details:**
    - It divides each spin vector by its magnitude (calculated using `__abs__`) in order to make sure all spins have unit length.
    - Note: This method fails the tests  if the `__abs__` method is not implemented.

## `randomise(self)`

- **Description:** This method initialises the lattice with random spins.
- **Details:**
    - It generates random spin vectors with components in the range of -1 to 1.
    - After the spins ahve been generated, they are normalised to ensure they all have unit length.

## `plot(self, **kwargs)`

- **Description:** This method visualises the spin lattice.
- **Details:**
    - It creates a quiver plot to represent the spins.
    - The size of the figure/frame is specified using the `figsize` argument.
    - Meshgrid creates a grid of coordinates that correspond to the positions of the spins. This is shown in the image below with a scale of 100.
    - Then x and y components from the spins are extracted (in variables u & v), and the quiver plot is created which shows arrows representing the spin directions.
    - The final visualisation displays the lattice structure and spin directions.
    - The plot is displayed with labels and a title.

    Note: Some commented-out code suggest adjustments that can be made to alter the plot's appearance, like changing the scale or colormap. For example, playing with the "scale" parameter and setting it to '100' reveals the spin/lattice formation.
  
  Additional information about quiverplots can be found here. https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.quiver.html

```
    def plot(self, **kwargs):
        figureSize = kwargs.get('figsize', (7,7))

        nx, ny, nz = self.array.shape
        x,y = np.meshgrid(range(nx), range(ny))

        u , v = self.array[:,:,0], self.array[:,:,1]
        plt.figure(figsize=figureSize)
        # After tinkering around, I discovered that at a scale of 100, we can see the uniform lattice arrangement of the
        # spins, but I have to scale it down in order to see the direction of the arrows.
        # plt.quiver(x, y, xx, yy, scale=100, cmap=cmap)

        plt.quiver(x, y, u, v, scale=10, cmap='Blues')

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Spin Lattice Visualisation')
        plt.show()
```
Scale at 10 | Scale at 100
--- | ---
![](Scale%20at%2010.png) | ![](Scale%20at%20100.png)


## Spins Test Results
![](Spin%20Tests.png)
All tests were successful

## Driver Class

The Driver class implements the Metropolis Monte Carlo algorithm whose main purpose is to minimise the system's energy and find all equilibrium states. It provides a method to drive the simulation.

#### Metropolis Monte Carlo Algorithm

The Metropolis Monte Carlo algorithm is used to minimise the system's energy. It randomly selects and modifies spins, this is then followed by acceptance or rejection based on energy changes, if the enerygy delta is higher than 0, we accept the changed spin, else we reject it and stay with the original spin. The process is repeated n times to reach equilibrium states.

#### Performance Considerations
The value of n in the Monte Carlo algorithm depends on the size of the system. Larger systems will require more iterations to reach equilibrium states. Hence optimisations may be necessary. The larger the system, the higher the runtime and test time. Testing my code took a long period because n = 10^6 for some test cases.


## System Class

The system's total energy is computed as the sum of individual energy terms, which include Zeeman energy, anisotropy energy, and exchange energy. These terms are evaluated for each spin configuration.
The systems class consists of these 5 functions:
#### `energy(self) -> float`

- Description: This function computes the total energy of the system as the sum of all individual energy terms which is returned from the energy() function.
- Returns: Total energy of the system.

#### `zeeman(self) -> float`

- Description: Calculates the Zeeman energy contribution.
- Returns: Total Zeeman energy.

#### `anisotropy(self) -> float`

- Description: Calculates the uniaxial anisotropy energy contribution.
- Returns: Total uniaxial anisotropy energy.

   - `self.u` is normalised to ensure that it's a unit vector, this sets the preferred spin orientation. Then`self.s.array` is dot-multiplied with the normalised `u` to align the spins with the anisotropy axis.
   - `sumOfNormalised` represents the aforementioned alignment of spins with the anisotropy axis.

   - `totalAnisotropyEnergy` is computed as the product of `-self.K` and the sum of aligned spins squared, with a negative sign indicating energy lowering. The result is returned.

```
 def anisotropy(self):
        # normalising u
        normalised = self.u / np.linalg.norm(self.u)
        sumOfNormalised = np.sum(self.s.array * normalised, axis=-1)

        totalAnisotropyEnergy = -self.K * np.sum(sumOfNormalised**2)
        # raise NotImplementedError
        return totalAnisotropyEnergy
```

#### `exchange(self) -> float`

- Description: Calculates the exchange energy contribution.
- Returns: Total exchange energy.

#### `dmi(self) -> float`

- Description: Computes the contribution of the Dzyaloshinskii-Moriya Interaction (DMI) energy in a magnetic system.
- Returns: Total Dzyaloshinskii-Moriya energy.
- 
- I iterate through all the spins in the lattice while considering their neighbours.
- The DMI energy between two spins s1 and s1 is edmi = −D · (s1 × s2), where D = Drij is the DMI vector. I implement this on line `if (i > nx-1) & (j > ny -2):
                    Edmi1 += np.sum((rX *self.s.array[i,j, 0] *self.s.array[i+rX,j,1]))`
  If there's a spin beside our current spin, I multiply it by the next spin. This gets the Edmi energy between 2 spins, I iterate through all the spins in the lattice structure, then get the dot product `totalEdmi += np.dot(Edmi1,Edmi2)` and multiply it by Dzyaloshinskii-Moriya energy constant.
- The rX and rY variables are adjusted based on the spin's position which influences subsequent calculations.
- In the nested for loop, I calculate the spin(s) energy using this formula :
` Edmi =-J [ Σ(i=0, to I =nx-1)(Σ(j=0, to j=ny-2)  rij*(Si,j* Si,j +1)  +  Σ(j=0, to j =ny-1)(Σ(i=0, to i=nx-2) rij*(Si,j* Si+1,j) ]`
The dot product of Edmi1 and Edmi2 is accumulated to obtain the total DMI energy (totalEdmi). The final DMI energy is ` -self.D *  totalEdmi` .
The calculated DMI energy is then returned.
  
```
    def dmi(self):
        nx, ny, _ = self.s.array.shape
        totalEdmi = 0
        rX, rY = 0, 0

        Edmi1, Edmi2 = 0, 0
        for i in range(nx):
            for j in range(ny):
                if i < nx - 1:
                    rX += 1
                if j < ny - 1:
                    rY += 1

                if (i > nx-1) & (j > ny -2):
                    Edmi1 += np.sum((rX *self.s.array[i,j, 0] *self.s.array[i+rX,j,1]))

                if (j > ny-1) & (i > nx-2):
                    Edmi2 += np.sum((rY * self.s.array[i,j, 0] * self.s.array[i+rY, j, 1] ))

                totalEdmi += np.dot(Edmi1,Edmi2)
        final= -self.D * totalEdmi

        return final
```
### System Test Results

![](Systems%20test.png)


## Acknowledgements
I want to thank Imperial College London ESE Department Professors and GTAs for their contributions to this project.
