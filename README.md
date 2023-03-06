# Wachspress2D3D


[![GitHub license](https://img.shields.io/github/license/tcherrie/Wachspress2D3D)](https://github.com/tcherrie/Wachspress2D3D) [![GitHub release](https://img.shields.io/github/release/tcherrie/Wachspress2D3D.svg)](https://github.com/tcherrie/Wachspress2D3D/releases/) [![GitHub stars](https://img.shields.io/github/stars/tcherrie/Wachspress2D3D)](https://github.com/tcherrie/Wachspress2D3D/stargazers) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6630214.svg)](https://doi.org/10.5281/zenodo.6630214)

The purpose of this code is to define interpolant and projector on 1D - 2D and 3D convex polytopes.

## Quickstart


Start by running `run_me2.mlx`, or if it doesn't work because of compatibility issues, `run_me.m`.
For Matlab versions earlier than R2020B, you have to install MMX on your system :
<https://github.com/yuvaltassa/mmx>.
Hereafter are some examples of projected points onto 2D and 3D polytopes obtained with this code :

![image](https://user-images.githubusercontent.com/72595712/172386691-7dddad3d-a374-476c-adab-eef12e27531a.png)

![image](https://user-images.githubusercontent.com/72595712/172387565-c65bc494-17af-4f78-9f18-595e04efa7c1.png)

## Contents



* `run_me.m` and **run_me2.mlx** are examples files.
* `Domain.m` is a class file which includes the definition of polytopes, the projection and the display methods. Among them, you can use :
    * `domain = Domain(type)` : constructor of a `Domain` instance ; `type` can be a string or a number, type `help Domain` for more details
    * `domain.projection(x)` : project the points `x` onto `domain`
    * `domain.plot()` : display the `domain`
    * `domain.plot_normalFan()` : display the normal fan of `domain`
    * `domain.plot_Wachspress(n)` : display the shape of the Wachspress basis function associated with the `n`-vertex
* `wachspress.m` is a function which returns the values of the Wachspress basis functions and their gradient a given points and a given domain. The points should be inside the domain. This file is an optimized version of a code provided by Floater *et. al.* [^1]
* `t.m` and `mult.m` are multithread transpose and matrix multiplication operators, respectively, to speed up the calculation when there are lots of points.

**NB** The code has been tested with success on Matlab R2020b and R2022b.

## Citation

Please use the following citation reference if you use the code:

    T. Cherrière and L. Laurent. tcherrie/wachspress2d3d: Wachspress2d3d (v1.0.0), June 2022. Zenodo. https://doi.org/10.5281/zenodo.6630215

Bibtex entry:

    @software{tcherrie_2022_6630215,
    author       = {Cherrière, Théodore and Laurent, Luc},
    title        = {tcherrie/Wachspress2D3D: Wachspress2D3D},
    month        = jun,
    year         = 2022,
    publisher    = {Zenodo},
    version      = {v1.0.0},
    doi          = {10.5281/zenodo.6630215},
    url          = {https://doi.org/10.5281/zenodo.6630215}
    }

**NB: version number and DOI must be adapted from [Zenodo's repository.](https://doi.org/10.5281/zenodo.6630214)** 

## License

Copyright (C) 2022 Théodore CHERRIERE (theodore.cherriere@ens-paris-saclay.fr)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

 
[^1]: M. Floater, A. Gillette, and N. Sukumar,
“Gradient bounds for Wachspress coordinates on polytopes,”
SIAM J. Numer. Anal., vol. 52, no. 1, pp. 515–532, 2014,
doi: 10.1137/130925712