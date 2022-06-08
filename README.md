# Wachspress2D3D

## 1) Quickstart

The purpose of this code is to define interpolant and l projector on 1D - 2D and 3D convex polytopes.
Start by running *run_me2.mlx*, or if it doesn't work because of compatibility issues, *run_me.m*.
For Matlab versions earlier than R2020B, you have to install MMX on your system :
<https://github.com/yuvaltassa/mmx>

![image](https://user-images.githubusercontent.com/72595712/172386691-7dddad3d-a374-476c-adab-eef12e27531a.png)

![image](https://user-images.githubusercontent.com/72595712/172387565-c65bc494-17af-4f78-9f18-595e04efa7c1.png)

## 2) Structure

* **run_me.m** and **run_me2.mlx** are examples files.
* **Domain.m** is a class file which includes the definition of polytopes, the projection and the display methods.
* **wachspress.m** is a function which returns the values of the Wachspress basis functions and their gradient a given points and a given domain. The points should be inside the domain.
* **t.m** and **mult.m** are multithread transpose and matrix multiplication operators to speed up the calculation in case there are lots of points.

## 3) License

Copyright (C) 2022 Théodore CHERRIERE (+mail si tu le souhaites)
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

[![GitHub license](https://img.shields.io/github/license/tcherrie/Wachspress2D3D)](https://github.com/tcherrie/Wachspress2D3D) [![GitHub release](https://img.shields.io/github/release/tcherrie/Wachspress2D3D.svg)](https://github.com/tcherrie/Wachspress2D3D/releases/) [![GitHub stars](https://img.shields.io/github/stars/tcherrie/Wachspress2D3D)](https://github.com/tcherrie/Wachspress2D3D/stargazers)
 
