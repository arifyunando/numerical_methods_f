# Numerical Methods for Engineers with Modern Fortran
> **Chapra S. C. & Canale R. P.** (2015). _Numerical methods for engineers (7th ed.)_. McGraw-Hill Higher Education.  

## Description
This repository contains implementations of the algorithm described in Chapra & Canale's book, "Numerical Methods for Engineers". Each folder contains different topics and hence can be studied separately.

### About `Numerical Methods`
> In numerical analysis, a numerical method is a mathematical tool designed to solve numerical problems. The implementation of a numerical method with an appropriate convergence check in a programming language is called a numerical algorithm. - [_Wikipedia_](https://en.wikipedia.org/wiki/Numerical_method#:~:text=In%20numerical%20analysis%2C%20a%20numerical,is%20called%20a%20numerical%20algorithm.)
> <br>---<br>
> Numerical methods are techniques by which mathematical problems are formulated so that they can be solved with arithmetic operations. Although there are many kinds of numerical methods, they have one common characteristic: they invariably involve large numbers of tedious arithmetic calculations. It is little wonder that with the development of fast, effi cient digital computers, the role of numerical methods in engineering problem solving has increased dramatically in recent years. - [_Chapra & Canale (2015)_](https://www.worldcat.org/nl/title/numerical-methods-for-engineers/oclc/897417371?referer=di&ht=edition)

## Compiling and Running the Codes
This repository is arranged such that every folder has its own topic and `CMakeLists.txt` file for easy compilation with [CMake](https://cmake.org/). To install CMake 3.15, please run this following commands
```sh
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt update
sudo apt upgrade
```

The code then can be compiled using GNU Fortran Compiler with this following commands.
```sh
sudo apt update 
sudo apt upgrade -y
sudo apt install -y build-essential 
cd <topic_folder>
mkdir build; cd build
cmake ..
make install
```

These commands will put the compiled binaries in a subfolder called bin which placed at `<project_parent>/<topic_folder>/build/` along with any input files that is needed for testing. 

For quick test, go into the binary folder and run `ctest` as follows. Note that `-VV` gives you very verbrose testing. Help information for the `ctest` command can be accessed through `--help` argument.

```sh
ctest -VV
```

Information on how to use the binary can be access through `--help, -H` argument when calling the binary, for example:

```sh
roots --help
```

## Table of Contents
1. [Roots of Equations](./1_Roots_of_Equations/)
2. [Linear Algebraic Equations](./2_Linear_Algebraic_Equations/)
3. [Optimizations](./3_Optimization/)
4. [Curve Fitting](./4_Curve_Fitting/)
5. [Numerical Differentiation and Integrations](./5_Numerical_Differentiation_and_Integration/)
6. [Ordinary Differential Equations](./6_Ordinary_Differential_Equation/)
7. [Partial Differential Equations](./7_Partial_Differential_Equation/)

## Disclaimers and Contribution
This is a simple study repository, thus the quality of the code is not the priority. Since it is basically a study log, contributions from third parties is not expected. Though, you are very welcome to fork this repository and give a GitHub star ⭐

## Author
**Arif Y. Sunanhadikusuma (Soen)** <br>
--- <br>
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://linkedin.com/in/arifyunando)
[![](https://img.shields.io/badge/Gmail-EA4335.svg?style=for-the-badge&logo=Gmail&logoColor=white)](mailto:arifyunando@gmail.com)<br>
_Civil Engineering (S.T.)_ <br>
Department of Civil Engineering <br>
Parahyangan Catholic University, Indonesia  <br> 
--- <br>
_Geotechnical Engineering (M.Sc)_ <br>
Civiel Techniek en Geowetenschappen (CiTG) <br>
TU Delft, The Netherlands