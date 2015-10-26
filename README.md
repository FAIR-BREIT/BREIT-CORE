#BREIT: Balance Rate Equations for Ion Transportation

## Introduction

BREIT: Balance Equations for Atomic Reactions is an open source software written in C++, which solve analytically the equilibrium and non-equilibrium charge state distributions equations (c.f. [HANS-DIETER BETZ Rev. Mod. Phys. 44, 465](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.44.465)). 
#### Method
The differential equations (non-equilibrium case) are solved using the eigenvalues decomposition method, and the asymptotic limits (equilibrium case) are solved by matrix inversion. In addition, a Runge-Kutta method can be used for cross-check.
#### Input
BREIT needs electron-loss and -capture cross-sections (as well as initial conditions) as inputs in order to solve the (non-equilibrium) Betz equations.
Only charge q greater or equal than zero are supported. 
An [INI-like format](https://en.wikipedia.org/wiki/INI_file) is used as input file. BREIT input file examples can be found [here](https://github.com/NicolasWinckler/BREIT/blob/master/data/input/Example-8lvl-system-bis.txt) or [there](https://github.com/NicolasWinckler/BREIT/blob/master/data/input/Example-15lvl-system.txt). JSON or XML input files will also be supported. Supported cross-section units are : "cm2", "1e-16 cm2", or "xb", where x is a [SI prefix](https://en.wikipedia.org/wiki/Metric_prefix). The supported units of the penetration depth, required for the non-equilibrium solutions, are : "10 mg/cm2" or "xg/cm2" where x is a [SI prefix](https://en.wikipedia.org/wiki/Metric_prefix).
#### Output


Output equilibrium results are provided on demand in the form of :

* 1-electron approximation (txt file)
* no-approximation obtained from matrrix inversion (txt file)
* figures (pdf or root files)

Output non-equilibrium results are provided on demand in the form of :

* analytical formulae (txt file)
* table (txt file)
* figures (pdf or root files)

#### Console command line (for local installation)

* --input-file path/to/inputfile (required)
* --output-directory path/to/output/dir (optional)
* --plot (optional)
* --save (optional)
* --save-equilibrium (optional)
* --save-analytic (optional)
* --save-approximation (optional)
* --save-table (optional)
* --save-fig-ne (optional)
* --save-fig-e (optional)
* --save-root-ne (optional)
* --save-root-e (optional)

#### TODO

* distance-to-equilibrium computation
* server 
* numerical error from matrix operations + error propagation from input cross-sections
* multiple input files for ion-optics
* implement XML and JSON input-file parser
* refactor the alternative Runge-Kutta method
* more documentations

## Step by step local installation




#### 1. Build Requirements

cmake >= 2.8.2

C++11 Compiler:

* GCC >= 4.6
* clang >= 3.2

Third party libraries:

* boost>=1.56 (required)
* Lapack>=3.5 (required for non-equilibrium solutions)
    * on linux use e.g. sudo apt-get install liblapack-dev
* Root 5.3x or Root >= 6.04(required for non-equilibrium solutions)

#### 2. Install [BREIT](https://github.com/NicolasWinckler/BREIT)

    # "Breit_path_install" is used here as a directory name, but can be renamed as you want.
    # Set the shell variable SIMPATH to the boost installation directory
    export SIMPATH=~/boost_install
    [setenv SIMPATH ~/boost_install]

    cd ~/Breit_path_install
    git clone https://github.com/FAIR-BREIT/BREIT-CORE.git
    cd BREIT-CORE
    mkdir build
    cd build
    cmake ..
    make

## Licence 
BREIT is distributed under the terms of the GNU Lesser General Public Licence version 3 (LGPL) version 3.
