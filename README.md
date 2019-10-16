# GPH
General purpose heuristics for mixed-integer programming

## Introduction
GPH implements general purpose heuristics for finding feasible solutions of good quality to problems of the form:

<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}%0A\min%20\quad%20%26%20c^Tx\\%0A\textrm{s.t.}%20\quad%20%26%20v%20\leq%20Ax%20\leq%20w%20\\%0A%26%20l%20\leq%20x%20\leq%20u%20\\%0A%20%20%26\forall%20i%20\in%20I%2C%20\quad%20x_i%20\in%20\mathbb{Z}%20%20%20%20\\%0A\end{aligned}">

The algorithms implemented are the heuristics usually embedded in branch-and-cut MIP solvers.
A detailed description of these algorithms can be found in the papers listed in the reference section.

The code is written in C++17 with a simple object-oriented design and provides data structures that allow for an efficient implementation of the algorithms.
Additionally, the code is parallelized using the Thread Building Blocks library.

These general purpose heuristics can be combined with problem-specific ones, that exploit the problem's structure.
This approach can be used to increase the likelihood of finding good quality solutions for hard problems.

## Design
* The problem is stored by the `MIP` class, with the constraint matrix stored in row-major and column-major order in packed storage.
* A heuristic is a class derived from `HeuristicMethod` and implements the method `run`.
* The `Search` class takes a list of heuristics in the constructor and runs them in parallel.

## Compilation
Some of the algorithms rely on a black-box LP solver, so one needs to be provided at compile time. Three LP solvers are supported: Cplex, SoPlex and GLPK.

### Requirements
* cmake
* gcc (>= 8) or clang (>= 5)
* TBB
* a supported LP solver

Example: compiling with SoPlex

```
mkdir build
cd build
cmake .. -DSOLVER=soplex -DTBB_DIR=/path/to/tbb/cmake
make
```

Additionally, if the LP solver is not installed system-wide, the path needs to be provided to cmake (by adding `-DSOPLEX_DIR=/path/to/soplex` in the example).

### Usage
The executable reads plain text files in `mps` format and writes the best solution, if any were found, in the `sol` format.

```
SYNOPSIS
        ./gph <input file> [-l <tlimit>] [-o <outfile>] [-v <verbosity>] [-t <nthreads>]

OPTIONS
        <tlimit>    time limit
        <outfile>   output file for the solution
        <nthreads>  number of threads to use
```

## Computational tests
TODO


## References
T. Achterberg. Constraint Integer Programming. PhD thesis, Technische Universität Berlin, 2007.

T. Achterberg and T. Berthold. Improving the Feasibility Pump. Technical Report 05-42, ZIB, 2005.

E. Balas, S. Ceria, M. Dawande, F. Margot, and G. Pataki. OCTANE: A New Heuristic for Pure 0-1 Programs.
