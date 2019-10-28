# GPH
General purpose heuristics for mixed-integer programming

## Introduction
GPH implements general purpose heuristics for finding feasible solutions of good quality to problems of the form:

<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}%0A\min%20\quad%20%26%20c^Tx\\%0A\textrm{s.t.}%20\quad%20%26%20v%20\leq%20Ax%20\leq%20w%20\\%0A%26%20l%20\leq%20x%20\leq%20u%20\\%0A%20%20%26x_i%20\in%20\mathbb{Z}, \forall%20i%20\in%20I%20\quad%20%20%20%20%20\\%0A\end{aligned}">

The algorithms implemented are the heuristics usually embedded in MIP solvers.
A description of these algorithms can be found [here](https://opus4.kobv.de/opus4-zib/files/1112/Achterberg_Constraint_Integer_Programming.pdf).

The code is written in C++17 with a simple object-oriented design and provides data structures that allow for an efficient implementation of the algorithms.
Additionally, the code is parallelized using the Thread Building Blocks library.

These general purpose heuristics can be combined with problem-specific ones, that exploit the problem's structure.
This approach can be used to increase the likelihood of finding good quality solutions for hard problems.

## Design
* The problem is stored by the class `MIP`, with the constraint matrix stored in row-major and column-major order as two sparse matrices.
* A heuristic is a class derived from `HeuristicMethod` and implements the method `run`.
* The class `Search` takes a list of heuristics in the constructor and runs them in parallel.

## Compilation
GPH depends on an external LP solver and the Thread Building Blocks library and uses cmake for compilation. Three LP solvers are supported: Cplex, SoPlex and GLPK. 

Example: compiling with SoPlex

```
mkdir build
cd build
cmake .. -DSOLVER=soplex -DTBB_DIR=/path/to/tbb/cmake
make
```

If the LP solver is not installed system-wide, the path needs to be provided to cmake (by adding `-DSOPLEX_DIR=/path/to/soplex` in the example).

### Usage
The executable reads plain text files in MPS format and writes the best solution to disk if at least one was found.

```
SYNOPSIS
        ./gph <input file> [-l <tlimit>] [-t <nthreads>]

OPTIONS
        <tlimit>    time limit
        <nthreads>  number of threads to use
```

## Computational tests
TODO
