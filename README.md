# GPH
General-purpose heuristics for mixed-integer programming

## Introduction
GPH implements approximation algorithms for quickly finding feasible solutions to problems of the form:

<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}%0A\min%20\quad%20%26%20c^Tx\\%0A\textrm{s.t.}%20\quad%20%26%20v%20\leq%20Ax%20\leq%20w%20\\%0A%26%20l%20\leq%20x%20\leq%20u%20\\%0A%20%20%26x_i%20\in%20\mathbb{Z}, \forall%20i%20\in%20I%20\quad%20%20%20%20%20\\%0A\end{aligned}">

The algorithms implemented are the heuristics usually embedded in MIP solvers.
A description of these algorithms can be found [here](https://opus4.kobv.de/opus4-zib/files/1112/Achterberg_Constraint_Integer_Programming.pdf).
The heuristics do not require the problem to have any particular structure.

The code is written in C++17 with a simple object-oriented design and provides data structures that allow for an efficient implementation of the algorithms.
Additionally, the code is parallelized using the Thread Building Blocks library. 

## Design
* The problem data is stored by the class `MIP`: the constraint matrix is stored in row-major and column-major order as two sparse matrices and the rest of the problem is stored as dense vectors.
* A feasibility heuristic is a class derived from `FeasibilityHeuristic` and implements the method `run`.
   Similarly, an improvement heuristic is a class derived from `ImprovementHeuristic` that implements the method `improve`.
* The class `Search` takes a list of heuristics as input and runs them in parallel. If the feasibility heuristics find a solution, it is passed to the improvement heuristics.

## Compilation
GPH depends on an external LP solver and the Thread Building Blocks library and uses CMake for compilation. Three LP solvers are supported: Cplex, SoPlex and GLPK.
zlib is an optional dependency to read input in gzip format.

Example: compiling with SoPlex

```
git clone https://github.com/anassmeskini/GPH
mkdir build
cd build
cmake .. -DSOLVER=soplex -DTBB_ROOT=/path/to/tbb
make
```

If the LP solver is not installed system-wide, the path needs to be provided to cmake (by adding `-DSOPLEX_DIR=/path/to/soplex` in the example).

### Usage
The executable reads plain text and gzip files in MPS format.

```
SYNOPSIS
        ./gph <input file> [-l <tlimit>] [-t <nthreads>] [-w] [-s <start_sol>] [-c <config>]

OPTIONS
        <tlimit>    time limit in seconds
        <nthreads>  number of threads to use
        -w          write solution to disk
        <start_sol> path to solution to improve
        <config>    configuration file
```

## Performance
Below is the result of running GPH on some instances of the set cover problem from [MIPLIB](https://miplib.zib.de) and [OR-LIBRARY](http://people.brunel.ac.uk/~mastjjb/jeb/info.html) with a time limit of three minutes on a standard desktop computer. The ability to find solutions and their quality is very sensitive to the type of problem.

| instance      | columns       | rows  | non-zero coefficients  | optimality gap  |
| ------------- |:-------------:| -----:|--------------:|-----:|
| core4872-1529 | 24k           | 5k    | 218k          | 29.3%|
| core4284-1064 |  21k          |  4k   | 245k          | 2.18%|
| core2586-950  | 13k           | 2K    | 104k          | 1.7% |
| core2536-691  | 15k           | 2K    | 177k          | 0.6% |
| n3seq24       | 119k          | 6k    | 3M            | 3.4% |
| scpk4         | 100K          | 2K    | 1M            | 62.9%|
| scpd4        | 4k           | 400    | 80k            | 14.6%|
| scpd5         | 4k           | 400    | 80k            | 9.1%|
| rail507      | 63K           | 507   | 409k          | 7.4% |
| rail582      | 55k           | 582   | 409k          | 3.9%|
| scpclr12     | 442           | 2k   | 126k           | 69.7%|
| scpclr13     | 715           | 4k   | 365k           | 102.8%|
