# Performance analysis on Traveling Salesman Problem solving techniques

The Travelling Salesman Problem (TSP) is a fundamental and widely studied optimization problem in the field of Operations Research. It was originally formulated to find the shortest possible route that visits a set of cities exactly once and returns to the starting point. Despite its simple formulation, the TSP is NP-hard, making it computationally challenging to solve for large instances.

This project implements and experimentally evaluates a selection of more accessible algorithmic strategies—including classical heuristics, metaheuristics, and exact methods via CPLEX, providing insights into the trade-offs between computational efficiency and solution quality. For a full description of the methods and experimental results, see [the report](thesis/main.pdf).

## Prerequisites

- A working installation of IBM CPLEX Optimizer
- Environment variable `CPLEX_HOME` set to the root directory of your CPLEX installation

## Compilation

```bash
# build using the provided Makefile in project root
make
```

## Execution

```bash
./tsp_solver [OPTIONS]
```

**Options:**

- `-s, --seed <int>`  
  Random seed

- `-r, --random`  
  Generate random nodes

- `-n, --nnodes <int>`  
  Number of nodes (required with `-r`)

- `-f, --file <filename>`  
  Input file to load instance

- `-m, --method <string>`  
  Method(s) to solve TSP (you can combine multiple methods separated by `+`). Available methods:

  - `n_n`: Nearest Neighbour heuristic
  - `two_opt`: 2-opt local improvement
  - `random`: Random path generation
  - `vns`: Variable Neighbourhood Search (starts with Nearest Neighbour, then VNS)
  - `tabu`: Tabu Search (starts with Nearest Neighbour, then Tabu)
  - `benders`: Benders Decomposition via CPLEX
  - `branch_and_cut`: Branch-and-Cut via CPLEX
  - `hard_fix`: Hard Fixing via CPLEX
  - `local_branch`: Local Branching via CPLEX

- `-p, --plot`  
  Show and save the solution plot

- `-tl, --time_limit <float>`  
  Time limit in seconds

- `-rs, --runs <int>`  
  Number of runs (default: 1)

- `-h, --help`  
  Show this help message

**Methods parameters:**

- **Greedy**  
  `-gs, --greedy_starts <int>`  
  Number of greedy starts (default: 10)

- **VNS**  
  `-vkmin, --vns_kmin <int>`  
  Minimum k for VNS (default: 1)  
  `-vkmax, --vns_kmax <int>`  
  Maximum k (default: 5)  
  `-vlr, --vns_lr <float>`  
  Learning rate for VNS (default: 1.0)  
  `-vjps, --vns_jumps <int>`  
  Number of jumps for VNS (default: 0)

- **Tabu Search**  
  `-tten, --tabu-tenure <int>`  
  Initial tabu tenure size (default: (min_tenure + max_tenure) / 2)  
  `-ttenmin, --tabu-min <int>`  
  Minimum tabu tenure size (default: nnodes / 10)  
  `-ttenmax, --tabu-max <int>`  
  Maximum tabu tenure size (default: nnodes + nnodes / 10)  
  `-tnoimpr, --tabu-noimprove <int>`  
  No-improve limit for Tabu (default: 50)

- **Hard Fixing CPLEX**  
  `-hf_p, --hard_fixing_percentage <float>`  
  Hard fixing percentage (default: 0.3)  
  `-hf_lt, --hard_fixing_local_time <float>`  
  Hard fixing local time limit per iteration (default: 15)

- **Local Branching CPLEX**  
  `-lb_k, --local_branch_k <int>`  
  Local Branching initial K (default: 5)  
  `-lb_s, --local_branch_step <int>`  
  Local Branching K step (default: 5)

**Examples:**

```bash
# random 50 nodes, seed 123, nearest neighbour + 2-opt + VNS, 60s time limit
./tsp_solver -r -n 50 -s 123 -m n_n+two_opt+tabu -tl 60 --tabu-tenure 84 --tabu-min 6  --tabu-max 84  --tabu-noimprove 250 -gs 5 -rs 20

# load from file, random path, 20s time limit
./tsp_solver -f instance.txt -m random -tl 20
```
