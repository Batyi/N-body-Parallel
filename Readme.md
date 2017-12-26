N-body Parallel project with MPI
Parallelization of SimulateWithBruteforce algorithm

By default:
item_per_process = 10 (number of planets in each process)
numberOfPlanets = 100
simulationTime = 1.0
You can change in main function;

Compilation:
mpicc project.c -o project -lm  (write -lm because we use math functions like (cos, sin))
mpiexec -n 10 ./project
