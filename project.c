#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#define PI 3.14159265358979323846
typedef struct _body {
        float x, y; // position
        float ax, ay; // acceleration
        float vx, vy; // velocity
        float mass; // mass
} body;
static bool IsItFirst = true;
static body *planets;
static body *planetsTemp;
static int items_per_process;
static int numberOfPlanets;
static const float DT = 0.1;
static double parallel_average_time = 0.0;
static double serial_average_time = 0.0;
void GenerateDebudData();
void SimulateWithBruteforce(int *argc, char ***argv);
void CalculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay) {
        float galacticPlaneRofX = b->x - a->x;
        float galacticPlaneRofY = b->y - a->y;
        float simulationSofteningLengthSquared = 10;
        float distanceSquared = (a->x * a->x + a->y * a->y) + simulationSofteningLengthSquared;
        float distanceSquaredCubed = distanceSquared * distanceSquared * distanceSquared;
        float inverse = 1.0 / sqrt(distanceSquaredCubed);
        float scale = b->mass * inverse;
        *ax = (galacticPlaneRofX * scale);
        *ay = (galacticPlaneRofY * scale);
}
void integrate(body *planet, float deltaTime) {
        planet->x += planet->vx * DT;
        planet->y += planet->vy * DT;
        planet->vx += planet->ax * DT;
        planet->vy += planet->ay * DT;
}
void GenerateDebugData () {
        const float accelerationScale = 100.0;
        srand(time(NULL));
        //printf("%d\n", rand());
        for (int i = 0; i < numberOfPlanets; i++) {
                float angle = ((float) i / numberOfPlanets) * 2.0 * PI + ((rand() % 1000 - 0.5) * 0.5);
                float initialMass = 5.0;
                body planet = {.x = rand() % 1000, .y = rand() % 1000, .vx = cos(angle) * accelerationScale * (rand() % 1000), .vy = sin(angle) * accelerationScale * (rand() % 1000), .mass = initialMass * (rand() % 1000) + initialMass * 0.5};
                float scale = (planet.mass / (initialMass * 1.5)) + 0.1;
                planets[i] = planet;
        }
}
void SimulateWithBruteforce(int *argc, char ***argv) {
        MPI_Init(argc, argv);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("Process %d out of %d\n", rank, world_size);
        if(rank == 0) {
                parallel_average_time -= MPI_Wtime();
                numberOfPlanets = world_size * items_per_process;
                printf("%d, %d, %d\n", numberOfPlanets, world_size, items_per_process);
                planets = (body*) malloc(numberOfPlanets * sizeof(*planets));
                planetsTemp = (body*) malloc(numberOfPlanets * sizeof(*planets));
                if(IsItFirst == true) {
                        IsItFirst = false;
                        GenerateDebugData();
                        planetsTemp = planets;
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Datatype planetType;
        MPI_Type_contiguous(7, MPI_FLOAT, &planetType);
        MPI_Type_commit(&planetType);
        body *local_numbers = (body *) malloc(sizeof(*local_numbers) * items_per_process);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(planets, items_per_process, planetType, local_numbers, items_per_process, planetType, 0, MPI_COMM_WORLD);
        printf("                Process %d out of %d\n", rank, world_size);
//      printf("AAAAAAAAAA                                                      %f\n", planets[0].x);
        for(int i = 0; i < items_per_process; i++) {
                float total_ax = 0, total_ay = 0;
                for (size_t j = 0; j < numberOfPlanets; j++) {
                        if (i == j) {continue;}
                        float ax = 0, ay = 0;
                        CalculateNewtonGravityAcceleration(&local_numbers[i], &planetsTemp[j], &ax, &ay);
                        printf("%f, %f\n", ax, ay);
                        total_ax += ax;
                        total_ay += ay;
                }
                local_numbers[i].ax = total_ax;
                local_numbers[i].ay = total_ay;
        }
        printf("Local number %f\n", local_numbers[0].ax);
        //printf("        Process %d out of %d\n", rank, world_size);
        MPI_Gather(local_numbers, items_per_process, planetType, planets, items_per_process, planetType, 0, MPI_COMM_WORLD);
        //printf("             Process %d out of %d\n", rank, world_size);
        if(rank == 0) {
                parallel_average_time += MPI_Wtime();
                printf("numberOfPlanets %d, DT %f\n", numberOfPlanets, DT);
                printf("Parallel Process Time: %f\n", parallel_average_time);
                for(int i = 0; i < numberOfPlanets; i++) {
                        integrate(&planets[i], DT);
                        printf("planets #%d acceleration of X = %f, of Y = %f\n", i + 1, planets[i].ax,  planets[i].ay);
                }
        }
        free(local_numbers);
        MPI_Finalize();
}

int main(int argc, char **argv) {
        parallel_average_time = 0.0;
        serial_average_time = 0.0;
        items_per_process = 10;

        //numberOfPlanets = world_size * items_per_process;
        //planets = (body*) malloc(numberOfPlanets * sizeof(*planets));
        //GenerateDebugData(numberOfPlanets);
        //printf("numberOfPlanets %d, DT %f\n", numberOfPlanets, DT);
        //for(double i = 0; i < 1; i += DT) {
                SimulateWithBruteforce(&argc, &argv);
        //      for(int j = 0; j < numberOfPlanets; j++) {
        //              printf("planets #%d acceleration of X = %f, of Y = %f\n", j + 1, planets[j].ax,  planets[j].ay);
        //      }
        //}

        return 0;
}
