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
static body *planets;
static const float DT = 0.1;
static double parallel_average_time = 0.0;
static double serial_average_time = 0.0;
void CalculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay) {
        float galacticPlaneRofX = b->x - a->x;
        float galacticPlaneRofY = b->y - a->y;
        float simulationSofteningLengthSquared = 10;
        float distanceSquared = sqrt(a->x * a->x + a->y * a->y) + simulationSofteningLengthSquared;
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
void GenerateDebugData (int numberOfPlanets) {
        const float accelerationScale = 100.0;
        srand(time(NULL));
        //printf("%d\n", rand());
        for (int i = 0; i < numberOfPlanets; i++) {
                float angle = ((float) i / numberOfPlanets) * 2.0 * PI + (((rand() % 10) / 10.0 - 0.5) * 0.5);
                float initialMass = 1.0;
                body planet = {.x = (rand() % 10) / 10.0, .y = (rand() % 10) / 10.0, .vx = cos(angle) * accelerationScale * ((rand() % 10) / 10.0), .vy = sin(angle) * accelerationScale * ((rand() % 10) / 10.0), .ax = 0.0, .ay = 0.0, .mass = initialMass * ((rand() % 10) / 10.0) + initialMass * 0.5};
                float scale = (planet.mass / (initialMass * 1.5)) + 0.1;
                planets[i] = planet;
        }
}
void SimulateWithBruteforce(int *argc, char ***argv, int items_per_process, int numberOfPlanets, float DT, float simulationTime) {
		MPI_Init(argc, argv);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) {
            parallel_average_time -= MPI_Wtime();
			printf("numberOfPlanets = %d\n", numberOfPlanets);
			printf("DT = %f\n", DT);
			for(int j = 0; j < numberOfPlanets; j++) {
				printf("planets #%d position of X = %f, of Y = %f\n", j + 1, planets[j].x,  planets[j].y);
				printf("planets #%d acceleration of X = %f, of Y = %f\n", j + 1, planets[j].ax,  planets[j].ay);
				printf("planets #%d velocity of X = %f, of Y = %f\n", j + 1, planets[j].vx,  planets[j].vy);
				printf("planets #%d mass = %f\n", j + 1, planets[j].mass);
			}
        }
        MPI_Datatype planetType;
        MPI_Type_contiguous(7, MPI_FLOAT, &planetType);
        MPI_Type_commit(&planetType);
        body *local_numbers = (body *) malloc(sizeof(*local_numbers) * items_per_process);
        MPI_Scatter(planets, items_per_process, planetType, local_numbers, items_per_process, planetType, 0, MPI_COMM_WORLD);
        
		
		for(double time = 0; time < 1; time += DT) {
			for(int i = 0; i < items_per_process; i++) {
					float total_ax = 0.0, total_ay = 0.0;
					for (size_t j = 0; j < numberOfPlanets; j++) {
							if (i + rank * items_per_process == j) {continue;}
							float ax = 0, ay = 0;
							CalculateNewtonGravityAcceleration(&local_numbers[i], &planets[j], &ax, &ay);
							total_ax += ax;
							total_ay += ay;
					}
					local_numbers[i].ax = total_ax;
					local_numbers[i].ay = total_ay;
					printf("Time = %f, planets #%d acceleration of X = %f, of Y = %f\n", time, i + rank * items_per_process + 1, local_numbers[i].ax, local_numbers[i].ay);
					integrate(&local_numbers[i], DT);
			}
		}
		
		body *planetsTemp = (body*) malloc(numberOfPlanets * sizeof(*planetsTemp));
        MPI_Gather(local_numbers, items_per_process, planetType, planetsTemp, items_per_process, planetType, 0, MPI_COMM_WORLD);
        if(rank == 0) {
                parallel_average_time += MPI_Wtime();
                printf("Parallel Process Time: %f\n", parallel_average_time);
        }
		if(planets != NULL) {
			free(planets);
		}
        free(local_numbers);
		if(planetsTemp != NULL) {
			free(planetsTemp);
		}
        MPI_Finalize();
}

int main(int argc, char **argv) {
        parallel_average_time = 0.0;
        //serial_average_time = 0.0;
        int items_per_process = 10;
		int numberOfPlanets = 100;
		float simulationTime = 1.0;
		planets = (body*) malloc(numberOfPlanets * sizeof(*planets));
		GenerateDebugData (numberOfPlanets);
        SimulateWithBruteforce(&argc, &argv, items_per_process, numberOfPlanets, DT, simulationTime);
        return 0;
}


























