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

void CalculateNewtonGravityAcceleration(body *a, body *b, float *ax, float *ay) {
	float galacticPlaneRofX = abs(b->x - a->x);
	float galacticPlaneRofY = abs(b->y - a->y);
	float simulationSofteningLengthSquared = 0.5;
	float distanceSquared = (a->x * a->x + a->y * a->y) + simulationSofteningLengthSquared;
	float distanceSquaredCubed = distanceSquared * distanceSquared * distanceSquared;
    float inverse = 1.0 / sqrt(distanceSquaredCubed);
    float scale = b->mass * inverse;
	*ax = (galacticPlaneRofX * scale);
	*ay = (galacticPlaneRofY * scale);
}
void integrate(body *planet, float deltaTime) {
	planet->x += planet->ax * DT * DT / 2.0 + planet->vx * DT;
	planet->y += planet->ay * DT * DT / 2.0 + planet->vy * DT;
	planet->vx += planet->ax * DT;
	planet->vy += planet->ay * DT;
}
void SimulateWithBruteforce(int numberOfPlanets) {
	for(size_t i = 0; i < numberOfPlanets; i++) {
		float total_ax = 0, total_ay = 0;
		for (size_t j = 0; j < numberOfPlanets; j++) {
			if (i == j) {continue;}
			float ax = 0, ay = 0;
			CalculateNewtonGravityAcceleration(&planets[i], &planets[j], &ax, &ay);
			//printf("%f, %f\n", ax, ay);
			total_ax += ax;
			total_ay += ay;
		}
		planets[i].ax = total_ax;
		planets[i].ay = total_ay;
		integrate(&planets[i], DT);
	}
}

void GenerateDebugData (int numberOfPlanets) {
	const float accelerationScale = 100.0;
	srand(time(NULL));
	for (int i = 0; i < numberOfPlanets; i++) {
		int randomValue = rand() % (1000) + 1;
		float angle = ((float) i / numberOfPlanets) * 2.0 * PI + ((randomValue - 0.5) * 0.5);
		float initialMass = 5.0;
		body planet = {.x = (rand() % (1000) + 1), .y = (rand() % (1000) + 1), /*.ax = , .ay = ,*/ .vx = cos(angle) * accelerationScale * (rand() % (1000) + 1), .vy = sin(angle) * accelerationScale * (rand() % (1000) + 1), .mass = initialMass * (rand() % (1000) + 1) + initialMass * 0.5};
        float scale = (planet.mass / (initialMass * 1.5)) + 0.1;
        planets[i] = planet;
    }
}

int main(int argc, char **argv) {
	int numberOfPlanets = 10;
	planets = (body*) malloc(numberOfPlanets * sizeof(*planets));
	GenerateDebugData(numberOfPlanets);
	printf("numberOfPlanets %d, DT %d\n", numberOfPlanets, DT);
	for(double i = 0; i < 10; i += DT) {
		SimulateWithBruteforce(numberOfPlanets);
		for(int j = 0; j < numberOfPlanets; j++) {
			printf("planets #%d acceleration of X = %f, of Y = %f\n", j + 1, planets[j].ax,  planets[j].ay);
		}
	}
	return 0;
}




























