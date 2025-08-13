// GPX.h - The header file for Generalized Partition Crossover (GPX) algorithm

#ifndef GPX_H
#define GPX_H

#define POPULATION_SIZE 50
#define MAX_GENERATIONS 1000
#define THRESHOLD_FITNESS 0.01


// Assuming you have a population array declared globally or passed as a parameter to the functions
int population[POPULATION_SIZE][MAX_NODES];

double evaluateFitness(struct Graph *graph, int *tour) {
    double tourLength = calculateTourLength(graph, tour);
    // If you have additional constraints or objectives, you can incorporate them here.

    return tourLength;
}

void gpcxCrossover(struct Graph *graph, int *parent1, int *parent2, int *offspring1, int *offspring2) {
    int numNodes = graph->numNodes;

    // Step 1: Initialize helper arrays and variables
    bool marked[numNodes];
    for (int i = 0; i < numNodes; i++) {
        marked[i] = false;
        offspring1[i] = -1;
        offspring2[i] = -1;
    }

    // Step 2: Choose the starting point randomly
    int startNode = rand() % numNodes;
    offspring1[0] = startNode;
    offspring2[0] = startNode;
    marked[startNode] = true;

    // Step 3: Assign nodes to partitions in alternating order
    int currNode1 = startNode;
    int currNode2 = startNode;
    bool isParent1Turn = true;

    for (int i = 1; i < numNodes; i++) {
        int nextNode;

        if (isParent1Turn) {
            // Find the unmarked node closest to currNode1 in parent2
            int closestNode = -1;
            double minDistance = -1;
            for (int j = 0; j < numNodes; j++) {
                if (!marked[j]) {
                    double distance = calculateDistance(graph->nodes[parent1[currNode1]], graph->nodes[parent2[j]]);
                    if (closestNode == -1 || distance < minDistance) {
                        closestNode = j;
                        minDistance = distance;
                    }
                }
            }
            nextNode = closestNode;
            marked[nextNode] = true;
            currNode1 = nextNode;
        } else {
            // Find the unmarked node closest to currNode2 in parent1
            int closestNode = -1;
            double minDistance = -1;
            for (int j = 0; j < numNodes; j++) {
                if (!marked[j]) {
                    double distance = calculateDistance(graph->nodes[parent2[currNode2]], graph->nodes[parent1[j]]);
                    if (closestNode == -1 || distance < minDistance) {
                        closestNode = j;
                        minDistance = distance;
                    }
                }
            }
            nextNode = closestNode;
            marked[nextNode] = true;
            currNode2 = nextNode;
        }

        // Assign the selected node to the offspring
        if (isParent1Turn) {
            offspring1[i] = nextNode;
        } else {
            offspring2[i] = nextNode;
        }

        // Switch the turn between parent1 and parent2
        isParent1Turn = !isParent1Turn;
    }
}
// Fisher-Yates shuffle για τυχαία αναδιάταξη
void shuffle(int *array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

// Αρχικοποίηση πληθυσμού με τυχαίες διαδρομές
void initializePopulation(struct Graph *graph) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        for (int j = 0; j < graph->numNodes; j++) {
            population[i][j] = j;
        }
        shuffle(population[i], graph->numNodes);
    }
}

/*
void gpcxAlgorithm(struct Graph *graph, int *tour) {
    int numNodes = graph->numNodes;
    // Αρχικοποίηση πληθυσμού
    initializePopulation(graph);
    for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Evaluate the fitness of the current tour (calculate total tour length)
        double tourLength = calculateTourLength(graph, tour);

        // Selection: Choose two parent tours from the population based on their fitness
        int parent1[MAX_NODES];
        int parent2[MAX_NODES];
        // For simplicity, you can randomly select two parents
        int parent1Idx = rand() % POPULATION_SIZE;
        int parent2Idx;
        do {
            parent2Idx = rand() % POPULATION_SIZE;
        } while (parent2Idx == parent1Idx);

        // Copy the selected parents to parent1 and parent2 arrays
        for (int i = 0; i < numNodes; i++) {
            parent1[i] = population[parent1Idx][i];
            parent2[i] = population[parent2Idx][i];
        }

        // Crossover: Create two offspring tours by applying GPCX on parent1 and parent2
        int offspring1[MAX_NODES];
        int offspring2[MAX_NODES];
        gpcxCrossover(graph, parent1, parent2, offspring1, offspring2);

        // Evaluate the fitness of the offspring tours
        double offspring1Length = evaluateFitness(graph, offspring1);
        double offspring2Length = evaluateFitness(graph, offspring2);

        // Replacement: Replace the least fit tours in the population with the offspring tours
        int minFit1Idx = -1, minFit2Idx = -1;
        double minFit1 = DBL_MAX, minFit2 = DBL_MAX;

        for (int i = 0; i < POPULATION_SIZE; i++) {
            double fitness = evaluateFitness(graph, population[i]);
            if (fitness < minFit1) {
                minFit2 = minFit1;
                minFit2Idx = minFit1Idx;
                minFit1 = fitness;
                minFit1Idx = i;
            } else if (fitness < minFit2) {
                minFit2 = fitness;
                minFit2Idx = i;
            }
        }

        // Replace the two least fit tours with the offspring tours
        for (int i = 0; i < numNodes; i++) {
            population[minFit1Idx][i] = offspring1[i];
            population[minFit2Idx][i] = offspring2[i];
        }

        // Termination condition (e.g., if a satisfactory solution is found)
        // For simplicity, we terminate if the best fitness reaches a threshold value
        double bestFitness = fmin(minFit1, minFit2);
        if (bestFitness < THRESHOLD_FITNESS) {
            break;
        }
    }

}*/
void gpcxAlgorithm(struct Graph *graph, int *tour) {
    int numNodes = graph->numNodes;
    initializePopulation(graph);

    int bestTour[MAX_NODES];
    double bestFitness = DBL_MAX;

    for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Evaluate current population and update best tour
        for (int i = 0; i < POPULATION_SIZE; i++) {
            double fitness = evaluateFitness(graph, population[i]);
            if (fitness < bestFitness) {
                bestFitness = fitness;
                for (int j = 0; j < numNodes; j++) {
                    bestTour[j] = population[i][j];
                }
            }
        }

        // Selection: Choose two parents randomly
        int parent1[MAX_NODES];
        int parent2[MAX_NODES];
        int parent1Idx = rand() % POPULATION_SIZE;
        int parent2Idx;
        do {
            parent2Idx = rand() % POPULATION_SIZE;
        } while (parent2Idx == parent1Idx);

        for (int i = 0; i < numNodes; i++) {
            parent1[i] = population[parent1Idx][i];
            parent2[i] = population[parent2Idx][i];
        }

        // Crossover
        int offspring1[MAX_NODES];
        int offspring2[MAX_NODES];
        gpcxCrossover(graph, parent1, parent2, offspring1, offspring2);

        // Evaluate offspring
        double offspring1Length = evaluateFitness(graph, offspring1);
        double offspring2Length = evaluateFitness(graph, offspring2);

        // Replacement: Find two worst individuals
        int worst1Idx = -1, worst2Idx = -1;
        double worst1Fitness = -1.0, worst2Fitness = -1.0;

        for (int i = 0; i < POPULATION_SIZE; i++) {
            double fitness = evaluateFitness(graph, population[i]);
            if (fitness > worst1Fitness) {
                worst2Fitness = worst1Fitness;
                worst2Idx = worst1Idx;
                worst1Fitness = fitness;
                worst1Idx = i;
            } else if (fitness > worst2Fitness) {
                worst2Fitness = fitness;
                worst2Idx = i;
            }
        }

        // Replace worst individuals with offspring
        for (int i = 0; i < numNodes; i++) {
            population[worst1Idx][i] = offspring1[i];
            population[worst2Idx][i] = offspring2[i];
        }

        // Optional: print progress
        if (generation % 100 == 0) {
            printf("Generation %d: Best fitness = %.4lf\n", generation, bestFitness);
        }

        // Termination condition
        if (bestFitness < THRESHOLD_FITNESS) {
            break;
        }
    }

    // Copy best tour to output
    for (int i = 0; i < numNodes; i++) {
        tour[i] = bestTour[i];
    }
}


#endif

