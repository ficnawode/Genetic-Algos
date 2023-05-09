#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// This used to be macros
int NUM_POINTS;
int POP_SIZE;
double MUTATION_PROBABILITY;
double CROSSOVER_PROBABILITY;
int NUM_ITERS;

#define IDX_POP(n, i) ((n)*NUM_POINTS + i)

typedef struct {
    double x, y;
} Point;

double calc_distance(const Point p1, const Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrtf(dx*dx + dy*dy);
}

double calc_path_length(const Point* pts, const int* idx) {
    double l = calc_distance(pts[idx[0]], pts[idx[NUM_POINTS-1]]);
    for (int i = 0; i < NUM_POINTS-1; ++i) {
		l += calc_distance(pts[idx[i]], pts[idx[i+1]]);
    }
    return l;
}

void swap(int* a, int* b) {
    const int c = *a;
    *a = *b;
    *b = c;
}

void shift(int* x, const int n) {
    int arr_shift[NUM_POINTS];
    for (int i = 0; i != NUM_POINTS; ++i)
		arr_shift[i] = x[i];
    for (int i = 0; i != NUM_POINTS; ++i)
		x[(i+n)%NUM_POINTS] = arr_shift[i];
}

// x - parents, y - childer:output
void crossover(const int* x1, const int* x2, int* y1, int* y2) {
    int a, b; // Crossover points
    a = rand() % NUM_POINTS;
    b = rand() % NUM_POINTS+1;
    if (b < a) swap(&a, &b);
    
    // Copy [a, b) from parents onto [0, b-a) of children
    for (int i = 0; i < b-a; ++i) {
		y1[i] = x2[i+a];
		y2[i] = x1[i+a];
    } 

    // Co tu sie kurwa dzieje
    int j1, j2, k;
    k = b;         // Point after b-a inteval
    j1 = j2 = b-a; // Where do we input next value
    for (int i = 0; i != NUM_POINTS; ++i) {
	    bool b1, b2;
	    b1 = b2 = true;
	for (int j = 0; j < b-a; ++j) {
	    if (x1[k%NUM_POINTS] == y1[j]) b1 = false;
	    if (x2[k%NUM_POINTS] == y2[j]) b2 = false;
	}
	if (b1) {
	    y1[j1] = x1[k%NUM_POINTS];
	    ++j1;
	}
	if (b2) {
	    y2[j2] = x2[k%NUM_POINTS];
	    ++j2;
	}
	++k;
    }
    shift(y1, a);
    shift(y2, a);
}

void mutate(int* in) {
    int a, b;
    do {
		a = rand() % NUM_POINTS;
		b = rand() % NUM_POINTS;
		if (b < a) swap(&a, &b);
    } while (a == b);
    while (a < b) {
		swap(&in[a], &in[b]);
		++a;
		--b;
    }
}

void init_population(int* pop) {
    for (int i = 0; i != POP_SIZE; ++i)
		for (int j = 0; j != NUM_POINTS; ++j)
	 		pop[IDX_POP(i, j)] = j;
    for (int i = 1; i < POP_SIZE; ++i) {
    	mutate(&pop[IDX_POP(i, 0)]);
    }
}

// TODO: There has got to be a better way to do this
void calc_fitness(const int* pop, const Point* points, double* f) {
    for (int i = 0; i < POP_SIZE; ++i) {
		// Longer the path, lower the fittness
		f[i] = 1.0 / calc_path_length(points, &pop[i*NUM_POINTS]);
    }
}

void roulette_selection(const int* pop, const double* f, int* selected_indexes) {
    double prob[POP_SIZE];
    double sum, prob_sum;
    sum = prob_sum = 0.0;
    for (int i = 0; i != POP_SIZE; ++i)
		sum += f[i];
    for (int i = 0; i != POP_SIZE; ++i) {
		prob[i] = prob_sum;
		prob_sum += f[i]/sum;
    }
    for (int i = 0; i != POP_SIZE; ++i) { // FIXME
		double num = (double) rand() / RAND_MAX;
		selected_indexes[i] = POP_SIZE-1;
		for (int j = 0; j < POP_SIZE-1; ++j) {
	    	if (num > prob[j] && num <= prob[j+1]) {
				// At j-th attempt of selection, i-th dude was selected
				selected_indexes[i] = j;
				break;
	    	}
		}
    }
}

void breed(const int* old_pop, int* new_pop, const int* selected_pairs) {
    for (int i = 0; i < POP_SIZE/2; ++i) {
		const int* x1 = &old_pop[IDX_POP(selected_pairs[2*i], 0)];
		const int* x2 = &old_pop[IDX_POP(selected_pairs[2*i+1], 0)];
		int* y1 = &new_pop[IDX_POP(2*i, 0)];
		int* y2 = &new_pop[IDX_POP(2*i+1, 0)];
		double r = (double) rand() / RAND_MAX;
		if (r <= CROSSOVER_PROBABILITY)
			crossover(x1, x2, y1, y2);
		else
	    	for (int i = 0; i != NUM_POINTS; ++i) {
				y1[i] = x1[i];
				y2[i] = x2[i];
	    	}
    }
}

int find_best_guy(const int* pop, const Point* pts) {
    int best = 0;
    double best_path = calc_path_length(pts, pop);
    for (int i = 1; i != POP_SIZE; ++i) {
		double curr_path = calc_path_length(pts, &pop[IDX_POP(i, 0)]);
	if (curr_path < best_path)
	    best = i;
    }
    return best;
}




int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, " Usage: %s <num_cities> <pop_size> <mut_prob> <cross_prob> <num_iters>\n", argv[0]);
		return -1;
	}
	NUM_POINTS = atoi(argv[1]);
	POP_SIZE = atoi(argv[2]);
	MUTATION_PROBABILITY = atof(argv[3]);
	CROSSOVER_PROBABILITY = atof(argv[4]);
	NUM_ITERS = atoi(argv[5]);

    Point points[NUM_POINTS];
    int population[NUM_POINTS*POP_SIZE]; // index permutations
    int tmp_population[NUM_POINTS*POP_SIZE];
    double f[POP_SIZE]; // Fittness
    int num_selections[POP_SIZE];
	int shortest_path_idxs[NUM_POINTS];
	double shortest_path = 1.0e10;
    srand(12);

    for (int i = 0; i != NUM_POINTS; ++i) {
		points[i].x = (double) rand() / RAND_MAX;
		points[i].y = (double) rand() / RAND_MAX;
    }
    init_population(population);
    
    int* curr_pop = population;
    int* new_pop = tmp_population;

	FILE* path_length_file = fopen("path_lengths.txt", "w");

    for (int i = 0; i != NUM_ITERS; ++i) {
		// gives f a fittness value for each guy
		calc_fitness(curr_pop, points, f);
		// num_selections - list of indexes of selected pop elements
		// 		    they can appear more than once
		roulette_selection(curr_pop, f, num_selections);
		breed(curr_pop, new_pop, num_selections);
		// Mutate every new guy with some probability
		for (int j = 0; j != POP_SIZE; ++j) {
		    const double r = (double) rand() / RAND_MAX;
	    	if (r <= MUTATION_PROBABILITY)
				mutate(&new_pop[j*NUM_POINTS]);
		}

		// Printing current shortest and all time shortest
		int curr_best_idx = find_best_guy(new_pop, points);
		double curr_best_length = calc_path_length(points, &new_pop[IDX_POP(curr_best_idx, 0)]);
		if (curr_best_length < shortest_path) {
			shortest_path = curr_best_length;
			for (int i = 0; i != NUM_POINTS; ++i)
				shortest_path_idxs[i] = new_pop[IDX_POP(curr_best_idx, i)];
		}
		fprintf(path_length_file, "%lf\t%lf\n", curr_best_length, shortest_path);
		int* tmp = curr_pop;
		curr_pop = new_pop;
		new_pop = tmp;
    }
    fclose(path_length_file);

    for (int i = 0; i != NUM_POINTS; ++i) {
		const Point p = points[shortest_path_idxs[i]];
		printf("%f\t%f\n", p.x, p.y);
    }

    return 0;
}
