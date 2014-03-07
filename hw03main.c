/* 
 * File:   hw03main.c
 * Author: jbao
 *
 * Created on February 16, 2014, 8:02 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double random_number();
void spin_random(int);

void print_intset(int size, int* set) {
    int i;
    for (i = 0; i < size; i++) {
        int entry = set[i];
        printf("%d: %d\n", 1 + i, entry);
    }
}

int* generate_lattice(int n, int m) {
    int size = n*m;
    int* lattice = malloc(size * sizeof (int));
    int i;
    for (i = 0; i < size; i++) { // randomly assign value of -1 or +1
        lattice[i] = (int) (random_number()*2)*2 - 1;
    }
    return lattice;
}

void h03test01() {
    int n = 6;
    int m = 6;
    int size = n*m;
    int* lattice = generate_lattice(n, m);
    print_intset(size, lattice);
    free(lattice);
}

int** get_connectivity(int n, int m) { // n = rows, m = cols
    int i, j;
    int** connectivity = malloc((n * m + 1) * sizeof (void*));
    connectivity[0] = malloc(sizeof (int)); // number of subarrays
    connectivity[0][0] = n*m;
    // index ordering scheme is #ofneighbors, then from west, north, east, south
    for (i = 1; i < n - 1; i++) { // top row, no corners
        int k = 3;
        int* toprow = malloc((1 + k) * sizeof (int));
        toprow[0] = k;
        toprow[1] = 1 + i - 1;
        toprow[2] = 1 + i + 1;
        toprow[3] = 1 + i + m;
        connectivity[1 + i] = toprow;
    }
    for (j = 1; j < m - 1; j++) { // body rows, no edges
        for (i = 1; i < n - 1; i++) {
            int id = n * j + i + 1;
            int k = 4;
            int* bodyrow = malloc((1 + k) * sizeof (int));
            bodyrow[0] = k;
            bodyrow[1] = id - 1;
            bodyrow[2] = id - m;
            bodyrow[3] = id + 1;
            bodyrow[4] = id + m;
            connectivity[1 + i + j * m] = bodyrow;
        }
    }
    for (i = (n - 1) * m; i < n * m - 1; i++) { // bottom row, no corners
        int k = 3;
        int* botrow = malloc((1 + k) * sizeof (int));
        botrow[0] = k;
        botrow[1] = 1 + i - 1;
        botrow[2] = 1 + i - m;
        botrow[3] = 1 + i + 1;
        connectivity[1 + i] = botrow;
    }
    for (i = 2; i < n; i++) { // edges
        int id = 1 + n * (i - 1);
        int k = 3;
        int* leftrow = malloc((k + 1) * sizeof (int));
        leftrow[0] = k;
        leftrow[1] = id - m;
        leftrow[2] = id + 1;
        leftrow[3] = id + m;
        connectivity[(i - 1) * m + 1] = leftrow;
        int* rightrow = malloc((k + 1) * sizeof (int));
        id = 1 + n * i - 1;
        rightrow[0] = k;
        rightrow[1] = id - 1;
        rightrow[2] = id - m;
        rightrow[3] = id + m;
        connectivity[i * m ] = rightrow;
    }

    int k = 2;
    int* topleft_corner = malloc((1 + k) * sizeof (int)); // 0
    topleft_corner[0] = k;
    topleft_corner[1] = 1 + 1;
    topleft_corner[2] = 1 + m;
    connectivity[1] = topleft_corner;
    int* topright_corner = malloc((1 + k) * sizeof (int)); // m
    topright_corner[0] = k;
    topright_corner[1] = 1 + m - 2;
    topright_corner[2] = 1 + m + m - 1;
    connectivity[m] = topright_corner;
    int* bottomleft_corner = malloc((1 + k) * sizeof (int)); // (m-1)*n
    bottomleft_corner[0] = k;
    bottomleft_corner[1] = 1 + (m - 2) * n;
    bottomleft_corner[2] = 1 + (m - 1) * n + 1;
    connectivity[1 + (m - 1) * n] = bottomleft_corner;
    int* bottomright_corner = malloc((1 + k) * sizeof (int)); // m*n-1
    bottomright_corner[0] = k;
    bottomright_corner[1] = 1 + m * n - 2;
    bottomright_corner[2] = 1 + (m - 1) * n - 1;
    connectivity[m * n] = bottomright_corner;

    return connectivity;
}

/**
 * array must have the first subarray be an array of length 1 with int entry equal to the number of content subarrays.
 * each subarray must start with the first value equal to the number of content entries in the subarray.
 * @param array
 */
void print_intint(int** array) {
    int i, j;
    int n = array[0][0];
    for (i = 1; i <= n; i++) {
        if (array[i] == NULL) {
            printf("error: missing entry at index %d\n", i);
        } else {
            printf("%d: ", i);
            int m = array[i][0];
            for (j = 1; j <= m; j++) {
                int entry = array[i][j];
                printf("%d ", entry);
            }
            printf("\n");
        }
    }
}

void h03test02() {
    int n = 6;
    int m = 6;
    int** connections = calloc(n*m, sizeof (void*));
    int i;
    connections[0] = malloc(sizeof (int));
    connections[0][0] = 5;
    for (i = 0; i < n; i++) {
        int k = 5;
        connections[i] = calloc(k + 1, sizeof (int));
        connections[i][0] = k;
    }
    connections[1][3] = 4;
    print_intint(connections);
    free(connections);
}

void h03test03() {
    int n = 6;
    int m = 6;
    int** connections = get_connectivity(n, m);
    print_intint(connections);
    free(connections);
}

int potential_e(int location, int* lattice, int** connectivity) {
    int num_neighbors = connectivity[location][0];
    int i;
    int neighbor_sum = 0;
    int own_value = lattice[location - 1];
    for (i = 1; i <= num_neighbors; i++) {
        int neighbor_id = connectivity[location][i];
        int neighbor_value = lattice[neighbor_id - 1];
        neighbor_sum += neighbor_value;
    }
    int result = -own_value * neighbor_sum;
    return result;

}

void h03test04() {
    int location = 2;
    int n = 6;
    int m = 6;
    int* lattice = generate_lattice(n, m);
    int** connections = get_connectivity(n, m);
    int u = potential_e(location, lattice, connections);

    printf("local potential at position %d: %d\n", location, u);
    print_intint(connections);
    print_intset(n*m, lattice);

    free(lattice);
    free(connections);
}

int potential(int* lattice, int** connectivity) {
    int size = connectivity[0][0];
    int location;
    int result = 0;
    for (location = 1; location <= size; location++) {
        int neighbor_sum = 0;
        int own_value = lattice[location - 1];
        int num_neighbors = connectivity[location][0];
        int j;
        for (j = 1; j <= num_neighbors; j++) {
            int neighbor_id = connectivity[location][j];
            int neighbor_value = lattice[neighbor_id - 1];
            neighbor_sum += neighbor_value;
        }
        int local_contribution = -own_value * neighbor_sum;
        result += local_contribution;
    }
    return result;
}

void h03test05() {
    int n = 6;
    int m = 6;
    int* lattice = generate_lattice(n, m);
    int** connections = get_connectivity(n, m);
    int u = potential(lattice, connections);

    printf("lattice potential: %d\n", u);
    print_intint(connections);
    print_intset(n*m, lattice);

    free(lattice);
    free(connections);
}

static double T = 1;

int flipstate(int potential_difference) {
    if (potential_difference < 0) {
        return -1; // flip
    } else {
        double probability = exp((double) (-potential_difference) / T);
        double random = random_number();
        if (random > probability) {
            return 1;
        } else {
            return -1;
        }
    }
}

void h03test06() {
    int flip;
    int i = 0;
    int times = 100;
    for (i = 0; i < times; i++) {
        flip = flipstate(1);
        printf("flip: %d\n", flip);
    }
}

double magnetization(int size, int* lattice) {
    int i;
    int sum = 0;
    for (i = 0; i < size; i++) {
        sum += lattice[i];
    }
    return sum / size;
}

void file_intout(int size, int* data, char* outpath) {
    FILE *f = fopen(outpath, "w");
    int i;
    for (i = 0; i < size; i++) {
        fprintf(f, "%d\n", data[i]);
    }
    fclose(f);
}

void montecarlo_v(int n, int m, int* lattice, int** connectivity, int num_steps) {
    int location, j;
    int size = n*m;
    double avg_potential = 0;
    //printf("mc inital potential: %d\n", potential(lattice, connectivity));
    for (j = 0; j < num_steps; j++) {
        for (location = 1; location <= size; location++) {
            int p1 = potential(lattice, connectivity);
            lattice[location - 1] *= -1;
            int p2 = potential(lattice, connectivity);
            int potential_difference = p2 - p1;
            int flip = flipstate(potential_difference);
            lattice[location - 1] *= -flip;
            avg_potential += potential(lattice, connectivity);
        }
    }

    avg_potential /= size * num_steps;
    double final_potential = (double) potential(lattice, connectivity);
    printf("mc avg potential: %3.8f\n", avg_potential);
    printf("mc final_potential: %3.8f\n\n", final_potential);
}

int count_spin(int size, int* lattice, int spintype) {
    int result = 0;
    int i;
    for (i = 0; i < size; i++) {
        if (lattice[i] == spintype) {
            result++;
        }
    }
    return result;
}

void montecarlo(int n, int m, int* lattice, int** connectivity, int num_steps) {
    int location, j;
    int size = n*m;
    double avg_potential = 0;
    double avg_magnetism = 0;
    int upspin = 1;
    //printf("mc inital potential: %d\n", potential(lattice, connectivity));
    for (j = 0; j < num_steps; j++) {
        for (location = 1; location <= size; location++) {
            int p1 = potential(lattice, connectivity);
            lattice[location - 1] *= -1;
            int p2 = potential(lattice, connectivity);
            int potential_difference = p2 - p1;
            int flip = flipstate(potential_difference);
            lattice[location - 1] *= -flip;
            avg_potential += potential(lattice, connectivity);

            int spinups = count_spin(size, lattice, upspin);
            int spindowns = size - spinups;
            double magnetism = (double) abs(spinups - spindowns) / size;
            avg_magnetism += magnetism;
        }
    }
    avg_potential /= size * num_steps;
    avg_magnetism /= size * num_steps;
    double final_potential = (double) potential(lattice, connectivity);

    printf("%3.8f, %3.8f\n", avg_potential, avg_magnetism);
}

void montecarlo_d(int n, int m, int* lattice, int** connectivity, int num_steps, char* outpath) {
    int sum_potential = 0;
    int location, j;
    int size = n*m;
    double avg_potential = 1;
    int initial_potential = (double) potential(lattice, connectivity);

    for (j = 0; j < num_steps; j++) {
        sum_potential = 0;
        for (location = 1; location <= size; location++) {
            int potential = potential_e(location, lattice, connectivity);
            //            printf("%d\n", potential);
            int potential_difference = -2 * potential;
            int flip = flipstate(potential_difference);
            lattice[location - 1] *= flip;
            int adder = (flip - 1) / 2;

            sum_potential += adder*potential_difference;
        }
        //        printf("sum potential: %d\n", sum_potential);
        avg_potential += (double) sum_potential / size;
    }
    avg_potential *= initial_potential;
    int final_potential = (double) potential(lattice, connectivity);
    printf("avg potential: %3.8f\n", avg_potential);
    printf("final_potential: %3.8f\n\n", final_potential);
}

void print_intgrid(int n, int m, int* lattice) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%3.d ", lattice[i * m + j]);
        }
        printf("\n");
    }
}

void h03test07() {
    int i, j;
    int n = 6;
    int m = 6;
    spin_random(100);
    int* lattice = generate_lattice(n, m);
    for (j = 0; j < n; j += 2) {
        for (i = 0; i < m; i += 2) {
            lattice[j * m + i] = 1;
            lattice[j * m + i + 1] = -1;
        }
    }
    for (j = 1; j < n; j += 2) {
        for (i = 0; i < m; i += 2) {
            lattice[j * m + i] = -1;
            lattice[j * m + i + 1] = +1;
        }
    }
    int** connections = get_connectivity(n, m);
    int u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    // config2
    for (j = 0; j < n; j += 2) {
        for (i = 0; i < m; i += 2) {
            lattice[j * m + i] = -1;
            lattice[j * m + i + 1] = 1;
        }
    }
    for (j = 1; j < n; j += 2) {
        for (i = 0; i < m; i += 2) {
            lattice[j * m + i] = 1;
            lattice[j * m + i + 1] = -1;
        }
    }
    connections = get_connectivity(n, m);
    u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    // config 3
    for (i = 0; i < m * n; i += 2) {
        lattice[ i] = 1;
        lattice[i + 1] = -1;
    }
    connections = get_connectivity(n, m);
    u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    // config 4
    for (i = 0; i < m * n; i += 2) {
        lattice[ i] = -1;
        lattice[i + 1] = 1;
    }
    connections = get_connectivity(n, m);
    u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    // config 5
    for (i = 0; i < m * n; i += 2) {
        lattice[ i] = 1;
        lattice[i + 1] = 1;
    }
    connections = get_connectivity(n, m);
    u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    // config 6
    for (i = 0; i < m * n; i += 2) {
        lattice[ i] = -1;
        lattice[i + 1] = -1;
    }
    connections = get_connectivity(n, m);
    u = potential(lattice, connections);
    printf("lattice potential: %d\n", u);
    print_intgrid(n, m, lattice);

    free(lattice);
    //print_intint(connections);
    free(connections);
}

int* copy_intarray(int size, int* array) {
    int* result = malloc(size * sizeof (int));
    int i;
    for (i = 0; i < size; i++) {
        result[i] = array[i];
    }
    return result;
}

void h03test08() {
    int n = 6;
    int m = 6;
    spin_random(100);
    int* lattice = generate_lattice(n, m);
    int** connections = get_connectivity(n, m);
    int u = potential(lattice, connections);
    int i;
    int num_temps = 20;
    int temp_step = 1;
    int num_loops = 1000;
    T = 1;
    int* lattice_copy;
    lattice_copy = copy_intarray(n*m, lattice);
    u = potential(lattice_copy, connections);
    printf("inital lattice potential: %d\n", u);

    print_intgrid(n, m, lattice_copy);
    printf("\n");
    for (i = 0; i < num_temps; i++) {
        printf("temperature: %d\n", T);
        lattice_copy = copy_intarray(n*m, lattice);
        montecarlo(n, m, lattice_copy, connections, num_loops);
        print_intgrid(n, m, lattice_copy);
        printf("\n\n");
        free(lattice_copy);
        T++;
    }
    free(lattice);
    free(connections);
}

void h03test09() {
    int n = 6;
    int m = 6;
    int size = n*m;
    spin_random(100);
    int* lattice = generate_lattice(n, m);
    int** connections = get_connectivity(n, m);
    int u = potential(lattice, connections);
    int i;
    int num_temps = 30;
    double temp_step = 0.2;
    int num_loops = 1000;
    T = 1;
    int* lattice_copy;
    lattice_copy = copy_intarray(n*m, lattice);
    u = potential(lattice_copy, connections);
    printf("inital lattice potential: %d\n", u);
    int spinups = count_spin(size, lattice, 1);
    int spindowns = size - spinups;
    double magnetism = (double) abs(spinups - spindowns) / size;
    printf("inital lattice magnetization: %3.8f\n", magnetism);
    print_intgrid(n, m, lattice_copy);
    printf("\nTemperature, Potential, Magnetization\n");
    for (i = 0; i < num_temps; i++) {
        printf("%3.2f, ", T);
        lattice_copy = copy_intarray(n*m, lattice);
        montecarlo(n, m, lattice_copy, connections, num_loops);
        free(lattice_copy);
        T+= temp_step;
    }
    free(lattice);
    free(connections);
}

/*
 * 
 */
int hw03main(int argc, char** argv) {
    //h03test01(); // populate lattice
    //    h03test02(); // test intint array print
    //    h03test03(); // test connectivity
    //h03test04(); // test potential
    //h03test05(); // test lattice potential
    // h03test06(); // flip the spin
    //h03test07(); // test configurations and potentials
    // h03test08(); // monte carlo verbose
    h03test09(); // monte carlo data form
    return (EXIT_SUCCESS);
}

