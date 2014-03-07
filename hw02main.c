/* 
 * File:   hw02main.c
 * Author: jbao
 *
 * Created on February 11, 2014, 11:05 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

static double a = 16807;
static double xm = 2147483647.0;
static double x0 = 2147483583.0;
static double seed = 127; // a prime not too large

/**
 * Lehmer modulo generator for random numberse
 * @return a psuedorandom number between 0 and 1
 */
double random_number() {
    seed = fmod(a * seed, xm);
    return seed / x0;
}

void spin_random(int n) {
    int i, num_runs = n; // initialize by removing first terms of sequence
    for (i = 0; i < num_runs; i++) {
        seed = fmod(a * seed, xm);
    }
}

void h02test01() { // print out some random numbers
    int num_randoms = 100;
    double random_num;
    int i;
    for (i = 0; i < num_randoms; i++) {
        random_num = random_number();
        printf("%3.8f\n", random_num);
    }
}

void random_number_set(int num_randoms, double* randoms) {
    int i;
    for (i = 0; i < num_randoms; i++) {
        randoms[i] = random_number();
    }
}

void print_set(int size, double* set) {
    int i;
    for (i = 0; i < size; i++) {
        double entry = set[i];
        printf("%3.4f\n", entry);
    }
}

void h02test02() { // print out the set
    int num_randoms = 50;
    double* set = malloc(num_randoms * sizeof (double));
    random_number_set(num_randoms, set);
    print_set(num_randoms, set);
    free(set);
}

double* histogram(int n, double* set, double bin_size, double min, double max, int* num_boxesp) {

    double range = max - min;
    int num_boxes = (int) (range / bin_size + 0.5);
    double* hist_table = calloc(num_boxes, sizeof (double));

    int i, j;
    for (i = 0; i < num_boxes; i++) {
        double a = i * bin_size;
        double b = a + bin_size;
        for (j = 0; j < n; j++) {
            double entry = set[j];
            if (entry > a && entry < b) {
                hist_table[i]++;
            }
        }
    }
    *num_boxesp = num_boxes;
    return hist_table;
}

void print_set_norm(int size, double* set, int normalization) {
    int i;
    double coef = 1 / (double) normalization;
    for (i = 0; i < size; i++) {
        double entry = set[i] * coef;
        printf("%3.4f\n", entry);
    }
}

void h02test03() {// make a histogram of n random numbers
    int num_randoms = 50000;
    double* set = malloc(num_randoms * sizeof (double));
    random_number_set(num_randoms, set);
    double bin_size = .05;
    int num_boxes = (int) (1.0 / bin_size + 0.5);
    double* hist_table = calloc(num_boxes, sizeof (double));

    int i, j;
    for (i = 0; i < num_boxes; i++) {
        double a = i * bin_size;
        double b = a + bin_size;
        for (j = 0; j < num_randoms; j++) {
            double entry = set[j];
            if (entry > a && entry < b) {
                hist_table[i]++;
            }
        }
    }
    print_set_norm(num_boxes, hist_table, num_randoms);
    free(set);
    free(hist_table);
}

void h02test04() {// make a histogram of n random numbers
    int num_randoms = 50000;
    double* set = malloc(num_randoms * sizeof (double));
    random_number_set(num_randoms, set);
    int* num_binsp;
    double bin_size = 0.1;
    double min = 0;
    double max = 1;
    double* hist_set = histogram(num_randoms, set, bin_size, min, max, num_binsp);
    print_set_norm(*num_binsp, hist_set, num_randoms);
    free(set);
    free(hist_set);
}

double* deviation(int size, double* set, double expectation) {
    int i;
    double coef = 100. / expectation;

    double* result = malloc(size * sizeof (double));
    for (i = 0; i < size; i++) {
        double error = fabs(expectation - set[i]);
        result[i] = error * coef;
    }
    return result;
}

double* histogram_n(int n, double* set, double bin_size, double min, double max, int* num_boxesp) {

    double range = max - min;
    int num_boxes = (int) (range / bin_size + 0.5);
    double* hist_table = calloc(num_boxes, sizeof (double));

    int i, j;
    for (i = 0; i < num_boxes; i++) {
        double a = i * bin_size;
        double b = a + bin_size;
        for (j = 0; j < n; j++) {
            double entry = set[j];
            if (entry > a && entry < b) {
                hist_table[i]++;
            }
        }
    }
    *num_boxesp = num_boxes;

    // normalize the histogram to sum to one
    double coef = 1.0 / n;
    for (i = 0; i < num_boxes; i++) {
        hist_table[i] = hist_table[i] * coef;
    }
    return hist_table;
}

void h02test05() { // test the percent deviation from the expected normalized bin frequency
    int num_randoms = 50000;
    double* set = malloc(num_randoms * sizeof (double));
    random_number_set(num_randoms, set);
    int* num_binsp;
    double bin_size = 0.05;
    double min = 0;
    double max = 1;
    double* hist_set = histogram_n(num_randoms, set, bin_size, min, max, num_binsp);
    double* dev = deviation(*num_binsp, hist_set, bin_size);
    print_set(*num_binsp, dev);
    print_set(*num_binsp, hist_set);
    free(set);
    free(hist_set);
    free(dev);
}

void file_out(int size, double* data, char* filepath) {
    FILE *f = fopen(filepath, "w");
    int i;
    for (i = 0; i < size; i++) {
        fprintf(f, "%3.8f\n", data[i]);
    }
    fclose(f);
}

h02test06() { // test the percent deviation from the expected normalized bin frequency
    int num_randoms = (int) pow(10, 6);
    double* set = malloc(num_randoms * sizeof (double));
    random_number_set(num_randoms, set);
    int* num_binsp;
    double bin_size = 0.05;
    double min = 0;
    double max = 1;
    double* hist_set = histogram_n(num_randoms, set, bin_size, min, max, num_binsp);
    double* dev = deviation(*num_binsp, hist_set, bin_size);
    file_out(*num_binsp, dev, "C:/Users/jbao/Dropbox/CS782/hw02/deviations.txt");
    file_out(*num_binsp, hist_set, "C:/Users/jbao/Dropbox/CS782/hw02/histogram.txt");
    free(set);
    free(hist_set);
    free(dev);
}

void randomwalk3d_a1() {
    int num_steps = (int) pow(10, 5);

    int dice_roll;
    int x = 0, y = 0, z = 0;
    int i;
    for (i = 0; i < num_steps; i++) {
        dice_roll = (int) (random_number()*6);
        switch (dice_roll) {
            case 0: x++;
                break;
            case 1: x--;
                break;
            case 2: y++;
                break;
            case 3: y--;
                break;
            case 4: z++;
                break;
            case 5: z--;
                break;
        }
    }
    printf("%d, %d, %d\n", x, y, z);
}

void h02test07() { // random walk
    clock_t cstart = clock();
    clock_t cend = 0;
    randomwalk3d_a1();
    cend = clock();
    printf("%.6f cpu sec\n", ((double) cend - (double) cstart)* 1.0e-6);
}

void randomwalk3d_a2() {
    int num_steps = (int) pow(10, 3);
    int num_walks = 100;
    int dice_roll;
    int x = 0, y = 0, z = 0;
    int i, j;
    FILE *f_position = fopen("C:/Users/jbao/Dropbox/CS782/hw02/rw_positions.txt", "w");
    FILE *f_displace = fopen("C:/Users/jbao/Dropbox/CS782/hw02/rw_displacement.txt", "w");
    fprintf(f_displace, "3d random walk         num_steps: %d  num_walks: %d\n", num_steps, num_walks);
    fprintf(f_position, "3d random walk         num_steps: %d  num_walks: %d\n", num_steps, num_walks);
    fclose(f_position);
    fclose(f_displace);

    double avg_displacement = 0;
    // the walk
    for (j = 0; j < num_walks; j++) {
        x = 0, y = 0, z = 0;
        for (i = 0; i < num_steps; i++) {
            dice_roll = (int) (random_number()*6);
            switch (dice_roll) {
                case 0: x++;
                    break;
                case 1: x--;
                    break;
                case 2: y++;
                    break;
                case 3: y--;
                    break;
                case 4: z++;
                    break;
                case 5: z--;
                    break;
            }
        }
        double displacement = sqrt(x * x + y * y + z * z);
        avg_displacement += displacement;
        FILE *f_position = fopen("C:/Users/jbao/Dropbox/CS782/hw02/rw_positions.txt", "a");
        FILE *f_displace = fopen("C:/Users/jbao/Dropbox/CS782/hw02/rw_displacement.txt", "a");
        fprintf(f_position, "%d, %d, %d\n", x, y, z);
        fprintf(f_displace, "%3.8f\n", displacement);
        fclose(f_position);
        fclose(f_displace);
    }
    f_displace = fopen("C:/Users/jbao/Dropbox/CS782/hw02/rw_displacement.txt", "a");
    avg_displacement /= num_walks;
    fprintf(f_displace, "avg_displacement:%3.8f\n", avg_displacement);
    fclose(f_displace);
}

void h02test08() {
    randomwalk3d_a2();
}

void randomwalk3d_a3() {
    int num_walks = 100;
    int dice_roll;
    int x = 0, y = 0, z = 0;
    int i, j, k;

    int report_steps[11] = {0, (int) pow(10, 3), 5000, 10000, 25000, 40000, 60000, 80000, 100000, 500000, 1000000};
    double report_displace[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // the walk
    for (j = 0; j < num_walks; j++) {
        x = 0, y = 0, z = 0;
        for (k = 1; k <= 10; k++) {
            for (i = report_steps[k - 1]; i < report_steps[k]; i++) {
                dice_roll = (int) (random_number()*6);
                switch (dice_roll) {
                    case 0: x++;
                        break;
                    case 1: x--;
                        break;
                    case 2: y++;
                        break;
                    case 3: y--;
                        break;
                    case 4: z++;
                        break;
                    case 5: z--;
                        break;
                }
            }
            double displacement = sqrt(x * x + y * y + z * z);
            report_displace[k - 1] += displacement;
        }
    }
    for (i = 0; i < 10; i++) {
        report_displace[i] /= num_walks;
    }
    file_out(10, report_displace, "C:/Users/jbao/Dropbox/CS782/hw02/rw_avgdisplace.txt");
}

void h02test09() {
    randomwalk3d_a3();
}

double* sum_set(int size, double* setA, double* setB) {
    int i;
    double* result = malloc(size * sizeof (double));
    for (i = 0; i < size; i++) {
        result[i] = setA[i] + setB[i];
    }
    return result;
}

void h02test10() {
    int size = (int) pow(10, 5);
    double* setA = malloc(size * sizeof (double));
    random_number_set(size, setA);
    double* setB = malloc(size * sizeof (double));
    random_number_set(size, setB);
    double* setXY = sum_set(size, setA, setB);
    free(setA);
    free(setB);
    int* num_binsp;
    double bin_size = 0.05;
    double min = 0;
    double max = 2;
    double* hist_set = histogram_n(size, setXY, bin_size, min, max, num_binsp);
    file_out(*num_binsp, hist_set, "C:/Users/jbao/Dropbox/CS782/hw02/XYhistogram.txt");
    free(setXY);
    free(hist_set);
}

/*
 * 
 */
int main02(int argc, char** argv) {
    //h02test01();
    //h02test02();
    //h02test03();
    //h02test04();
    //h02test05();
    //h02test06(); // histogram
    //h02test07();
    //h02test08(); // large random walk
    //h02test09(); // large random walk with intermediate averages
    //h02test10(); // distribution of sum of two random numbers
    return (EXIT_SUCCESS);
}

void junk() {
    /*
         double* hist_table = histogram(num_randoms, set, 0.1, 0, 1, num_bins);
        print_set(*num_bins, hist_table);
     * 
     */
}