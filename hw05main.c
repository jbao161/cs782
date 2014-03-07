/* 
 * File:   hw04main.c
 * Author: jbao
 *
 * Created on February 22, 2014, 2:51 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static double PI = 3.14159265358979323846;
double random_number();
void spin_random(int);

/*
 * 
 */

double potential_sphere(double* sphere) {
    double result = 0;
    int size = (int) sphere[0];
    int i, j;
    for (i = 1; i <= size; i++) {
        int xyz1_index = 1 + (-1 + i) * 3;
        double x1 = sphere[xyz1_index];
        double y1 = sphere[xyz1_index + 1];
        double z1 = sphere[xyz1_index + 2];
        for (j = i + 1; j <= size; j++) {
            int xyz2_index = 1 + (-1 + j) * 3;
            double x2 = sphere[xyz2_index];
            double y2 = sphere[xyz2_index + 1];
            double z2 = sphere[xyz2_index + 2];
            double distance = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
            result += 1.0 / pow(distance, 12) - 1.0 / pow(distance, 6);
        }
    }
    result *= 4;
    return result;
}

double* create_sphere_polar(int num_points, double radius) {
    double* result = malloc(1 + 2 * num_points * sizeof (double));
    result[0] = (double) num_points;
    int i;

    for (i = 1; i <= 2 * num_points; i += 2) {
        double phi = random_number()*2; // 0 to 2 in units of PI
        double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
        result[i] = phi;
        result[i + 1] = theta;
    }
    return result;
}

double* create_sphere_euclidean(int num_points, double radius) {

    double* result = malloc(1 + 3 * num_points * sizeof (double));
    result[0] = (double) num_points;
    double potential;
    int i;
start:
    for (i = 1; i < 3 * num_points; i += 3) {
        double phi = random_number()*2; // 0 to 2 in units of PI
        double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
        double x = radius * sin(theta * PI) * cos(phi * PI);
        double y = radius * sin(theta * PI) * sin(phi * PI);
        double z = radius * cos(theta * PI);
        result[i] = x;
        result[i + 1] = y;
        result[i + 2] = z;
    }
    // last point is at origin
    result[3 * num_points - 2] = 0;
    result[3 * num_points - 1] = 0;
    result[3 * num_points ] = 0;

    potential = potential_sphere(result);
    if (potential > 0) {
        goto start;
    }
    return result;
}

double* perfect_sphere(int num_points, double radius) {

}

void print_doubleset(double* set) {
    int size = (int) set[0];
    printf("number of points: %d\n", size);
    int i;
    for (i = 1; i <= size; i++) {
        printf("%3.8f\n", set[i]);
    }
}

void out_sphere_polar(double radius, double* set, char* filepath) {
    int size = (int) set[0];
    int i;
    FILE* f = fopen(filepath, "w");
    for (i = 1; i <= 2 * size; i += 2) {
        double phi = set[i];
        double theta = set[i + 1];
        double x = radius * sin(theta * PI) * cos(phi * PI);
        double y = radius * sin(theta * PI) * sin(phi * PI);
        double z = radius * cos(theta * PI);
        fprintf(f, "%3.8f,%3.8f,%3.8f\n", x, y, z);
    }
    fclose(f);
}

void out_sphere_euclidean(double* set, char* filepath) {
    int size = (int) set[0];
    int i;
    FILE* f = fopen(filepath, "w");
    if (f == 0) {
        //something went wrong while opening
        printf("%s file could not be opened", filepath);
        return;
    }
    fprintf(f, "%d\n\n", size);
    for (i = 1; i <= 3 * size; i += 3) {
        double x = set[i];
        double y = set[i + 1];
        double z = set[i + 2];
        fprintf(f, "H %3.8f %3.8f %3.8f\n", x, y, z);
    }
    if (fclose(f) != 0) {
        //something went wrong while opening
        printf("%s file could not be closed", filepath);
    }
}

void h04test01() {
    double* sphere = create_sphere_polar(120, 10);
    print_doubleset(sphere);
    free(sphere);
}

void h04test02() {
    int size = 13;
    double radius = 2.5;
    double* sphere = create_sphere_euclidean(size, radius);
    out_sphere_euclidean(sphere, "C:/Users/jbao/Dropbox/CS782/hw04/sphere_positions.xyz");
    int i;
    for (i = 1; i <= size; i++) {
        int index = 1 + (i - 1)*3;
        double x = sphere[index];
        double y = sphere[index + 1];
        double z = sphere[index + 2];
        double distance = sqrt(x * x + y * y + z * z);
        printf("%3.7f \n", distance);
    }
    free(sphere);
}

void h04test03() {
    double* sphere = create_sphere_euclidean(4, 1.3);
    out_sphere_euclidean(sphere, "C:/Users/jbao/Dropbox/CS782/hw04/sphere_positions.xyz");
    double potential = potential_sphere(sphere);
    printf("potential: %3.8f", potential);
    free(sphere);
}

double potential_difference(int particle_id, double x, double y, double z, double* sphere) {
    double potential1 = 0;
    int size = (int) sphere[0];
    int j;
    int particle_index = (-1 + particle_id) * 3 + 1;
    double x1 = sphere[particle_index];
    double y1 = sphere[particle_index + 1];
    double z1 = sphere[particle_index + 2];
    for (j = 1; j < particle_id; j++) {
        int xyz2_index = 1 + (j - 1) * 3;
        double x2 = sphere[xyz2_index];
        double y2 = sphere[xyz2_index + 1];
        double z2 = sphere[xyz2_index + 2];
        double distance = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
        potential1 += 1.0 / pow(distance, 12) - 1.0 / pow(distance, 6);
    }
    for (j = particle_id + 1; j <= size; j++) {
        int xyz2_index = 1 + (j - 1) * 3;
        double x2 = sphere[xyz2_index];
        double y2 = sphere[xyz2_index + 1];
        double z2 = sphere[xyz2_index + 2];
        double distance = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
        potential1 += 1.0 / pow(distance, 12) - 1.0 / pow(distance, 6);
    }
    potential1 *= 4;

    double potential2 = 0;
    for (j = 1; j < particle_id; j++) {
        int xyz2_index = 1 + (j - 1) * 3;
        double x2 = sphere[xyz2_index];
        double y2 = sphere[xyz2_index + 1];
        double z2 = sphere[xyz2_index + 2];
        double distance = sqrt((x2 - x)*(x2 - x) + (y2 - y)*(y2 - y) + (z2 - z)*(z2 - z));
        potential2 += 1.0 / pow(distance, 12) - 1.0 / pow(distance, 6);
    }
    for (j = particle_id + 1; j <= size; j++) {
        int xyz2_index = 1 + (j - 1) * 3;
        double x2 = sphere[xyz2_index];
        double y2 = sphere[xyz2_index + 1];
        double z2 = sphere[xyz2_index + 2];
        double distance = sqrt((x2 - x)*(x2 - x) + (y2 - y)*(y2 - y) + (z2 - z)*(z2 - z));
        potential2 += 1.0 / pow(distance, 12) - 1.0 / pow(distance, 6);
    }
    potential2 *= 4;

    double result = potential2 - potential1;
    return result;
}

void h04test04() {
    int num_particles = 13;
    int particle_id = 1;
    double radius = 1.3;
    double* sphere = create_sphere_euclidean(num_particles, radius);
    out_sphere_euclidean(sphere, "C:/Users/jbao/Dropbox/CS782/hw04/before_positions.xyz");
    double potential1 = potential_sphere(sphere);
    printf("before potential: %3.8f\n", potential1);
    double phi = random_number()*2; // 0 to 2 in units of PI
    double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
    double x = radius * sin(theta * PI) * cos(phi * PI);
    double y = radius * sin(theta * PI) * sin(phi * PI);
    double z = radius * cos(theta * PI);

    int particle_index = 1 + (-1 + particle_id) * 3;
    double pd = potential_difference(particle_id, x, y, z, sphere);

    sphere[particle_index] = x;
    sphere[particle_index + 1] = y;
    sphere[particle_index + 2] = z;
    out_sphere_euclidean(sphere, "C:/Users/jbao/Dropbox/CS782/hw04/after_positions.xyz");
    double potential2 = potential_sphere(sphere);
    printf("after potential: %3.8f\n", potential2);
    printf("single element potential difference: %3.8f\n", pd);
    pd = potential2 - potential1;
    printf("all elements potential difference: %3.8f\n", pd);
    pd = potential2 - potential1;
    free(sphere);
}
static double T = 0.3;
static double mc_stepsize = 0.02;

double monte_carlo_brute(double* sphere, int num_loops) {
    int size = sphere[0];
    double sum_potential = 0;
    int particle_id;
    int i, k;
    double initial_potential = potential_sphere(sphere);

    double sd_mean = 0;
    double sd_avg = 0;
    double sq_potential = 0;

    for (i = 0, k = 0; i < num_loops; i++) {
        for (particle_id = 1; particle_id <= size; particle_id++) {
            k++;
            int particle_index = 1 + 3 * (particle_id - 1);
            double phi = random_number()*2; // 0 to 2 in units of PI
            double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
            double dx = mc_stepsize * sin(theta * PI) * cos(phi * PI);
            double dy = mc_stepsize * sin(theta * PI) * sin(phi * PI);
            double dz = mc_stepsize * cos(theta * PI);
            double x = sphere[particle_index];
            double y = sphere[particle_index + 1];
            double z = sphere[particle_index + 2];

            //            printf("%3.6f ,%3.6f ,%3.6f --> %3.6f ,%3.6f ,%3.6f \n", x, y, z, x + dx, y + dy, z + dz);
            x += dx;
            y += dy;
            z += dz;
            double pd = potential_difference(particle_id, x, y, z, sphere);
            //            printf("pd: %3.8f\n", pd);

            if (pd < 0) {// new step was accepted]
accept_new:
                sphere[particle_index] = x;
                sphere[particle_index + 1] = y;
                sphere[particle_index + 2] = z;

                sum_potential += initial_potential + pd;
                initial_potential += pd;
            } else {
                double accept_probability = exp(-pd / T);
                double random_chance = random_number();
                if (random_chance < accept_probability) {
                    goto accept_new;
                }
            }
            double potential = potential_sphere(sphere);
            sq_potential += potential * potential;
            double sd_mean2 = sd_mean + (potential - sd_mean) / k;
            sd_avg = sd_avg + (potential - sd_mean)*(potential - sd_mean2);
            sd_mean = sd_mean2;
        }
    }
    double avg_potentialsq = sq_potential / (num_loops * size);

    sd_avg = sqrt(sd_avg / (k - 1));
    printf("sd %3.8f\n", sd_avg);
    double avg_potential = sum_potential / (num_loops * size);
    double brute_sd = sqrt(avg_potentialsq - avg_potential * avg_potential);
    printf("bsd %3.8f\n", brute_sd);
    return avg_potential;
}

double monte_carlo(double* sphere, int num_loops) {
    int size = sphere[0];
    double sum_potential = 0;
    int particle_id;
    int i, k;
    double potential = potential_sphere(sphere);

    double sd_mean = 0;
    double sd_avg = 0;
    double sq_potential = 0;

    for (i = 0, k = 0; i < num_loops; i++) {
        for (particle_id = 1; particle_id <= size; particle_id++) {
            k++;
            int particle_index = 1 + 3 * (particle_id - 1);
            double phi = random_number()*2; // 0 to 2 in units of PI
            double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
            double dx = mc_stepsize * sin(theta * PI) * cos(phi * PI);
            double dy = mc_stepsize * sin(theta * PI) * sin(phi * PI);
            double dz = mc_stepsize * cos(theta * PI);
            double x = sphere[particle_index];
            double y = sphere[particle_index + 1];
            double z = sphere[particle_index + 2];

            x += dx;
            y += dy;
            z += dz;
            double pd = potential_difference(particle_id, x, y, z, sphere);

            if (pd < 0) {// new step was accepted]
accept_new:
                sphere[particle_index] = x;
                sphere[particle_index + 1] = y;
                sphere[particle_index + 2] = z;

                sum_potential += potential + pd;
                potential += pd;
            } else {
                double accept_probability = exp(-pd / T);
                double random_chance = random_number();
                if (random_chance < accept_probability) {
                    goto accept_new;
                }
                sum_potential += potential;
            }
            double sd_mean2 = sd_mean + (potential - sd_mean) / k;
            sd_avg = sd_avg + (potential - sd_mean)*(potential - sd_mean2);
            sd_mean = sd_mean2;
        }
    }
    sd_avg = sqrt(sd_avg / (k - 1));
    printf("sd %3.8f\n", sd_avg);
    return sd_mean;
}

double monte_carlo_outpotentials(double* sphere, int num_loops, char* filepath) {
    int size = sphere[0];
    double sum_potential = 0;
    int particle_id;
    int i, k;
    double potential = potential_sphere(sphere);

    double sd_mean = 0;
    double sd_avg = 0;
    double sq_potential = 0;

    FILE* f = fopen(filepath, "w");

    for (i = 0, k = 0; i < num_loops; i++) {
        for (particle_id = 1; particle_id <= size; particle_id++) {
            k++;
            int particle_index = 1 + 3 * (particle_id - 1);
            double phi = random_number()*2; // 0 to 2 in units of PI
            double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
            double dx = mc_stepsize * sin(theta * PI) * cos(phi * PI);
            double dy = mc_stepsize * sin(theta * PI) * sin(phi * PI);
            double dz = mc_stepsize * cos(theta * PI);
            double x = sphere[particle_index];
            double y = sphere[particle_index + 1];
            double z = sphere[particle_index + 2];

            x += dx;
            y += dy;
            z += dz;
            double pd = potential_difference(particle_id, x, y, z, sphere);

            if (pd < 0) {// new step was accepted]
accept_new:
                sphere[particle_index] = x;
                sphere[particle_index + 1] = y;
                sphere[particle_index + 2] = z;

                sum_potential += potential + pd;
                potential += pd;
            } else {
                double accept_probability = exp(-pd / T);
                double random_chance = random_number();
                if (random_chance < accept_probability) {
                    goto accept_new;
                }
                sum_potential += potential;
            }
            double sd_mean2 = sd_mean + (potential - sd_mean) / k;
            sd_avg = sd_avg + (potential - sd_mean)*(potential - sd_mean2);
            sd_mean = sd_mean2;
            fprintf(f, "%3.8f, ", potential);
        }
    }
    sd_avg = sqrt(sd_avg / (k - 1));
    printf("sd %3.8f\n", sd_avg);
    fclose(f);
    return sd_mean;
}

void h04test05() {
    T = 0.3;
    mc_stepsize = 0.02;
    int size = 13;
    double radius = 1.3;
    int num_loops = 1 * 10e3;
    double* sphere = create_sphere_euclidean(size, radius);
    double avg_potential = monte_carlo(sphere, num_loops);
    double potential = potential_sphere(sphere);
    printf("final potential: %3.8f\n", potential);
    char filepath[256];
    sprintf(filepath, "C:/Users/jbao/Dropbox/CS782/hw04/mc_final_%04d.xyz", num_loops);
    out_sphere_euclidean(sphere, filepath);
    printf("avg potential: %3.8f\n", avg_potential);

    free(sphere);
}

void h04test06() {
    T = 0.3;
    mc_stepsize = 0.02;
    int size = 13;
    double radius = 1.3;
    int num_loops = 7 * 10e0;
    double* sphere = create_sphere_euclidean(size, radius);
    char filepath[256];
    sprintf(filepath, "C:/Users/jbao/Dropbox/CS782/hw04/mc_potentials_%04d.txt", num_loops);
    double avg_potential = monte_carlo_outpotentials(sphere, num_loops, filepath);
    double potential = potential_sphere(sphere);
    printf("final potential: %3.8f\n", potential);
    sprintf(filepath, "C:/Users/jbao/Dropbox/CS782/hw04/mc_final_%04d.xyz", num_loops);
    out_sphere_euclidean(sphere, filepath);
    printf("avg potential: %3.8f\n", avg_potential);

    free(sphere);
}

double monte_carlo_report(double* sphere, int num_loops, double* sd) {
    int size = sphere[0];
    double sum_potential = 0;
    int particle_id;
    int i, k;
    double potential = potential_sphere(sphere);

    double sd_mean = 0;
    double sd_avg = 0;
    double sq_potential = 0;

    for (i = 0, k = 0; i < num_loops; i++) {
        for (particle_id = 1; particle_id <= size; particle_id++) {
            k++;
            int particle_index = 1 + 3 * (particle_id - 1);
            double phi = random_number()*2; // 0 to 2 in units of PI
            double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
            double dx = mc_stepsize * sin(theta * PI) * cos(phi * PI);
            double dy = mc_stepsize * sin(theta * PI) * sin(phi * PI);
            double dz = mc_stepsize * cos(theta * PI);
            double x = sphere[particle_index];
            double y = sphere[particle_index + 1];
            double z = sphere[particle_index + 2];

            x += dx;
            y += dy;
            z += dz;
            double pd = potential_difference(particle_id, x, y, z, sphere);

            if (pd < 0) {// new step was accepted]
accept_new:
                sphere[particle_index] = x;
                sphere[particle_index + 1] = y;
                sphere[particle_index + 2] = z;

                sum_potential += potential + pd;
                potential += pd;
            } else {
                double accept_probability = exp(-pd / T);
                double random_chance = random_number();
                if (random_chance < accept_probability) {
                    goto accept_new;
                }
                sum_potential += potential;
            }
            double sd_mean2 = sd_mean + (potential - sd_mean) / k;
            sd_avg = sd_avg + (potential - sd_mean)*(potential - sd_mean2);
            sd_mean = sd_mean2;

        }
    }
    sd_avg = sqrt(sd_avg / (k - 1));
    *sd = sd_avg;
    return sd_mean;
}

double monte_carlo_diagnostic(double* sphere, int num_loops, double* sd, char* filepath) {
    int size = sphere[0];
    double sum_potential = 0;
    int particle_id;
    int i, k;
    double potential = potential_sphere(sphere);

    double sd_mean = 0;
    double sd_avg = 0;
    double sq_potential = 0;

    FILE* f = fopen(filepath, "w");

    for (i = 0, k = 0; i < num_loops; i++) {
        for (particle_id = 1; particle_id <= size; particle_id++) {
            k++;
            int particle_index = 1 + 3 * (particle_id - 1);
            double phi = random_number()*2; // 0 to 2 in units of PI
            double theta = acos(random_number()*2 - 1) / PI; // 0 to 1 in units of PI
            double dx = mc_stepsize * sin(theta * PI) * cos(phi * PI);
            double dy = mc_stepsize * sin(theta * PI) * sin(phi * PI);
            double dz = mc_stepsize * cos(theta * PI);
            double x = sphere[particle_index];
            double y = sphere[particle_index + 1];
            double z = sphere[particle_index + 2];

            x += dx;
            y += dy;
            z += dz;
            double pd = potential_difference(particle_id, x, y, z, sphere);

            if (pd < 0) {// new step was accepted]
accept_new:
                sphere[particle_index] = x;
                sphere[particle_index + 1] = y;
                sphere[particle_index + 2] = z;

                sum_potential += potential + pd;
                potential += pd;
            } else {
                double accept_probability = exp(-pd / T);
                double random_chance = random_number();
                if (random_chance < accept_probability) {
                    goto accept_new;
                }
                sum_potential += potential;
            }
            double sd_mean2 = sd_mean + (potential - sd_mean) / k;
            sd_avg = sd_avg + (potential - sd_mean)*(potential - sd_mean2);
            sd_mean = sd_mean2;
            fprintf(f, "%3.8f, ", potential);
        }
    }
    fclose(f);
    sd_avg = sqrt(sd_avg / (k - 1));
    *sd = sd_avg;
    return sd_mean;
}

double* copy_sphere(double* array) {
    int num_atoms = (int) array[0];
    int size = 1 + 3 * num_atoms;
    double* result = malloc(size * sizeof (double));
    int i;
    for (i = 0; i <= size; i++) {
        result[i] = array[i];
    }
    return result;
}

void h04test07() {
    int size = 13;
    double radius = 1.3;
    double* initial_sphere = create_sphere_euclidean(size, radius);
    print_doubleset(initial_sphere);
    double* sphere = copy_sphere(initial_sphere);
    print_doubleset(sphere);
    double* sphere2 = copy_sphere(initial_sphere);
    print_doubleset(sphere2);
}

void h04test08() {
    int num_temps = 3;
    double temps[3] = {0.3, 0.1, 0.05};

    mc_stepsize = 0.02;
    int size = 13;
    double radius = 1.3;
    int num_loops = 1e4;
    int i;

    spin_random(1000);
    double* initial_sphere = create_sphere_euclidean(size, radius);
    double initial_potential = potential_sphere(initial_sphere);

    char report_path[256];
    sprintf(report_path, "C:/Users/jbao/Dropbox/CS782/hw04/%04d_mc_report.txt", num_loops);
    out_sphere_euclidean(initial_sphere, report_path);

    char initialpotential_path[256];
    sprintf(initialpotential_path, "C:/Users/jbao/HW/CS782/hw04/%04d_initial.xyz", num_loops);
    out_sphere_euclidean(initial_sphere, initialpotential_path);

    FILE* f = fopen(report_path, "a");
    fprintf(f, "\ninitial potential: %3.8f\n\n", initial_potential);
    fclose(f);

    for (i = 0; i < num_temps; i++) {
        T = temps[i];
        char potential_path[256];
        sprintf(potential_path, "C:/Users/jbao/HW/CS782/hw04/%04d_potentials_%2.2f.txt", num_loops, T);
        double* sd = malloc(sizeof (double));
        double* sphere = copy_sphere(initial_sphere);

        double avg_potential = monte_carlo_diagnostic(sphere, num_loops, sd, potential_path);
        double potential = potential_sphere(sphere);

        char filepath[256];
        sprintf(filepath, "C:/Users/jbao/HW/CS782/hw04/%04d_final_T%2.2f.xyz", num_loops, T);
        out_sphere_euclidean(sphere, filepath);

        FILE* f = fopen(report_path, "a");
        fprintf(f, "T: %3.2f\n", T);
        fprintf(f, "final potential: %3.8f\n", potential);
        fprintf(f, "avg potential: %3.8f\n", avg_potential);
        fprintf(f, "avg standard deviation: %3.8f\n\n", *sd);
        fclose(f);
        free(sd);
        sd= NULL;
        free(sphere);
        sphere = NULL;
    }
    free(initial_sphere);
}

int main(int argc, char** argv) {
    //h04test01(); // print polar sphere
    //h04test02(); // output euclidean sphere
    //h04test03(); // potential
    //h04test04(); // potential difference
    //h04test05(); // monte carlo
    //h04test06(); // monte carlo output positions (to check sd, mean calculations)
    //h04test07(); // copy sphere
    h04test08(); // report monte carlo results
    return (EXIT_SUCCESS);
}

void junk4() {

}