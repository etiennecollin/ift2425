//------------------------------------------------------
// Module  : TP3-IFT2425-III.c
// Author  : Etienne Collin & Justin Villeneuve
// Date    : 2025-01-23
// Version : 1.0
// Language: C++
// Note    :
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS
//------------------------------------------------
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//-------------------------//
//---- Fonction Pour TP ---//
//-------------------------//

float calculate_x_float(float mu, float x0, int n) {
    float x = x0;

    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        x = mu * x * (1 - x);
        sum += sqrt(x);  // TODO: naive sum; maybe change the way we sum
    }

    return (2.0 * n) / sum;
}

double calculate_x_double(double mu, double x0, int n) {
    double x = x0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        x = mu * x * (1 - x);
        sum += sqrt(x);  // TODO: naive sum; maybe change the way we sum
    }

    return (2.0 * n) / sum;
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char** argv) {
    float mu = 4.0;
    int n = pow(10, 7);

    // With floats
    float x_float;
    printf("With floats\n");
    x_float = calculate_x_float(mu, 0.2, n);
    printf("[0.20:>%10.10f]\n", x_float);
    x_float = calculate_x_float(mu, 0.4, n);
    printf("[0.40:>%10.10f]\n", x_float);
    x_float = calculate_x_float(mu, 0.6, n);
    printf("[0.60:>%10.10f]\n", x_float);

    // With doubles
    double x_double;
    printf("With doubles\n");
    x_double = calculate_x_double(mu, 0.2, n);
    printf("[0.20:>%10.10f]\n", x_double);
    x_double = calculate_x_double(mu, 0.4, n);
    printf("[0.40:>%10.10f]\n", x_double);
    x_double = calculate_x_double(mu, 0.6, n);
    printf("[0.60:>%10.10f]\n", x_double);

    return 0;
}
