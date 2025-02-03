//------------------------------------------------------
// Module  : Tp-IFT2425-I.1.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cstdio>
#include <iostream>
#include <new>

//------------------------------------------------
// DEFINITIONS
//------------------------------------------------
#define CARRE(X) ((X) * (X))
#define CUBE(X) ((X) * (X) * (X))

//----------------------------------------------------------
//----------------------------------------------------------
// AUXILARY FUNCTIONS --------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
double f(double x) {
    const double array[] = {0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01};
    const size_t array_size = sizeof(array) / sizeof(array[0]);

    double sum1 = 0;  // sum1 : y_i^x ln(y_i)
    double sum2 = 0;  // sum2 : y_i^x
    double sum3 = 0;  // sum3 : ln(y_i)

    for (int i = 0; i < array_size; i++) {
        sum1 += pow(array[i], x) * log(array[i]);
        sum2 += pow(array[i], x);
        sum3 += log(array[i]);
    }

    return (sum1 / sum2) - (1 / x) - (sum3 / array_size);
}

double derivative_f(double x) {
    double h = 1e-5;
    return (f(x + h) - f(x)) / h;
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char **argv) {
    //---------------------------
    // Algorithme NEWTON
    //---------------------------
    double x1 = 0.25;
    double x2 = x1;
    double delta = 1e-6;
    double epsilon = 1e-6;

    int counter = 1;

    do {
        x1 = x2;
        x2 = x1 - f(x1) / derivative_f(x1);

        printf("Iteration: %d, x1 = %.6f, x2 = %.6f, f(x2) = %.6e, |x1 - x2| = %.6e\n", counter, x1, x2, f(x2),
               fabs(x2 - x1));

        counter++;
    } while (fabs(x2 - x1) >= delta && fabs(f(x2)) >= epsilon && derivative_f(x2) != 0);

    // Print the root
    printf("r=%f\n", x2);

    printf("f(r)=%e\n", f(x2));

    // Retour sans probleme
    printf("Fini...\n");
    return 0;
}
