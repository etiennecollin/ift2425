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
double f(double c_mv, double *array, size_t array_size) {
  // sum1 : y_i^(c_mv) ln(y_i)
  double sum1 = 0;
  for (int i = 0; i < array_size; i++) {
    sum1 += pow(array[i], c_mv) * log(array[i]);
  }

  // sum2 : y_i^(c_mv)
  double sum2 = 0;
  for (int i = 0; i < array_size; i++) {
    sum2 += pow(array[i], c_mv);
  }

  // sum3 : ln(y_i)
  double sum3 = 0;
  for (int i = 0; i < array_size; i++) {
    sum3 += log(array[i]);
  }

  return (sum1 / sum2) - (1 / c_mv) - (sum3 / array_size);
}

double derivativef(double x, double *array, size_t array_size) {
  double h = pow(10, -5);
  double num =
      -f(x + 2 * h, array, array_size) + 8 * f(x + h, array, array_size) -
      8 * f(x - h, array, array_size) + f(x - 2 * h, array, array_size);
  double denom = 12 * h;

  return num / denom;
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

  double y[] = {0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01};
  int n = 10;

  double x1 = 0.25;
  double x2 = 0;

  // TODO(grosjuice): Maybe change the tolerance (was not specified in the
  // homework)
  double delta = pow(10, -5);
  double epsilon = pow(10, -6);

  while (fabs(x2 - x1) >= delta && fabs(f(x1, y, n)) >= epsilon &&
         derivativef(x1, y, n) != 0) {
    x2 = x1;
    x1 = x1 - f(x1, y, n) / derivativef(x1, y, n);
  }

  // Print the root
  printf("Root : %f\n", x1);

  // Retour sans probleme
  printf("\n Fini... \n\n\n");
  return 0;
}
