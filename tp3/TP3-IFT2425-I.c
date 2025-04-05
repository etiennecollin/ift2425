//------------------------------------------------------
// Module  : TP3-IFT2425-I.c
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
#include <string.h>
#include <time.h>

#include <cstddef>

//-------------------------//
//---- Fonction Pour TP ---//
//-------------------------//

///////////////////////////////////////////////////////////////////////////////
// Question 1.1
///////////////////////////////////////////////////////////////////////////////

// Function to integrate
float f(float x) { return 4 * sqrt(1 - pow(x, 2)); }

// Uses the trapezoidal rule for integration.
// Uses naive summation.
float calculateIntegral(float a, float b, int NBINTERV) {
    float h = (b - a) / NBINTERV;
    float integral = 0.0;

    // We rewrote the sum. As is, it looks like this:
    // float xi_prev = a;
    // for (int i = 1; i <= NBINTERV; i++) {
    //     float xi = a + i * h;
    //     integral += (h * (f(xi) + f(xi_prev))) / 2;
    //     xi_prev = xi;
    // }
    //
    // Instead, the first and last terms are calculated separately
    // and the middle terms are calculated in the loop.
    for (int i = 1; i <= NBINTERV - 1; i++) {
        float xi = a + i * h;
        integral += f(xi);
    }
    integral = h * ((f(a) + f(b)) / 2.0 + integral);

    return integral;
}

///////////////////////////////////////////////////////////////////////////////
// Question 1.2
///////////////////////////////////////////////////////////////////////////////

// Determines the integral on each subinterval and puts each value in an array
float* calculate_array_subintegrals(float a, float b, int NBINTERV) {
    float* array = new float[NBINTERV];
    float h = (b - a) / NBINTERV;
    float xi_prev = a;

    for (int i = 0; i < NBINTERV; i++) {
        float xi = a + i * h;
        array[i] = (h * (f(xi_prev) + f(xi))) / 2;
        xi_prev = xi;
    }

    return array;
}

float pairwise_sum(float* array, int left, int right) {
    // Base case
    if (left == right || left > right) return array[left];

    // Recursive case
    int middle = (left + right) / 2;
    float left_sum = pairwise_sum(array, left, middle);
    float right_sum = pairwise_sum(array, middle + 1, right);
    return left_sum + right_sum;
}

///////////////////////////////////////////////////////////////////////////////
// Question 1.3
///////////////////////////////////////////////////////////////////////////////

// n is the length of the array
float kahan_sum(float* array, size_t n) {
    // Here is the original Kahan algorithm as asked in the question
    // float sum = 0.0;
    // float compensation = 0.0;
    //
    // for (int i = 0; i < n; i++) {
    //     float y = array[i] + compensation;
    //     float temp = sum;
    //     sum += y;
    //     compensation = (temp - sum) + y;
    // }
    //
    // return sum;

    // Here is a "second-order iterative Kahan–Babuška algorithm" by Klein
    // https://link.springer.com/article/10.1007/s00607-005-0139-x
    float sum = 0.0;
    float compensation_sum_1 = 0.0;
    float compensation_sum_2 = 0.0;

    for (int i = 1; i < n; i++) {
        float compensation_1, compensation_2;
        float temp = sum + array[i];

        // First order correction
        if (fabs(sum) >= fabs(array[i])) {
            compensation_1 = (sum - temp) + array[i];
        } else {
            compensation_1 = (array[i] - temp) + sum;
        }

        sum = temp;
        temp = compensation_sum_1 + compensation_1;

        // Second order correction
        if (fabs(compensation_sum_1) >= fabs(compensation_1)) {
            compensation_2 = (compensation_sum_1 - temp) + compensation_1;
        } else {
            compensation_2 = (compensation_1 - temp) + compensation_sum_1;
        }

        // Update the compensation sums
        compensation_sum_1 = temp;
        compensation_sum_2 = compensation_sum_2 + compensation_2;
    }

    // Correct error at the end instead of at each step
    return sum + (compensation_sum_1 + compensation_sum_2);
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char** argv) {
    const double PI = 3.14159265358979323846264338;
    int NBINTERV = 5000000;

    // Question 1.1
    float a = 0;
    float b = 1;
    float integral = calculateIntegral(a, b, NBINTERV);
    float error = fabs(integral - PI);
    float logError = log10(error);
    printf("[1>Given_Order:]  Pi=%10.10f  Er=%10.10f  LogEr=%2.2f  \n", integral, error, logError);

    // Question 1.2
    float* array = calculate_array_subintegrals(a, b, NBINTERV);
    integral = pairwise_sum(array, 0, NBINTERV - 1);
    error = fabs(integral - PI);
    logError = log10(error);
    printf("[2>PairwiseSum:]  Pi=%10.10f  Er=%10.10f  LogEr=%2.2f  \n", integral, error, logError);

    // Question 1.3
    integral = kahan_sum(array, NBINTERV);
    error = fabs(integral - PI);
    logError = log10(error);
    printf("[3>KahanSummat:]  Pi=%10.10f  Er=%10.10f  LogEr=%2.2f  \n", integral, error, logError);

    return 0;
}
