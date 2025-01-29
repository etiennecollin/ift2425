#include <iostream>

//------------------------------------------------------
// module  : Tp-IFT2425-I.1.c
// author  : Etienne Collin & Justin Villeneuve
// date    : 2025-01-23
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//

//----------------------------------------------------------
//----------------------------------------------------------
// AUXILARY FUNCTIONS --------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------

double f(double c_mv, double *y, size_t n) // n is the size of array y
{
  // sum1 : y_i^(c_mv) ln(y_i)
  double sum1 = 0;
  for (int i = 0; i < n; i++) {
    sum1 += pow(y[i], c_mv) * log(y[i]);
  }

  // sum2 : y_i^(c_mv)
  double sum2 = 0;
  for (int i = 0; i < n; i++) {
    sum2 += pow(y[i], c_mv);
  }

  // sum3 : ln(y_i)
  double sum3 = 0;
  for (int i = 0; i < n; i++) {
    sum3 += log(y[i]);
  }

  return (sum1 / sum2) - (1 / c_mv) - (sum3 / n);
}

double derivativef(double x, double *y, size_t n) {
  double sum1 = 0; // sum  yi^x * (log(yi))^2
  double sum2 = 0; // sum yi^x
  double sum3 = 0; // sum yi^x * log(yi)

  for (int i = 0; i < n; i++) {
    sum1 += pow(y[i], x) * pow(log(y[i]), 2);
    sum2 += pow(y[i], x);
    sum3 += pow(y[i], x) * log(y[i]);
  }

  return (sum1 * sum2 - pow(sum2, 2)) / pow(sum2, 2);
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

  double delta = pow(10, -5); // TODO : maybe change the tolerance (was not
                              // specified in the homework)
  double epsilon = pow(10, -6);

  while (fabs(x2 - x1) >= delta && fabs(f(x1, y, n)) >= epsilon &&
         derivativef(x1, y, n) != 0) // TODO : add last condition
  {
    x2 = x1;
    x1 = x1 - f(x1, y, n) / derivativef(x1, y, n);
  }

  return x1; // return the root
}
