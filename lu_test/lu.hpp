#include <iostream>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <algorithm>

//-------------------------------------------------------------------------------------
namespace algelin {

void print(std::ostream& os, double const* a, int rows, int cols, int precision = 4);
void lu(double *a, int *p, int n);
void eye(double* a, int n);
void copy(double* a, double *b, int n);
void permuta_rows(double* pa, double* const a, int rows, int cols, int* per);
void low_solve(double* y, double const* a, double* b, int n);
void up_solve(double* x, double const* a, double* b, int n);



} // namespace algelin
//-------------------------------------------------------------------------------------