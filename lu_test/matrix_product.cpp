#include<iostream>
#include<memory>
#include<cstdlib>

//---------------------------------------------------------------------------------------------------------
namespace algelin {

void matrix_product(const double* A, const double* B, double* C, int rows, int mid, int cols)
{//const vai manter o ponteiro no inicio
    for (int i = 0; i < rows; ++i, A += mid) 
    {
      for(int j = 0; j < cols; ++j)
      {
        double cij = 0;
        double const* va = A;
        double const* vb = B + j;
        for(int k = 0; k < mid; ++k, ++va, vb += cols)
          cij += *va * *vb;
        
        *C++ = cij;
      }
    }
}


} // namespace algelin
//-----------------------------------------------------------------------------------