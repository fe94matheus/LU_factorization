#include <iostream>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <algorithm>

//-------------------------------------------------------------------------------------
namespace algelin {

void print(std::ostream& os, double const* a, int rows, int cols, int precision)
{
   os << std::setprecision(precision) << std::showpos << std::scientific;
   for(int i = 0; i < rows; ++i)
   {
      for(int j = 0; j < cols; ++j)
         os << *a++ << "  ";
      os << std::endl;
   }
}

void lu(double *a, int *p, int n)
 {
    std::iota(p, p + n, 0) ;
    
    double* d = a;
    for (int k = 0; k < n - 1; ++k)//pivoteamento parcial
    {
        double max = fabs( *d ); 

        int imax = 0;
        double* c = d;

        for (int i = 1; i < n-k; ++i)
        {
            double v = fabs( *(c += n) );
            if(v > max)
            {
                max = v;
                imax = i;
            }
        }

        if( imax )
        {
                std::swap(p[k],p[k+imax]);
                //troca a linha l_0(aponta para o inicio da linha l_0) 
                //pela linha l_n(aponta para o inicio da linha l_n)

                double* l_0 = d - k;
                double* l_n = l_0  + n * imax;
                for (int j = 0; j < n; ++j)
                {
                    std::swap(*l_0++,*l_n++);
                }
                
        }

        c = d;
        for(int i = 0; i < n -1 - k; ++i)
        {
            double q = *(c+=n)/ *d;
            *c=q;
            double *l_0=d;
            double *l_i=c;
            for (int j = k+1; j < n; ++j)
            {
                *(++l_i)-= q * *(++l_0);
            }
            
        }
        d+=(n+1);    
    }
    
    
}

void eye(double* a, int n)
{
    for (int i = 0; i < n-1; i++)
    {
        *a++=1;
        for (int j = 0; j < n; j++)
          *a++=0;
    }
    *a=1;
  
  
}

void permuta_rows(double* pa, double* const a, int rows, int cols, int* per)
{
    for (int i = 0; i < rows; ++i)
    {
        double const* li=a+per[i]*cols;
        for (int j = 0; j < cols; ++j)
        {
            *(pa++) = *(li++);

        }
    }   
}

void low_solve(double* y, double* l, double* b, int n){
    *y = *b;
    double* d = l;
    double* e = b;
    

    for (int k = n; k > 1; --k)
    {
        double sum = 0;
        double* c = d;
        double* w=y;
        c+=k;
        for (int i = 0; i <= n-k; ++i)
        {
            sum += *(c++) * *(w++);
        }
        *(++w) = *(++e) - sum;
        d+=(n+1);
    }
    
}
//Ux=y
void up_solve(double* x, double* u, double* y, int n){
    *(x+(n-1)) = (*(y+(n-1)))/(*(u+((n*n)-1)));
    double* d = u;
    double* e = y;
    

    for (int k = n; k > 1; --k)
    {
        double sum = 0;
        double* c = d;
        double* w=x;
        c-=k;
        for (int i = 0; i <= n-k; ++i)
        {
            sum += *(c--) * *(w--);
        }
        d-=(n+1);
        *(--w) = *(--e) - sum;
    }
}

}  // namespace algelin
//-------------------------------------------------------------------------------------------