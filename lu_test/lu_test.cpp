#include<iostream>
#include<memory>
#include<cstdlib>
#include<gtest/gtest.h>
#include "lu.cpp"
#include "matrix_product.cpp"

//-------------------------------------------------------------------------------------
namespace algelin {



void Hilbert(double* a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; ++j, ++a)
        {
            *a = 1.0/ (i + j + 1);
        }  
    }
}

TEST(LUFactorization, PALU)
{
    srand(1110);
    uint64_t track = 0;
    uint64_t counter = 0;
    
    for(int max = 1; max < 10; ++max)
    {
        std::cout  << max << std::endl;

        std::unique_ptr<double[]> ma(new double[max * max]);
        std::unique_ptr<double[]> mb(new double[max * max]);
        std::unique_ptr<double[]> pa(new double[max * max]);
        std::unique_ptr<double[]> ml(new double[max * max]);
        std::unique_ptr<double[]> mu(new double[max * max]);
        std::unique_ptr<double[]> mp(new double[max * max]);
        std::unique_ptr<int[]> vp(new int[max]);

        for(int sample = 0; sample < 1; ++sample)
        {
            double* a = ma.get();
            double* b = mb.get();
            double* p_a = pa.get();
  
            double* l = ml.get();
            double* u = mu.get();
            double* prod= mp.get();
            int* p = vp.get();
        
            double* pa2 = a;
            double* pb2 = b;
        
            lu(a,p,max);
            permuta_rows(p_a,b,max,max,p);

            eye(l,max);
            std::copy_n(a, max * max, u);

            double* d=a;
            double* dd=l;
            double* ddd=u;
            for (int k = 0; k< max-1; ++k)
            {
                double*c =d;
                double* ll=dd;
                double* uu=ddd;
                for (int i = 0; i < max-1-k ; ++i)
                {
                    double q=*(c+=max);
                    *(ll+=max)=q;
                    *(uu+=max)=0.0;
                }
                d+=(max+1);
                dd+=(max+1);
                ddd+=(max+1);
            }
            
            std::cout << "------------------------------" << std::endl;
            print(std::cout, l, max, max, 4);
            std::cout << std::endl;
            print(std::cout, u, max, max, 4);

            double tol = 1.0e-8;
            matrix_product(l,u,prod,max,max,max);
            for (int i = 0; i < max * max; ++i)
            {
                double ei = std::abs(p_a[i] - prod[i]);
                ASSERT_LE(ei, tol);
            }
        }
    }
}

TEST(LU, Solve)
{
    constexpr uint64_t max = 10;
    std::srand(10);
    double b[max];
    double b2[max];
    double x[max];
    double y[max];
    double x2[max];
    double a[max * max];
    int p[max];
    
    for(uint64_t n = 0; n < max; ++n)
    {
      double tol = 1.0e-3;
      Hilbert(a, n);
      lu(a, p, n);  // pa = lu
 
      for(int k = 0; k < 100; ++k)
      {
        std::generate_n(x, n, std::rand);
        matrix_product(a, x, b, n, n, 1);
        // a x = b => pa x = p b
        permuta_rows(b2, b, n, n, p);
        // l u x = p b
        low_solve(y, a, b, n); // resolve l y  = b
        up_solve(x2, a, y, n); // resole  u x =  y
        for(int i = 0; i < n; ++i)
        {
          double ei = fabs(x[i] - x2[i]);
          ASSERT_LE(ei, tol);  
        }

      }
    }

}

} // namespace algelin 
//--------------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
    
}
