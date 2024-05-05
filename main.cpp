#include <cstddef>
#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include <chrono>
#include <complex>

int main()
{
    // start clock
    auto start = std::chrono::high_resolution_clock::now();


    algebra::Matrix<double, algebra::Order>::uncompressed data{ {{1,2,3, 4}, 
                                                                    {5,6,7, 8},
                                                                    {0,0,0, 0},
                                                                    {9, 10, 11, 12}} };

    //! read and print from full
    if (false)
    {
        std::cout << "*** GENERIC MATRIX ***" << std::endl;
        
        algebra::Matrix<double, algebra::Order>::uncompressed data{ {{1,2,3}, {4,5,6}, {7,8,9}} };
        for (size_t i=0; i<3; i++)
        {
            for (size_t j=0; j<data[i].size(); j++)
            {
                std::cout << i << " " << j << ": ";
                std::cout << data[i][j] << std::endl;
            }
        }

        std::cout << "*** READ INPUT SPARSE ***" << std::endl;
        algebra::Matrix<double, algebra::Order>::uncompressed data2{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };

        algebra::Matrix<double, algebra::Order> M1(data2, algebra::Order::Row_major); 
        std::cout << "*** OUTPUT UNCOMPRESSED ***" << std::endl;
        M1.print();
    }

    //! read and print from file
    if (false)
    {
        algebra::Matrix<double, algebra::Order> M4("data/zenios.mtx");
        M4.print();
    }

    //! compress()
    if (false)
    {
        algebra::Matrix<double, algebra::Order>::uncompressed data3{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };

        algebra::Matrix<double, algebra::Order> M1(data3); 
        std::cout << "*** COMPRESS ***" << std::endl;
        M1.compress(algebra::Compression::CSR);

        std::cout << "*** OUTPUT COMPRESSED ***" << std::endl;
        M1.print();
    }

    //! uncomrpess()
    if (false)
    {
        algebra::Matrix<double, algebra::Order> M_uncomp(data);
        std::cout << "*** COMPRESS ***" << std::endl;
        M_uncomp.compress(algebra::Compression::CSR);

        std::cout << "*** OUTPUT COMPRESSED ***" << std::endl;
        M_uncomp.print();

        std::cout << "*** UNCOMPRESS ***" << std::endl;
        M_uncomp.uncompress();
        std::cout << "*** OUTPUT UNCOMPRESSED ***" << std::endl;
        M_uncomp.print();
    }

    //! norm()
    if (false)
    {
        algebra::Matrix<double, algebra::Order>::uncompressed data5{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };
        algebra::Matrix<double, algebra::Order> M_norm(data5);
        std::cout << "Norm-1: " << M_norm.norm(algebra::Norm::One) << std::endl;
        std::cout << "Infinity: " << M_norm.norm(algebra::Norm::Infinity) << std::endl;
        std::cout << "Frobenius: " << M_norm.norm(algebra::Norm::Frobenius) << std::endl;
    }

    //! norm() with complex
    if (true)
    {
        std::complex<double> a(1.0, 0.0);
        std::complex<double> b(0.0, 1.0);
        algebra::Matrix<std::complex<double>, algebra::Order>::uncompressed data_complex{ {{a,0,0}, 
                                                                        {0,b,0},
                                                                        {0,0,a+b} } };
        algebra::Matrix<std::complex<double>, algebra::Order> M_norm_complex(data_complex);
        std::cout << "Norm-1: " << M_norm_complex.norm(algebra::Norm::One) << std::endl;
        std::cout << "Infinity: " << M_norm_complex.norm(algebra::Norm::Infinity) << std::endl;
        std::cout << "Frobenius: " << M_norm_complex.norm(algebra::Norm::Frobenius) << std::endl;
    }

    //! subscript operator[]
    if (false)
    {
        algebra::Matrix<double, algebra::Order>::uncompressed data6{ {{1.,0,0, 0}, 
                                                                        {0,0.,0, 0},
                                                                        {0,2.,3., 0},
                                                                        {4., 0, 0, 5.}} };

        algebra::Matrix<double, algebra::Order> M3(data6);
        M3.print();
        //std::cout << "assign new value" << std::endl;
        M3[ {1,2} ] = 15;
        M3[ {6,4} ] = 36;
        std::cout << M3[ {1,2} ] << std::endl;
        std::cout << M3[ {3,4} ] << std::endl;
        M3.print();

        std::cout << "const matrix" << std::endl;
        const algebra::Matrix<double, algebra::Order> Mc(data6);
        auto p = Mc[ {6,4} ];
        std::cout << p << std::endl; 
        Mc.print();
    }

    //! matrix mult
    if(false)
    {
        algebra::Matrix<double, algebra::Order> M_mult(data);
        M_mult.print();
        M_mult.compress(algebra::Compression::CSR);
        std::vector<double> v( {1, 1, 1, 1} );
        //std::vector<double> v_res;
        
        auto v_res = M_mult * v;

        //v_res = M_mult * v;

        for (auto it = v_res.cbegin(); it!= v_res.cend(); ++it)
            std::cout << *(it) << std::endl;
    }

    // end clock
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

    return 0;
}