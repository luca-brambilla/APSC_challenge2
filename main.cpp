#include <cstddef>
#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include <chrono>
#include <complex>

int main()
{
    algebra::Matrix<double, algebra::Order>::fullmatrix data{ {{1,2,3, 4}, 
                                                                    {5,6,7, 8},
                                                                    {0,0,0, 0},
                                                                    {9, 10, 11, 12}} };

    //! test fullmatrix
    if (false)
    {
        std::cout << "*** GENERIC MATRIX ***" << std::endl;
        
        algebra::Matrix<double, algebra::Order>::fullmatrix data0{ {{1,2,3}, {4,5,6}, {7,8,9}} };
        for (size_t i=0; i<3; i++)
        {
            for (size_t j=0; j<data[i].size(); j++)
            {
                std::cout << i << " " << j << ": ";
                std::cout << data[i][j] << std::endl;
            }
        }
    }

    //! read fullmatrix and print
    if (false)
    {
        std::cout << "*** READ INPUT SPARSE ***" << std::endl;
        algebra::Matrix<double, algebra::Order>::fullmatrix data_full{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };

        algebra::Matrix<double, algebra::Order> M_full(data_full, algebra::Order::Row_major); 
        std::cout << "*** OUTPUT UNCOMPRESSED ***" << std::endl;
        M_full.print();
    }

    //! read and print from file
    if (false)
    {
        algebra::Matrix<double, algebra::Order> M_file("data/zenios.mtx");
        M_file.print();
    }

    //! compress()
    if (false)
    {
        algebra::Matrix<double, algebra::Order>::fullmatrix data_compress{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };

        algebra::Matrix<double, algebra::Order> M_compress(data_compress); 
        std::cout << "*** COMPRESS ***" << std::endl;
        M_compress.compress(algebra::Compression::CSR);

        std::cout << "*** OUTPUT COMPRESSED ***" << std::endl;
        M_compress.print();
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
        algebra::Matrix<double, algebra::Order>::fullmatrix data_norm{ {{1,2,3, 4}, 
                                                                        {5,6,7, 8},
                                                                        {0,0,0, 0},
                                                                        {9, 10, 11, 12}} };
        algebra::Matrix<double, algebra::Order> M_norm(data_norm);
        std::cout << "Norm-1: " << M_norm.norm(algebra::Norm::One) << std::endl;
        std::cout << "Infinity: " << M_norm.norm(algebra::Norm::Infinity) << std::endl;
        std::cout << "Frobenius: " << M_norm.norm(algebra::Norm::Frobenius) << std::endl;
    }

    //! norm() with complex
    if (false)
    {
        std::complex<double> a(1.0, 0.0);
        std::complex<double> b(0.0, 1.0);
        algebra::Matrix<std::complex<double>, algebra::Order>::fullmatrix data_complex{ {{a,0,0}, 
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
        algebra::Matrix<double, algebra::Order>::fullmatrix data_sub{ {{1.,0,0, 0}, 
                                                                        {0,0.,0, 0},
                                                                        {0,2.,3., 0},
                                                                        {4., 0, 0, 5.}} };

        algebra::Matrix<double, algebra::Order> M_sub(data_sub);
        M_sub.print();
        //std::cout << "assign new value" << std::endl;
        M_sub[ {1,2} ] = 15;
        M_sub[ {6,4} ] = 36;
        std::cout << M_sub[ {1,2} ] << std::endl;
        std::cout << M_sub[ {3,4} ] << std::endl;
        M_sub.print();

        std::cout << "const matrix" << std::endl;
        const algebra::Matrix<double, algebra::Order> M_const(data_sub);
        auto p = M_const[ {6,4} ];
        std::cout << p << std::endl; 
        M_const.print();
    }

    //! vector-matrix product
    if(false)
    {
        algebra::Matrix<double, algebra::Order> M_mult(data);
        M_mult.print();
        M_mult.compress(algebra::Compression::CSR);
        M_mult.print();
        std::vector<double> v( {1, 1, 1, 1} );
        //std::vector<double> v_res;
        
        auto v_res = M_mult * v;

        //v_res = M_mult * v;

        for (auto it = v_res.cbegin(); it!= v_res.cend(); ++it)
            std::cout << *(it) << std::endl;
    }

    //! test chrono
    if (false)
    {
        algebra::Matrix<double, algebra::Order> M_chrono("data/lnsp_131.mtx");
        //M_chrono.print();
        std::vector<double> v_chrono(M_chrono.ncols(), 1.);

        // start clock
        auto start = std::chrono::high_resolution_clock::now();
        // matrix-vector multiplication
        std::vector<double> res_chrono = M_chrono * v_chrono;
        // end clock
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

        // compression
        M_chrono.compress(algebra::Compression::CSR);
        //M_chrono.print();

        // start clock
        start = std::chrono::high_resolution_clock::now();
        //matrix-vector multiplication
        res_chrono = M_chrono * v_chrono;
        // end clock
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

    }
    return 0;
}