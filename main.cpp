#include <cstddef>
#include <iostream>
#include <vector>
#include "Matrix.hpp"

int main()
{

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

        algebra::Matrix<double, algebra::Order> M1(data2); 
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
        M1.compress();

        std::cout << "*** OUTPUT COMPRESSED ***" << std::endl;
        M1.print();
    }

    //! uncomrpess()
    if (true)
    {
        algebra::Matrix<double, algebra::Order> M_uncomp(data);
        std::cout << "*** COMPRESS ***" << std::endl;
        M_uncomp.compress();

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
        algebra::Matrix<double, algebra::Order> M2(data5); 
        std::cout << "Norm-1: " << M2.norm(algebra::Norm::One) << std::endl;
        std::cout << "Infinity: " << M2.norm(algebra::Norm::Infinity) << std::endl;
        std::cout << "Frobenius: " << M2.norm(algebra::Norm::Frobenius) << std::endl;
    }

    //! subscript operator[]
    if (false)
    {
        algebra::Matrix<double, algebra::Order>::uncompressed data6{ {{1.,0,0, 0}, 
                                                                        {0,2.,0, 0},
                                                                        {0,0,3., 0},
                                                                        {0, 0, 0, 4.}} };

        algebra::Matrix<double, algebra::Order> M3(data6);
        M3.print();
        //std::array<std::size_t,2> i({1,1});
        M3[ {1,2} ] = 5;
        std::cout << M3[ {1,2} ] << std::endl;
        std::cout << M3[ {3,4} ] << std::endl;
        M3.print();

        std::cout << "const matrix" << std::endl;
        const algebra::Matrix<double, algebra::Order> Mc(data6);
        auto p = Mc[ {6,4} ];
        std::cout << p << std::endl; 
        Mc.print();
    }

    return 0;
}