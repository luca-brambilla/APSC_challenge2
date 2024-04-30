#include <iostream>
#include "Matrix.hpp"

int main()
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

    algebra::Matrix<double, algebra::Order> M(data2); 
    std::cout << "*** OUTPUT UNCOMPRESSED ***" << std::endl;
    M.print();

    std::cout << "*** COMPRESS ***" << std::endl;
    M.compress();

    std::cout << "*** OUTPUT COMPRESSED ***" << std::endl;
    M.print();

    return 0;
}