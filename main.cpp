#include <cstddef>
#include <iostream>
#include "Matrix.hpp"

int main()
{
    std::cout << "*** INPUT ***" << std::endl;
    algebra::Matrix<double, algebra::Order>::uncompressed data{ {{1,2,3}, {4,5,6}, {7,8,9}} };
    for (size_t i=0; i<3; i++)
    {
        for (size_t j=0; j<data[i].size(); j++)
        {
            std::cout << i << " " << j << ": ";
            std::cout << data[i][j] << std::endl;
        }
    }

    algebra::Matrix<double, algebra::Order> M(data); 
    std::cout << "*** OUTPUT ***" << std::endl;
    M.print();

    std::cout << "Hello" << std::endl;

    return 0;
}