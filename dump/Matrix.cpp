#include "Matrix.hpp"
#include <array>
#include <cstddef>
#include <iostream>

namespace algebra{

// costructor
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(const std::vector<std::vector<T>> &v)
{
    // loop over matrix (vector of vectors)
    // for (auto iti=v.cbegin(); iti!=v.cend(); ++iti)
    // {
    //     for (auto itj = iti->cbegin(); itj != iti->cend(); ++itj)
    //     {
            
    //     }
    // }

    std::size_t ncol = v[0].size();
    std::size_t nrow = v.size();
    std::array<T, 2> ind;
    for (std::size_t i=0; i<nrow; ++i)
    {
        //! error in output
        ind[0]=i;
        if (v[i].size() != ncol)
        {
            std::cerr << "number of elements in columns is not consistent" << std::endl;
        }
        for (std::size_t j = 0; j<ncol; ++j)
        {
            ind[1]=j;
            dynamic_data.at(ind) = v[i][j];
        }
    }

}


template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::print()
{
    for (auto it=dynamic_data.cbegin(); it!=dynamic_data.cend(); ++it)
    {
        std::cout << it->first[0] << " " << it->first[1] << ": " << it->second << std::endl;
    }
}



} // namespace algebra