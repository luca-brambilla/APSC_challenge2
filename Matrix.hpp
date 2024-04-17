#include <map>
#include <array>
#include <vector>
#include <iostream>

#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace algebra{

// enumerator for storage order
enum Order {Column, Row};

// Matrix class template
template <typename T, typename StorageOrder>
class Matrix
{
public:
    typedef std::array<std::size_t,2> indexes;
    typedef std::vector<std::vector<T>> uncompressed;
    typedef std::map<indexes,T> coo_matrix;
    typedef std::array<std::vector<T>,3> cs_matrix;

public:
    // constructor
    Matrix();
    Matrix(uncompressed const &v);
    Matrix(Matrix const &r);

    // compression utilities
    void compress();
    void uncompress();
    const bool is_compressed();

    void print();
    double norm();

    // operations
    Matrix matmul();

    // access operator
    T operator[] (std::size_t const &i) const;
    T& operator[] (std::size_t const &i);

    
private:
    bool compressed = false;
    coo_matrix dynamic_data;
    cs_matrix compressed_data;
};




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
    indexes ind;
    for (std::size_t i=0; i<nrow; ++i)
    {
        //! error in output
        ind[0]=i;
        if (v[i].size() != ncol)
        {
            std::cerr << "number of elements in columns is not consistent" << std::endl;
        }

        //! initialize
        for (std::size_t j = 0; j<ncol; ++j)
        {
            ind[1]=j;
            dynamic_data.insert( {ind, v[i][j]} );
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

#endif