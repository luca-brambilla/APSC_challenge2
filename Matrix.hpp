#include <cstddef>
#include <map>
#include <array>
#include <vector>
#include <iostream>

#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace algebra{

// tolerance to consider a number 0 for sparse matrix
constexpr double ZERO_TOL = 1e-6;

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
    bool is_compressed() const;

    void print() const;
    double norm() const;

    // operations
    Matrix matmul();

    // access operator
    T operator[] (std::size_t const &i) const;
    T& operator[] (std::size_t const &i);

    
private:
    bool compressed = false;
    coo_matrix dynamic_data;
    cs_matrix compressed_data;

    std::size_t ncol = 0;
    std::size_t nrow = 0;
};




// costructor starting from a full matrix
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

    ncol = v[0].size();
    nrow = v.size();
    indexes ind;
    for (std::size_t i=0; i<nrow; ++i)
    {
        //! error in output
        ind[0]=i;
        if (v[i].size() != ncol)
        {
            std::cerr << "number of elements in columns is not consistent" << std::endl;
        }

        //! initialize - loop over all elements
        //! could be done with swap? - no... threshold
        for (std::size_t j = 0; j<ncol; ++j)
        {
            if (std::abs(v[i][j]) > ZERO_TOL)
            {
                ind[1]=j;
                dynamic_data.insert( {ind, v[i][j]} );
            }
        }
    }

}

// is_compressed
template<typename T, typename StorageOrder>
bool Matrix<T, StorageOrder>::is_compressed() const
{
    if (compressed)
        return true;
    else
        return false; 
}

// compress
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::compress()
{
    if (compressed)
    {
        std::cout << "Matrix is already compressed" << std::endl;
        return;
    }

    // compressed flag
    compressed = true;
    // keep track of row index change and element
    std::size_t temp_ind = dynamic_data.cbegin()->first[0];
    // first element row
    compressed_data[2].push_back(temp_ind);

    for (auto it=dynamic_data.cbegin(); it!=dynamic_data.cend(); ++it)
    {
        //std::cout << it->first[0] << " " << std::next(it)->first[0];
        //* value
        compressed_data[0].push_back(it->second);
        //* column
        compressed_data[1].push_back(it->first[1]);
        
        // index of new row
        if(it->first[0] != temp_ind)
        {
            temp_ind = it->first[0];

            // add zeros to empty rows
            auto it_tmp = std::prev(it);
            // std::cout << "**" << it_tmp->first[0];
            for( std::size_t k=it_tmp->first[0]+1; k<temp_ind; ++k )
            {
                // std::cout << ".";
                compressed_data[2].push_back(compressed_data[0].size()-1);
            }
            // std::cout<< "**";

            //* rows
            temp_ind = it->first[0];
            compressed_data[2].push_back(compressed_data[0].size()-1);
        }
        // std::cout << " - ";
    }
    // std::cout<<std::endl;
    compressed_data[2].push_back(compressed_data[0].size()-1);

    // for(auto it=compressed_data[2].cbegin(); it!=compressed_data[2].cend(); ++it)
    // {
    //     std::cout << *it << " ";
    // }
    // std::cout << std::endl;
}

// print
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::print() const
{

    // print if compressed format
    if (compressed)
    {

        std::cout << "print compressed" << std::endl;
        // std::cout << compressed_data[0].size() << " "
        //         << compressed_data[1].size() << " " 
        //         << compressed_data[2].size() << std::endl;
        // std::cout << "start print" << std::endl;


        //* i = row index
        std::size_t i=0;

        //* i_sz = dimension of row index vector
        std::size_t i_sz = compressed_data[2].size();

        //* j index along element vector
        for(std::size_t j=0; j<compressed_data[1].size(); )
        {

            // if index of element vector equal to next element, then new rob
            // index of element vector must be different from last element
            if ( (j == compressed_data[2][i+1]) and (j!=compressed_data[2][i_sz-1]) ) 
            {
                ++i;
            }

            // if next index is the same, then it is an empty row
            if (compressed_data[2][i+1] == compressed_data[2][i])
            {
                ++i;
                continue;
            }

            std::cout << i << " " << j << ": " << compressed_data[0][j] << std::endl;

            ++j;
        }

        return;
    }
    // print if not compressed
    else
    {
        std::cout << "print uncompressed" << std::endl;
        for (auto it=dynamic_data.cbegin(); it!=dynamic_data.cend(); ++it)
        {
            std::cout << it->first[0] << " " << it->first[1] << ": " << it->second << std::endl;
        }
    }


}


} // namespace algebra

#endif