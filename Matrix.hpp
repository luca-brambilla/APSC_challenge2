#include <cstddef>
#include <cstdio>
#include <map>
#include <array>
#include <vector>
#include <iostream>
#include <cmath>

#include <string>
#include <fstream>
#include <sstream>

#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace algebra{

// tolerance to consider a number 0 for sparse matrix
constexpr double ZERO_TOL = 1e-8;

// enumerator for storage order
enum Order {Column, Row};
enum Norm {One, Infinity, Frobenius};

// forward declaration matrix class
template <typename T, typename StorageOrder>
class Matrix;

// forward declaration of friend function inside class template
template<typename T, typename StorageOrder>
std::vector<T> operator*(const Matrix<T, StorageOrder> &m, const std::vector<T> &v );


// Matrix class template
template <typename T, typename StorageOrder>
class Matrix
{
public:
    typedef std::array<std::size_t,2> indexes;
    typedef std::vector<std::vector<T>> uncompressed;
    typedef std::map<indexes,T> coo_matrix;
    //! not all vector<T>, indices are size_t
    typedef std::array<std::vector<T>,3> cs_matrix;

public:
    // constructor
    Matrix();
    Matrix(uncompressed const &v);
    Matrix(Matrix const &r);
    Matrix(std::string const &name);

    // getters
    std::size_t ncols() { return ncol; };
    std::size_t nrows() { return nrow; };

    // compression utilities
    void compress();
    void uncompress();
    bool is_compressed() const;

    void print() const;
    double norm(Norm const &n) const;

    // operations
    friend std::vector<T> operator*<T,StorageOrder>(Matrix<T,StorageOrder> const &m, std::vector<T> const &v );
    friend Matrix<T,StorageOrder> operator*( Matrix<T,StorageOrder> const &m1, Matrix const &m2);

    // access operator
    T operator[] (indexes const &i) const;
    T& operator[] (indexes const &i);

    
private:
    bool compressed = false;
    coo_matrix dynamic_data;
    cs_matrix compressed_data;

    std::size_t ncol = 0;
    std::size_t nrow = 0;
};




//! costructor starting from a full matrix
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


//! constructor reading from file
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(std::string const &name)
/*
read data from file
*/
{
    // open the file
    std::ifstream file(name);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    
    // read first commented lines
    std::string line;
    getline(file, line);
    while(line[0] == '%' )
    {
        getline(file, line);
    }

    // read number of rows and columns
    std::istringstream iss(line);
    iss >> nrow >> ncol;
    std::cout << nrow << " " << ncol << std::endl;

    // read each line from the file
    while (getline(file, line))
    {
        // string stream from the line
        std::istringstream iss(line);
        
        // hold data for each line
        std::size_t i;
        std::size_t j;
        T num;

        // read data from the line
        if (iss >> i >> j >> num) {
            // store value if above tolerance
            if (std::abs(num) > ZERO_TOL)
            {
                dynamic_data.insert( { {i,j}, num} );
            }
        }
        else {
            std::cerr << "Error reading line: " << line << std::endl;
        }
    }

    // close the file
    file.close();
}

//! is_compressed
template<typename T, typename StorageOrder>
bool Matrix<T, StorageOrder>::is_compressed() const
{
    if (compressed)
        return true;
    else
        return false; 
}

//! compress
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

            // add same element to empty rows
            auto it_tmp = std::prev(it);
            for( std::size_t k=it_tmp->first[0]+1; k<temp_ind; ++k )
            {
                compressed_data[2].push_back(compressed_data[0].size()-1);
            }

            //* rows
            temp_ind = it->first[0];
            compressed_data[2].push_back(compressed_data[0].size()-1);
        }
    }
    compressed_data[2].push_back(compressed_data[0].size()-1);
    dynamic_data.clear();
}


//! uncompress

template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::uncompress()
{
    if (!compressed)
    {
        std::cout << "Matrix is already uncompressed" << std::endl;
        return;
    }

    // compressed flag
    compressed = false;

    // sizes of arrays
    std::size_t i_sz = compressed_data[2].size();
    std::size_t j_sz = compressed_data[1].size();

    // insert elements in map
    for (std::size_t i=0; i<i_sz; ++i)
    {
        // use I vector to loop from index i to i+1 in vector J and data
        for (std::size_t j=compressed_data[2][i]; j<compressed_data[2][i+1]; ++j)
        {
            // {col,     row}, data
            dynamic_data.insert( { {i, compressed_data[1][j]}, compressed_data[0][j]} );
        }
    }

    // insert last element
    dynamic_data.insert( { {nrow-1, compressed_data[1][j_sz-1]}, compressed_data[0][j_sz-1]} );
    compressed_data[0].clear();
    compressed_data[1].clear();
    compressed_data[2].clear();

}


//! print
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::print() const
{

    // print if compressed format
    if (compressed)
    {

        std::cout << "print compressed" << std::endl;

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


//! norm
// NOT working for both compressed and uncompressed
template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm(Norm const &n) const
{
    double res = 0.;
    double sum = 0.;
    std::size_t id_tmp = dynamic_data.begin()->first[0];
    std::vector<T> sums(ncol);


    switch (n) {
        case Norm::One:
            // max of sum by columns

            // save sum of each row
            for(auto it = dynamic_data.cbegin(); it != dynamic_data.cend(); ++it)
            {
                sums[it->first[1]] += std::abs(it->second);
            }

            // find maximum among all elements of sums
            for(std::size_t i = 0; i < sums.size(); ++i)
            {
                if (sums[i] > res)
                {
                    res = sums[i];
                }
            }

            break;

        case Norm::Infinity:
            // max of sum by rows
            for(auto it = dynamic_data.cbegin(); it != dynamic_data.cend(); ++it)
            {
                // change row - reset sum - change maximum
                if (it->first[0] != id_tmp)
                {
                    if (sum > res)
                    {
                        res = sum;
                    }
                    sum = 0.;
                    id_tmp = it->first[0];
                }
                
                // update sum
                sum += std::abs(it->second);
            }

            // check last sum for maximum
            if (sum > res)
            {
                res = sum;
            }
            break;

        case Norm::Frobenius:
            // sum of all elements squared
            for(auto it = dynamic_data.cbegin(); it != dynamic_data.cend(); ++it)
            {
                res += it->second * it->second;
            }

            res = std::sqrt(res);
            break;
    }

    return res;
}


// operator[] access copy
template<typename T, typename StorageOrder>
T Matrix<T, StorageOrder>::operator[] (indexes const &i) const
{
    if (compressed)
    {
        // first element of selected row
        std::size_t elem_first = compressed_data[2][ i[0] ];
        return compressed_data[0][ elem_first + i[1]];
    }

    if (dynamic_data.count(i))
    {
        return dynamic_data.at(i);
    }

    if (i[0]>nrow or i[1]>ncol)
    {
        std::cerr << "out of bound index" << std::endl; 
    }
    
    return 0;

}

// operator[] access assign
template<typename T, typename StorageOrder>
T& Matrix<T, StorageOrder>::operator[] (indexes const &i)
{
    if (compressed)
    {
        // first element of selected row
        std::size_t elem_first = compressed_data[2][ i[0] ];
        return compressed_data[0][ elem_first + i[1]];
    }

    // if (dynamic_data.count(i))
    // {
    //     return dynamic_data.at(i);
    // }

    // create new element
    // returns 0 if value is not assigned in usage
    return dynamic_data[i];
}


template<typename T, typename StorageOrder>
std::vector<T> operator*(Matrix<T,StorageOrder> const &m, std::vector<T> const &v )
{

    std::size_t v_sz = v.size();
    if ( m.nrow != v_sz )
    {
        std::cerr << "sizes are not compatible for multiplication: (" 
            << m.nrow << ", " << m.ncol << ") * (" << v.size() << ", 1)"
            << std::endl;
        
        std::vector<T> res;
        return res;
    }

    std::vector<T> res(m.nrow);


    //* compressed format
    if (m.compressed)
    {

        //&Matrix<T, StorageOrder>::cs_matrix data = m.compressed_data;
        // sizes of arrays
        std::size_t i_sz = m.compressed_data[2].size();
        std::size_t j_sz = m.compressed_data[1].size();

        // insert elements in map
        for (std::size_t i=0; i<i_sz; ++i)
        {
            std::size_t k = m.compressed_data[2][i];
            // use I vector to loop from index i to i+1 in vector J and data
            for (std::size_t j=0; k+j<m.compressed_data[2][i+1]; ++j)
            {
                // {col, row}
                res[i] += m.compressed_data[0][k+j] * v[j];

                //std::cout << i << ", " << j
                //        << m.compressed_data[0][k+j] << " " << v[j] << std::endl;
            }

        }
        res[ res.size()-1 ] += m.compressed_data[0][j_sz-1] * v[v_sz-1];
        
        return res;
    }

    //* uncompressed format
    for (auto it=m.dynamic_data.cbegin(); it!=m.dynamic_data.cend(); ++it)
    {
        // row index
        size_t i = it->first[0];
        // column index
        size_t j = it->first[1];
        // matrix value
        T val_ij = it->second;

        // partial multiplication
        res[i] += val_ij * v[j];
    }

    return res;
}

} // namespace algebra

#endif