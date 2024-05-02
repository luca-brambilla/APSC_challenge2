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
#ifndef ZERO_TOL
constexpr double ZERO_TOL = 1e-8;
#endif

// enumerator for storage order

// Coordinates, Compressed Sparse Row, Compressed Sparse Column, Modified Sparse Row, Modified Sparse Column
enum Order {COO, CSR, CSC, MSR, MSC};
enum Norm {One, Infinity, Frobenius};

// forward declaration matrix class
template <typename T, typename StorageOrder>
class Matrix;

// forward declaration of friend function inside class template
template<typename T, typename StorageOrder>
std::vector<T> operator*( Matrix<T, StorageOrder> const &m, std::vector<T> const &v );

template<typename T, typename StorageOrder>
Matrix<T,StorageOrder> operator*( Matrix<T,StorageOrder> const &m1, Matrix<T,StorageOrder> const &m2);

// Matrix class template
template <typename T, typename StorageOrder>
class Matrix
{
public:
    typedef std::array<std::size_t,2> indexes;
    typedef std::vector<std::vector<T>> uncompressed;
    typedef std::map<indexes,T> coo_matrix;

public:
    // constructor
    Matrix();

    Matrix(std::size_t const& r, size_t const& c);

    Matrix(uncompressed const &m);
    Matrix(uncompressed const &m, Order const &o);
    
    Matrix(Matrix const &m);
    
    Matrix(std::string const &name);
    Matrix(std::string const &name, Order const &o);

    // getters
    std::size_t ncols() { return ncol; };
    std::size_t nrows() { return nrow; };

    //
    void resize(std::size_t const& r, size_t const& c);

    // compression utilities
    void compress(Order const &o);
    void uncompress();
    bool is_compressed() const;

    void print() const;

    // norms
    double norm_one() const;
    double norm_infty() const;
    double norm_frob() const;
    double norm(Norm const &n) const;

    // operations
    friend std::vector<T> operator*<T,StorageOrder>(Matrix<T,StorageOrder> const &m, std::vector<T> const &v );
    friend Matrix<T,StorageOrder> operator*<T,StorageOrder>( Matrix<T,StorageOrder> const &m1, Matrix const &m2);

    // access operator
    T operator[] (indexes const &i) const;
    T& operator[] (indexes const &i);

    
private:
    Order ordering = Order::COO;

    bool compressed = false;

    coo_matrix dynamic_data;

    std::vector<size_t> IA;
    std::vector<size_t> JA;
    std::vector<T> AA;

    std::size_t ncol = 0;
    std::size_t nrow = 0;
};

/**
 * @brief Construct a new Matrix< T,  Storage Order>:: Matrix object
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param r 
 * @param c 
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(std::size_t const& r, size_t const& c) :
    nrow(r), ncol(c) {}

//! costructor starting from a full matrix
/**
 * @brief Construct a new Matrix< T,  Storage Order>:: Matrix object
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param v 
 * @param o 
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(const std::vector<std::vector<T>> &m)
{
    // loop over matrix (vector of vectors)
    // for (auto iti=v.cbegin(); iti!=v.cend(); ++iti)
    // {
    //     for (auto itj = iti->cbegin(); itj != iti->cend(); ++itj)
    //     {
            
    //     }
    // }

    ncol = m[0].size();
    nrow = m.size();
    for (std::size_t i=0; i<nrow; ++i)
    {
        //! error in output
        if (m[i].size() != ncol)
        {
            std::cerr << "number of elements in columns is not consistent" << std::endl;
            return;
        }

        //! initialize - loop over all elements
        //! could be done with swap? - no... threshold
        for (std::size_t j = 0; j<ncol; ++j)
        {
            if (std::abs(m[i][j]) > ZERO_TOL)
            {
                dynamic_data.insert( { {i,j}, m[i][j]} );
            }
        }
    }

}

/**
 * @brief Construct a new Matrix< T,  Storage Order>:: Matrix object
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param m 
 * @param o 
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(const std::vector<std::vector<T>> &m, Order const &o) :
    Matrix(m)
{
    ordering = o;
    switch (ordering)
    {
    case Order::COO:
        std::cout << "COO initialization" << std::endl;
        break;
    case Order::CSR:
        std::cout << "CSR initialization" << std::endl;
        break;
    case Order::CSC:
        std::cout << "CSC initialization" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR initialization" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC initialization" << std::endl;
        break;
    }
}

//! constructor reading from file
/**
 * @brief Construct a new Matrix< T,  Storage Order>:: Matrix object
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param name 
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(std::string const &name)
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

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param r 
 * @param c 
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::resize(std::size_t const& r, size_t const& c)
{
    if (compressed)
    {
        std::cerr << "cannot resize compressed matrix" << std::endl;
        return;
    }
    
    ncol = r;
    nrow = c;
}

//! is_compressed
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @return true 
 * @return false 
 */
template<typename T, typename StorageOrder>
bool Matrix<T, StorageOrder>::is_compressed() const
{
    if (compressed)
        return true;
    else
        return false; 
}

//! compress
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::compress(Order const &o)
{
    if (compressed)
    {
        std::cout << "Matrix is already compressed" << std::endl;
        return;
    }

    // compressed flag
    compressed = true;
    ordering = o;
    switch (ordering)
    {
    case Order::COO:
    {
        compressed = false;
        return;
    }
    case Order::CSR:
    {
        // keep track of row index change and element
        std::size_t temp_ind = dynamic_data.cbegin()->first[0];
        // first element row
        IA.push_back(temp_ind);

        for (auto it=dynamic_data.cbegin(); it!=dynamic_data.cend(); ++it)
        {
            //std::cout << it->first[0] << " " << std::next(it)->first[0];
            //* value
            AA.push_back(it->second);
            //* column
            JA.push_back(it->first[1]);
            
            // index of new row
            if(it->first[0] != temp_ind)
            {
                temp_ind = it->first[0];

                // add same element to empty rows
                auto it_tmp = std::prev(it);
                for( std::size_t k=it_tmp->first[0]+1; k<temp_ind; ++k )
                {
                    IA.push_back(AA.size()-1);
                }

                //* rows
                temp_ind = it->first[0];
                IA.push_back(AA.size()-1);
            }
        }
        IA.push_back(AA.size()-1);

        break;
    }

    case Order::CSC:
        break;

    case Order::MSC:
        break;

    case Order::MSR:
        break;
    } // switch(orderign)

    // clear memory
    dynamic_data.clear();
}


//! uncompress
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 */
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

    switch (ordering) 
    {
    case Order::COO:
        compressed = true;
        return;

    case Order::CSC:
        {
        // sizes of arrays
        std::size_t i_sz = IA.size();
        std::size_t j_sz = JA.size();

        // insert elements in map
        for (std::size_t i=0; i<i_sz; ++i)
        {
            // use I vector to loop from index i to i+1 in vector J and data
            for (std::size_t j=IA[i]; j<IA[i+1]; ++j)
            {
                // {col,     row}, data
                dynamic_data.insert( { {i, JA[j]}, AA[j]} );
            }
        }

        // insert last element
        dynamic_data.insert( { {nrow-1, JA[j_sz-1]}, AA[j_sz-1]} );
        }

    case Order::CSR:
        break;

    case Order::MSC:
        break;

    case Order::MSR:
        break;
    }

    // sizes of arrays
    std::size_t i_sz = IA.size();
    std::size_t j_sz = JA.size();

    // insert elements in map
    for (std::size_t i=0; i<i_sz; ++i)
    {
        // use I vector to loop from index i to i+1 in vector J and data
        for (std::size_t j=IA[i]; j<IA[i+1]; ++j)
        {
            // {col,     row}, data
            dynamic_data.insert( { {i, JA[j]}, AA[j]} );
        }
    }

    // insert last element
    dynamic_data.insert( { {nrow-1, JA[j_sz-1]}, AA[j_sz-1]} );
    AA.clear();
    JA.clear();
    IA.clear();

}


//! print
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 */
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
        std::size_t i_sz = IA.size();

        //* j index along element vector
        for(std::size_t j=0; j<JA.size(); )
        {

            // if index of element vector equal to next element, then new rob
            // index of element vector must be different from last element
            if ( (j == IA[i+1]) and (j!=IA[i_sz-1]) ) 
            {
                ++i;
            }

            // if next index is the same, then it is an empty row
            if (IA[i+1] == IA[i])
            {
                ++i;
                continue;
            }

            std::cout << i << " " << JA[j] << ": " << AA[j] << std::endl;

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

template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_one() const
{
    double res=0.0;
    std::vector<T> sums(ncol);
    // max of sum by columns

    switch (ordering)
    {
    case Order::COO:
    {
        std::cout << "COO norm-1" << std::endl;

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
    }
    case Order::CSR:
        std::cout << "CSR norm-1" << std::endl;
        break;
    case Order::CSC:
        std::cout << "CSC norm-1" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR norm-1" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC norm-1" << std::endl;
        break;
    }

    return res;
}

template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_infty() const
{
    double res=0.0;
    double sum=0.0;
    std::size_t id_tmp = dynamic_data.begin()->first[0];
    // max of sum by rows

    switch (ordering)
    {
    case Order::COO:
    {
        std::cout << "COO norm-infinity" << std::endl;

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
    }
    case Order::CSR:
        std::cout << "CSR norm-infinity" << std::endl;
        break;
    case Order::CSC:
        std::cout << "CSC norm-infinity" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR norm-infinity" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC norm-infinity" << std::endl;
        break;
    }

    return res;
}

template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_frob() const
{
    double res=0.0;

    switch (ordering)
    {
    case Order::COO:
    {
        std::cout << "COO norm-frobenius" << std::endl;

        // sum of all elements squared
        for(auto it = dynamic_data.cbegin(); it != dynamic_data.cend(); ++it)
        {
            res += std::abs(it->second) * std::abs(it->second);
        }

        break;
    }
    case Order::CSR:
        std::cout << "CSR norm-frobenius" << std::endl;
        break;
    case Order::CSC:
        std::cout << "CSC norm-frobenius" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR norm-frobenius" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC norm-frobenius" << std::endl;
        break;
    }

    return std::sqrt(res);
}

//! norm
// NOT working for both compressed and uncompressed
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param n 
 * @return double 
 */
template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm(Norm const &n) const
{
    double res = 0.;

    switch (n)
    {
    case Norm::One:
    {
        return norm_one();
    }
    case Norm::Infinity:
    {
        return norm_infty();
    }
    case Norm::Frobenius:
    {
        return norm_frob();
    }

    } // switch(n)

    return res;
}


// operator[] access copy
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param i 
 * @return T 
 */
template<typename T, typename StorageOrder>
T Matrix<T, StorageOrder>::operator[] (indexes const &i) const
{
    switch (ordering)
    {
    case Order::COO:
    {
        std::cout << "COO subscript copy" << std::endl;

        // check if present
        if (dynamic_data.count(i))
        {
            return dynamic_data.at(i);
        }

        // if out of bounds error
        if (i[0]>nrow or i[1]>ncol)
        {
            std::cerr << "out of bound index" << std::endl; 
        }
        
        return 0;
    }
    case Order::CSR:
    {
        std::cout << "CSR subscript copy" << std::endl;
        
        // first element of selected row
        std::size_t elem_first = IA[ i[0] ];
        return AA[ elem_first + i[1]];
    }
    case Order::CSC:
        std::cout << "CSC subscript copy" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR subscript copy" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC subscript copy" << std::endl;
        break;
    }

}

// operator[] access assign
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param i 
 * @return T& 
 */
template<typename T, typename StorageOrder>
T& Matrix<T, StorageOrder>::operator[] (indexes const &i)
{

    switch (ordering)
    {
    case Order::COO:
    {
        std::cout << "COO subscript reference" << std::endl;
        // cannot tell if access for assignment or copy for non-const Matrix

        // if out of bounds error
        if (i[0]>nrow or i[1]>ncol)
        {
            std::cerr << "out of bound index - assigning out of bounds" << std::endl;
        }

        // add new element if not present
        return dynamic_data[i];
    }
    case Order::CSR:
    {
        std::cout << "CSR subscript reference" << std::endl;

        // first element of selected row
        std::size_t elem_first = IA[ i[0] ];
        return AA[ elem_first + i[1]];
    }
    case Order::CSC:
        std::cout << "CSC subscript reference" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR subscript reference" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC subscript reference" << std::endl;
        break;
    }

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

    switch (m.ordering)
    {
    case Order::COO:
    {
        std::cout << "COO matrix-vector multiplication" << std::endl;

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
        break;
    }
    case Order::CSR:
    {
        std::cout << "CSR matrix-vector multiplication" << std::endl;

        // sizes of arrays
        std::size_t i_sz = m.IA.size();
        std::size_t j_sz = m.JA.size();

        // insert elements in map
        for (std::size_t i=0; i<i_sz; ++i)
        {
            std::size_t k = m.IA[i];
            // use I vector to loop from index i to i+1 in vector J and data
            for (std::size_t j=0; k+j<m.IA[i+1]; ++j)
            {
                // {col, row}
                res[i] += m.AA[k+j] * v[j];

                //std::cout << i << ", " << j
                //        << m.AA[k+j] << " " << v[j] << std::endl;
            }
        }
        res[ res.size()-1 ] += m.AA[j_sz-1] * v[v_sz-1];
        break;
    }
    case Order::CSC:
        std::cout << "CSC matrix-vector multiplication" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR matrix-vector multiplication" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC matrix-vector multiplication" << std::endl;
        break;
    }

    return res;
}

//! matrix-matrix multiplication
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam StorageOrder 
 * @param m1 
 * @param m2 
 * @return Matrix<T,StorageOrder> 
 */
template<typename T, typename StorageOrder>
Matrix<T,StorageOrder> operator*(Matrix<T,StorageOrder> const &m1, Matrix<T,StorageOrder> const &m2 )
{
    Matrix<T,StorageOrder> res;

    //! check all different orderings if equal for m1 and m2
    if (m1.ordering == m2.ordering)
    {

    }

    switch (m1.ordering)
    {
    case Order::COO:
        std::cout << "COO matrix-matrix multiplication" << std::endl;
        break;
    case Order::CSR:
        std::cout << "CSR matrix-matrix multiplication" << std::endl;
        break;
    case Order::CSC:
        std::cout << "CSC matrix-matrix multiplication" << std::endl;
        break;
    case Order::MSR:
        std::cout << "MSR matrix-matrix multiplication" << std::endl;
        break;
    case Order::MSC:
        std::cout << "MSC initialization" << std::endl;
        break;
    }

    return res;
}

} // namespace algebra

#endif