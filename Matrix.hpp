/**
 * @file
 *
 * @author Luca Brambilla <luca13.brambilla@mail.polimi.it>
 */

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
enum Order {Column_major, Row_major};
enum Norm {One, Infinity, Frobenius};
enum Compression {CSR, CSC};

// forward declaration matrix class
template <typename T, typename StorageOrder>
class Matrix;

// forward declaration of friend function inside class template
template<typename T, typename StorageOrder>
std::vector<T> operator*( Matrix<T, StorageOrder> const &m, std::vector<T> const &v );

template<typename T, typename StorageOrder>
Matrix<T,StorageOrder> operator*( Matrix<T,StorageOrder> const &m1, Matrix<T,StorageOrder> const &m2);

/**
 * @brief Template class for sparse matrices. Template parameters are the data type 
 * and internal ordering (column major or column major).
 * 
 * It is possible to pass from an uncompressed state (coordinate representation) to
 * a compressed state (compressed sparse row or compressed sparse column representation).
 *
 * Methods are available to compute the matrix norm and matrix-vector product. 
 * 
 * @tparam T                Data type
 * @tparam StorageOrder     Enumerator for storage ordering: column major or row major
 */
template <typename T, typename StorageOrder>
class Matrix
{
public:
    typedef std::array<std::size_t,2> indexes;
    typedef std::vector<std::vector<T>> uncompressed;
    typedef std::map<indexes,T> coo_matrix;

public:
    // constructors
    Matrix();

    Matrix(std::size_t const& r, size_t const& c);

    Matrix(uncompressed const &m, Order const &o=Row_major);
    
    Matrix(Matrix const &m);
    
    Matrix(std::string const &name);
    Matrix(std::string const &name, Order const &o=Row_major);

    // getters
    std::size_t ncols() { return ncol; };
    std::size_t nrows() { return nrow; };

    // utilities
    void resize(std::size_t const& r, size_t const& c);
    void print() const;

    // compression utilities
    void compress(Compression const &c);
    void uncompress();
    bool is_compressed() const;

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
    Order ordering = Order::Column_major;
    Compression compression = Compression::CSR;

    bool compressed = false;

    coo_matrix dynamic_data;

    std::vector<size_t> IA;
    std::vector<size_t> JA;
    std::vector<T> AA;

    std::size_t ncol = 0;
    std::size_t nrow = 0;
};

/**
 * @brief Construct a new empty Matrix object assigning the shape
 * 
 * @param r         number of rows
 * @param c         number of columns
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(std::size_t const& r, size_t const& c) :
    nrow(r), ncol(c) {}



/**
 * @brief Construct a new Matrix starting from a full matrix (vector of vector
 * representation) with a given ordering (default is row-major).
 * 
 * @param m         input matrix
 * @param o         ordering
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(const std::vector<std::vector<T>> &m,
                                Order const &o)
{

    ncol = m[0].size();
    nrow = m.size();

    switch (ordering){
    case Row_major:
    {
        for (std::size_t i=0; i<nrow; ++i)
        {
            // error in output
            if (m[i].size() != ncol)
            {
                std::cerr << "number of elements in columns is not consistent" << std::endl;
                return;
            }

            // initialize - loop over all elements
            for (std::size_t j = 0; j<ncol; ++j)
            {
                // rows first, columns second
                if (std::abs(m[i][j]) > ZERO_TOL)
                {
                    dynamic_data.insert( { {i,j}, m[i][j]} );
                }
            }
        }
        break;
    }

    case Column_major:
    {
        for (std::size_t i=0; i<nrow; ++i)
        {
            // error in output
            if (m[i].size() != ncol)
            {
                std::cerr << "number of elements in columns is not consistent" << std::endl;
                return;
            }

            // initialize - loop over all elements
            for (std::size_t j = 0; j<ncol; ++j)
            {
                // columns first, rows second
                if (std::abs(m[j][i]) > ZERO_TOL)
                {
                    dynamic_data.insert( { {j,i}, m[j][i]} );
                }
            }
        }
        break;
    }
    
    } // switch(ordering)

}


/**
 * @brief Construct a new Matrix object reading from a file.
 * 
 * @param name        String containing the path to the file to read.
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
 * @brief Construct a new Matrix object starting from an existing Matrix.
 * 
 * @param m             Matrix object
 */
template<typename T, typename StorageOrder>
Matrix<T, StorageOrder>::Matrix(Matrix const &m) :
    ordering(m.ordering), compression(m.compression), compressed(m.compressed),
    dynamic_data(m.dynamic_data), IA(m.IA), JA(m.JA), AA(m.AA),
    ncol(m.ncol), nrow(m.nrow)
{}


/**
 * @brief Reshape the matrix passing the new number of rows and columns
 * 
 * @param r         new number of rows
 * @param c         new number of columns
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::resize(std::size_t const& r, size_t const& c)
{
    if (compressed)
    {
        std::cerr << "cannot resize compressed matrix" << std::endl;
        return;
    }
    
    //! delete elements if reshape with smaller size
    if ( (r<nrow) or (c<ncol) )
    {

    }

    ncol = r;
    nrow = c;

}

/**
 * @brief Checks if the matrix is compressed and returns true if compressed, false 
 * otherwise.
 * 
 * @return bool
 */
template<typename T, typename StorageOrder>
bool Matrix<T, StorageOrder>::is_compressed() const
{
    if (compressed)
        return true;
    else
        return false; 
}

/**
 * @brief Pass from a coordinate representation to a compressed representation.
 * Possible representations are:
 * - Compressed Sparse Row (CSR)
 * - Compressed Sparse Column (CSC) 
 * 
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::compress(Compression const &c)
{
    if (compressed)
    {
        std::cout << "Matrix is already compressed" << std::endl;
        return;
    }


    switch (c)
    {
    case Compression::CSR:
    {
        if (ordering != Order::Row_major)
        {
            std::cerr << "only compress to CSR if row-major ordering" << std::endl;
            return;
        }

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

    case Compression::CSC:
    {
        if (ordering != Order::Column_major)
        {
            std::cerr << "only compress to CSC if column-major ordering" << std::endl;
            return;
        }


        break;
    }

    } // switch(compression)

    // compressed flag
    compressed = true;
    compression = c;
    // clear memory
    dynamic_data.clear();
}


/**
 * @brief Pass from a compressed representation to the coordinate representation.
 * 
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::uncompress()
{
    if (!compressed)
    {
        std::cout << "Matrix is already uncompressed" << std::endl;
        return;
    }

    switch (compression) 
    {
    case Compression::CSR:
        {

        if (ordering != Order::Row_major)
        {
            std::cerr << "only decompress from CSR if column-major ordering" << std::endl;
            return;
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
        }

    case Compression::CSC:
    {
        if (ordering != Order::Column_major)
        {
            std::cerr << "only decompress from CSC if column-major ordering" << std::endl;
            return;
        }

        break;
    }

    } //switch(compression)

    AA.clear();
    JA.clear();
    IA.clear();

}


/**
 * @brief Print in coordinate representation
 * 
 * @tparam T 
 * @tparam StorageOrder 
 */
template<typename T, typename StorageOrder>
void Matrix<T, StorageOrder>::print() const
{
    // print if not compressed
    if (!compressed)
    {
        std::cout << "print uncompressed" << std::endl;
        for (auto it=dynamic_data.cbegin(); it!=dynamic_data.cend(); ++it)
        {
            std::cout << it->first[0] << " " << it->first[1] << ": " << it->second << std::endl;
        }
    }

    // print if compressed format
    std::cout << "print compressed" << std::endl;

    switch(compression)
    {
    case Compression::CSR:
    {
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
    case Compression::CSR:
    {
        return;
    }

    } //switch(ordering)



}


/**
 * @brief Compute the 1-norm of the matrix.
 * 
 * @return double 
 */
template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_one() const
{
    double res=0.0;
    std::vector<T> sums(ncol);
    // max of sum by columns


    if (!compressed)
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
    }


    switch (compression)
    {
    case Compression::CSR:
    {
        std::cout << "CSR norm-1" << std::endl;
        break;
    }
    case Compression::CSC:
    {
        std::cout << "CSC norm-1" << std::endl;
        break;
    }

    } //switch(compression)

    return res;
}


/**
 * @brief Compute the infinity norm of the matrix.
 * 
 * @return double 
 */
template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_infty() const
{
    double res=0.0;
    double sum=0.0;
    std::size_t id_tmp = dynamic_data.begin()->first[0];
    // max of sum by rows

    if (!compressed)
    {
        std::cout << "COO infinity norm" << std::endl;

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
    }


    switch(compression)
    {
    case Compression::CSR:
    {
        std::cout << "CSR infinity norm" << std::endl;
        break;
    }
    case Compression::CSC:
    {
        std::cout << "CSC infinity norm" << std::endl;
        break;
    }

    } //switch(compression)

    return res;
}


/**
 * @brief Compute the Frobenius norm of the matrix.
 * 
 * @return double 
 */
template<typename T, typename StorageOrder>
double Matrix<T, StorageOrder>::norm_frob() const
{
    double res=0.0;

    if (!compressed)
    {
        std::cout << "COO Frobenius norm" << std::endl;

        // sum of all elements squared
        for(auto it = dynamic_data.cbegin(); it != dynamic_data.cend(); ++it)
        {
            res += std::abs(it->second) * std::abs(it->second);
        }

    }
    switch (ordering) {
    
    case Compression::CSR:
        std::cout << "CSR Frobenius norm" << std::endl;
        break;
    case Compression::CSC:
        std::cout << "CSC Frobenius norm" << std::endl;
        break;

    } // switch(ordering)

    return std::sqrt(res);
}


//! NOT working for both compressed and uncompressed
/**
 * @brief Given an enumerator, return the desired matrix norm.
 * Possible norms are:
 * - 1-norm
 * - Infinity norm
 * - Frobenius norm
 * 
 * @param n             enum
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
 * @brief Const version of subscript operator. If a
 * 
 * @param i 
 * @return T 
 */
template<typename T, typename StorageOrder>
T Matrix<T, StorageOrder>::operator[] (indexes const &i) const
{
    if (!compressed)
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


    switch (compression)
    {
    case Compression::CSR:
    {
        std::cout << "CSR subscript copy" << std::endl;
        
        // first element of selected row
        std::size_t elem_first = IA[ i[0] ];
        return AA[ elem_first + i[1]];
    }
    case Compression::CSC:
    {
        std::cout << "CSC subscript copy" << std::endl;
        break;
    }

    } //switch(compression)

}

// operator[] access assign
/**
 * @brief 
 * 
 * @param i 
 * @return T& 
 */
template<typename T, typename StorageOrder>
T& Matrix<T, StorageOrder>::operator[] (indexes const &i)
{

    if (!compression)
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


    switch (ordering)
    {
    
    case Compression::CSR:
    {
        std::cout << "CSR subscript reference" << std::endl;

        // first element of selected row
        std::size_t elem_first = IA[ i[0] ];
        return AA[ elem_first + i[1]];
    }
    case Compression::CSC:
    {
        std::cout << "CSC subscript reference" << std::endl;
        break;
    }

    } // switch(ordering)

}

/**
 * @brief Matrix-vector multiplication.
 * 
 * @param m             Matrix object
 * @param v             Standard vector
 * @return std::vector<T> 
 */
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

    if (!m.compressed)
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
    }

    switch (m.compression) {
    

    case Compression::CSR:
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
    case Compression::CSC:
        std::cout << "CSC matrix-vector multiplication" << std::endl;
        break;
    }

    return res;
}


/**
 * @brief Matrix-Matrix multiplication. Only works if objects have the same ordering
 * and the same compression.
 * 
 * @param m1            First matrix object
 * @param m2            Second matrix object
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

    if (!m1.compressed)
    {

    }

    switch (m1.compression)
    {
    case Compression::CSR:
        std::cout << "CSR matrix-matrix multiplication" << std::endl;
        break;
    case Compression::CSC:
        std::cout << "CSC matrix-matrix multiplication" << std::endl;
        break;
    }

    return res;
}

} // namespace algebra

#endif