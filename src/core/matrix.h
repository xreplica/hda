/*
 * Matrix class header:
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/matrix.h $
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <string>
#include <cstring>

using namespace std;

class matrixException:public exception
{
private:
    const string msg;

public:
    matrixException(string details = "") : msg(details) {}
    ~matrixException() throw() {};

    virtual const char* what() const throw()
    {
        if(msg.length() > 0)
            return ("Illegal operation on a matrix: " + msg).c_str();
        else
            return "Illegal operation on a matrix";
    }
};

class Matrix
{
private:
    typedef vector<double>::size_type st;

    unsigned long       _rows, _cols;
    std::vector<double> _coeffs;
    st                  indexof(int i, int j) const {return st(i)*_cols+j;} // row major storage

public:
    Matrix();
    Matrix(unsigned long);
    Matrix(unsigned long, unsigned long);
    Matrix(vector<double>);

    unsigned long rows() const {return _rows;}
    unsigned long cols() const {return _cols;}

    static Matrix identity(unsigned long);

    void          reset(unsigned long, unsigned long);
    void          reset(unsigned long);
    void          swap(Matrix&);
    inline void   swap(Matrix&, Matrix&);
    Matrix        transposed() const;
    void          transpose();
    Matrix        inverse() const;
    void          invert();
    double        minorMatrix(unsigned long, unsigned long) const;
    double        cofactor(unsigned long, unsigned long) const;
    Matrix        cofactorsMatrix() const;
    Matrix        adjoint() const;
    Matrix        adjugate() const;
    double        determinant()const;
    Matrix        cholesky() const;

    inline double const& operator()(int i, int j) const
    {
        if((unsigned long)i>=_rows || (unsigned long)j>= _cols)
            throw(matrixException("index out of bounds"));
        return _coeffs[indexof(i,j)];
    }

    inline double& operator()(int i, int j)
    {
        if((unsigned long)i>=_rows || (unsigned long)j>= _cols)
            throw(matrixException("index out of bounds"));
        return _coeffs[indexof(i,j)];
    }

    inline double const* operator[](int i) const
    {
        if((unsigned long)i>=_rows)
            throw(matrixException("index out of bounds"));
        return &_coeffs[indexof(i,0)];
    }

    inline double* operator[](int i)
    {
        if((unsigned long)i>=_rows)
            throw(matrixException("index out of bounds"));
        return &_coeffs[indexof(i,0)];
    }

    Matrix  operator+(const Matrix&) const;
    Matrix& operator+=(const Matrix&);
    Matrix  operator-(const Matrix&) const;
    Matrix& operator-=(const Matrix&);
    Matrix  operator*(const Matrix&) const;
    Matrix  operator*(const double&) const;
    Matrix  operator*(const int&) const;
    Matrix& operator*=(const Matrix&);
    Matrix& operator*=(const double&);
    Matrix& operator*=(const int&);
    Matrix  operator/(const Matrix&) const;
    Matrix  operator/(const double&) const;
    Matrix  operator/(const int&) const;
    Matrix& operator/=(const Matrix&);
    Matrix& operator/=(const double&);
    Matrix& operator/=(const int&);
    Matrix  operator^(unsigned int const& rhs) const;
};

#endif
