/*
 * Matrix class implementation
 *
 * $Author: francois $
 * $Date: 2015-08-18 20:31:22 -0400 (Tue, 18 Aug 2015) $
 * $Revision: 38 $
 * $URL: https://127.0.0.1:10000/svn/hda/trunk/src/core/matrix.cpp $
 */

#ifndef __MATRIX_CPP__
#define __MATRIX_CPP__

#include "src/core/matrix.h"

//CONSTRUCTORS

Matrix::Matrix() : _rows(1), _cols(1) , _coeffs(1, 0.0){}

Matrix::Matrix(unsigned long d) : _rows(d), _cols(d), _coeffs(st(d)*d, 0.0) {}

Matrix::Matrix(unsigned long r, unsigned long c) : _rows(r), _cols(c), _coeffs(st(r)*c, 0.0) {}

/*Matrix::Matrix(const Matrix& m) : _rows(m._rows), _cols(m._cols), _coeffs(m._coeffs) {}

Matrix::Matrix(Matrix&& m) : _rows(move(m._rows)), _cols(move(m._cols)), _coeffs(move(m._coeffs)) {}*/

Matrix::Matrix(vector<double> vect) : _rows(vect.size()), _cols(1), _coeffs(vect.size(), 0.0)
{
    for(unsigned long i=0; i<vect.size(); ++i)
        _coeffs.at(i) = vect.at(i);
}

//FUNCTIONS

Matrix Matrix::identity(unsigned long s)
{
    Matrix identity(s, s);
    for(unsigned long i=0; i<s; ++i)
        identity[i][i] = 1.0;

    return identity;
}

void Matrix::reset(unsigned long r, unsigned long c)
{
    _rows = r;
    _cols = c;
    _coeffs.clear();
    _coeffs.resize(st(r)*c, 0.0);
}

void Matrix::reset(unsigned long s)
{
    _rows = s;
    _cols = s;
    _coeffs.clear();
    _coeffs.resize(st(s)*s, 0.0);
}

void Matrix::swap(Matrix& that)
{
    std::swap(this->_rows,that._rows);
    std::swap(this->_cols,that._cols);
    this->_coeffs.swap(that._coeffs);
}

inline void Matrix::swap(Matrix& a, Matrix& b)
{
    a.swap(b);
}

Matrix Matrix::transposed() const
{
    Matrix tmp(_cols, _rows);
    for(unsigned long i=0; i< _rows; ++i)
        for(unsigned long j=0; j<_cols; ++j)
            tmp[j][i] = _coeffs[indexof(i,j)];
    return tmp;
}

void Matrix::transpose()
{
    Matrix tmp = transposed();
    swap(tmp);
}

double Matrix::minorMatrix(unsigned long i, unsigned long j) const
{
    Matrix tmp(_rows-1, _cols-1);
    unsigned long i2 = 0;
    for(unsigned long i1=0; i1<_rows; ++i1)
    {
        if(i1 == i)
            continue;
        unsigned long j2 = 0;
        for(unsigned long j1=0; j1<_cols; ++j1)
        {
            if(j1 == j)
                continue;
            tmp[i2][j2] = _coeffs[indexof(i1,j1)];
            ++j2;
        }
        ++i2;
    }
    return tmp.determinant();
}

double Matrix::cofactor(unsigned long i, unsigned long j) const
{
    return ((i+j%2==0) ? minorMatrix(i, j) : minorMatrix(i, j)*-1);
}

Matrix Matrix::cofactorsMatrix() const
{
    if(_rows != _cols)
        throw(matrixException("cannot calculate the matrix of cofactors of non-square Matrix"));

    Matrix tmp(_rows, _cols);
    for(unsigned long i=0; i<_rows; ++i)
        for(unsigned long j=0; j<_cols; ++j)
            tmp[i][j] = cofactor(i, j);

    return tmp;
}

Matrix Matrix::adjoint() const
{
    return cofactorsMatrix().transposed();
}

Matrix Matrix::adjugate() const
{
    return cofactorsMatrix().transposed();
}

double Matrix::determinant() const
{
    double det = 0.0;
    if(_rows != _cols)
        throw(matrixException("cannot calculate determinant of non-square Matrix"));
    else if(_rows == 1)
        det = _coeffs[indexof(0,0)];
    else if(_rows == 2)
        det = (_coeffs[indexof(0,0)] * _coeffs[indexof(1,1)]) - (_coeffs[indexof(1,0)] * _coeffs[indexof(0,1)]);
    else
    {
        //only works for symetrical positive-definite matrices
        Matrix L = cholesky();
        det = 1.0;
        for(unsigned long i=0; i<_rows; ++i)
            det *= (L[i][i] * L[i][i]);
        /*for(unsigned long j=0; j<_cols; ++j)
            det += _coeffs[indexof(0,j)] * cofactor(0,j);
        */
    }

    return det;
}

Matrix Matrix::inverse() const
{
    if(_rows != _cols)
        throw(matrixException("cannot calculate inverse of non-square Matrix"));
    if(_rows < 2)
        throw(matrixException("cannot calculate inverse of a matrix with dimension less than 2"));

    for(unsigned long i=0; i<_rows; ++i)
        if(_coeffs[indexof(i,i)] == 0)
            throw(matrixException("cannot calculate inverse of a matrix with 0 on its diagonal"));

    /*double det = determinant();
    if(det == 0.0)
        throw(matrixException("cannot calculate inverse of a Matrix where determinant = 0"));*/

    Matrix tmp(_rows,_cols);
    if(_rows == 2)
    {
        tmp[0][0] = _coeffs[indexof(1,1)];
        tmp[0][1] = _coeffs[indexof(0,1)] * -1;
        tmp[1][0] = _coeffs[indexof(0,0)];
        tmp[1][1] = _coeffs[indexof(1,1)] * -1;
        tmp *= (1.0/abs(determinant()));
    }
    else if(_rows == 3)
    {
        for(unsigned long i=0; i<_rows; ++i)
        {
            for(unsigned long j=0; j<_cols; ++j)
            {
                tmp[i][j] = ((i+j%2==0) ? minorMatrix(i, j) : minorMatrix(i, j)*-1);
//              if((i+j)%2 == 0)
//                  tmp[i][j] = minorMatrix(i,j);
//              else
//                  tmp[i][j] = minorMatrix(i, j) * -1;
            }
        }
        tmp.transpose();
        tmp /= determinant();
    }
    else
    {
        //somehow still faulty
        //tmp = adjugate()*(1 / det);

        /*only works for symetrical positive definite matrices, but gives correct results --not anymore?*/
        /*
        Matrix L = cholesky();
        double Ldet = L.determinant();
        Matrix Ladj = L.adjugate();
        Matrix Linv = Ladj/Ldet;
        Matrix Lt = L.transposed();
        double Ltdet = Lt.determinant();
        Matrix Ltadj = Lt.adjugate();
        Matrix Ltinv = Ltadj/Ltdet;
        tmp = Ltinv * Linv;
        */
        //============================================================
        for(unsigned long i=0; i<_rows; ++i)
            for(unsigned long j=0; j<_rows; ++j)
                tmp[i][j] = _coeffs[indexof(i,j)];

        for (unsigned long i=1; i < _rows; ++i)
            tmp[0][i] /= tmp[0][0]; // normalize row 0

        for (unsigned long i=1; i < _rows; i++)
        {
            for (unsigned long j=i; j < _rows; j++)
            { // do a column of L
                double sum = 0.0;
                for (unsigned long k = 0; k < i; k++)
                    sum += tmp[j][k] * tmp[k][i];
                tmp[j][i] -= sum;
            }
            if (i == _rows-1)
                continue;
            for (unsigned long j=i+1; j < _rows; j++)
            {  // do a row of U
                double sum = 0.0;
                for (unsigned long k = 0; k < i; k++)
                    sum += tmp[i][k]*tmp[k][j];
                tmp[i][j] = (tmp[i][j]-sum) / tmp[i][i];
            }
        }

        for (unsigned long i = 0; i < _rows; i++ )
        { // invert L
            for (unsigned long j = i; j < _rows; j++ )
            {
                double x = 1.0;
                if ( i != j )
                {
                    x = 0.0;
                    for (unsigned long k = i; k < j; k++ )
                        x -= tmp[j][k]*tmp[k][i];
                }
                tmp[j][i] = x / tmp[j][j];
            }
        }

        for (unsigned long i = 0; i < _rows; i++ )
        { // invert U
            for (unsigned long j = i; j < _rows; j++ )
            {
                if ( i == j )
                    continue;
                double sum = 0.0;
                for (unsigned long k = i; k < j; k++ )
                    sum += tmp[k][j]*( (i==k) ? 1.0 : tmp[i][k] );
                tmp[i][j] = -sum;
            }
        }

        for (unsigned long i = 0; i < _rows; i++ )
        {   // final inversion
            for (unsigned long j = 0; j < _rows; j++ )
            {
                double sum = 0.0;
                for (unsigned long k = ((i>j)?i:j); k < _rows; k++ )
                    sum += ((j==k)?1.0:tmp[j][k])*tmp[k][i];
                tmp[j][i] = sum;
            }
        }

        //==========================================================
    }

    return tmp;
}

void Matrix::invert() {
    Matrix tmp = inverse();
    swap(tmp);
}

Matrix Matrix::cholesky() const
{
    if(_rows != _cols)
        throw(matrixException("cannot calculate Cholesky decomposition of non-square Matrix"));
    if(_rows < 2)
        throw(matrixException("cannot caluclate Cholesky decomposition of a matrix with dimension less than 2"));

    /*double det = determinant();
    if(det == 0)
        throw(matrixException("cannot calculate Cholesky decomposition of a Matrix where determinant = 0"));*/

    //CHECK FOR POSITIVE DEFINITE-NESS

    /*int n = _rows;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    Matrix l(n);
    l[0][0] = sqrt(_coeffs[indexof(0,0)]);
    for (int j = 1; j <= n-1; ++j)
        l[j][0] = _coeffs[indexof(j,0)] / l[0][0];
    for (int i = 1; i <= (n-2); ++i)
    {
        for (int k = 0; k <= (i-1); ++k)
            sum1 += pow(l[i][k], 2);
        l[i][i]= sqrt(_coeffs[indexof(i,i)] - sum1);
        for (int j = (i+1); j <= (n-1); ++j)
        {
            for (int k = 0; k <= (i-1); ++k)
                sum2 += l[j][k]*l[i][k];
                l[j][i]= (_coeffs[indexof(j,i)] - sum2) / l[i][i];
        }
    }
    for (int k = 0; k <= (n-2); ++k)
        sum3 += pow(l[n-1][k], 2);
    l[n-1][n-1] = sqrt(_coeffs[indexof(n-1,n-1)] - sum3);

    return l;*/

    Matrix L(_rows);

    for (unsigned long i = 0; i < _rows; ++i)
    {
        for (unsigned long j = 0; j < (i+1); ++j)
        {
            double s = 0;
            for (unsigned long k = 0; k < j; ++k)
                s += L[i][k] * L[j][k];
            L[i][j] = (i == j) ? sqrt(_coeffs[indexof(i,i)] - s) : (1.0 / L[j][j] * (_coeffs[indexof(i,j)] - s));
        }
    }

    return L;
}

//OPERATORS

Matrix Matrix::operator+(const Matrix& rhs) const {
    if(_cols!=rhs.cols()||_rows!=rhs.rows())
        throw(matrixException("Cannot add matrices of different dimensions"));

    Matrix tmp(_rows, _cols);
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            tmp[i][j] = _coeffs[indexof(i,j)] + rhs[i][j];

    return tmp;
}

Matrix& Matrix::operator+=(const Matrix& rhs)
{
    Matrix tmp = (*this);
    tmp = tmp + rhs;
    swap(tmp);

    return *this;
}

Matrix Matrix::operator-(Matrix const& rhs) const
{
    if(_cols!=rhs.cols()||_rows!=rhs.rows())
        throw(matrixException("Cannot subtract matrices of different dimensions"));

    Matrix tmp(_rows, _cols);
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            tmp[i][j] = _coeffs[indexof(i,j)] - rhs[i][j];

    return tmp;
}

Matrix& Matrix::operator-=(Matrix const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp - rhs;
    swap(tmp);

    return *this;
}

Matrix Matrix::operator*(Matrix const& rhs) const
{
    if(_cols!=rhs.rows())
        throw(matrixException("cannot multiply matrices of unmatching dimensions"));

    Matrix tmp(_rows,rhs.cols());
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            for(unsigned long k=0; k<_cols; ++k)
                tmp[i][j] += _coeffs[indexof(i,k)] * rhs[k][j];

    return tmp;
}

Matrix Matrix::operator*(double const& rhs) const
{
    Matrix tmp(_rows,_cols);
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            tmp[i][j] = _coeffs[indexof(i,j)] * rhs;

    return tmp;
}

Matrix Matrix::operator*(int const& rhs) const
{
    Matrix tmp = *this;

    return tmp * (double)rhs;
}

Matrix& Matrix::operator*=(Matrix const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp * rhs;
    swap(tmp);

    return *this;
}

Matrix& Matrix::operator*=(double const& rhs) {
    Matrix tmp = *this;
    tmp = tmp * rhs;
    swap(tmp);

    return *this;
}

Matrix& Matrix::operator*=(int const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp * rhs;
    swap(tmp);

    return *this;
}

Matrix Matrix::operator/(Matrix const& rhs) const
{
    if(_cols!=rhs.rows())
        throw(matrixException("cannot divide matrices of unmatching dimensions"));

    Matrix tmp(_rows,rhs.cols());
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            for(unsigned long k=0; k<_cols; ++k)
                tmp[i][j] += _coeffs[indexof(i,k)] / rhs[k][j];

    return tmp;
}

Matrix Matrix::operator/(double const& rhs) const
{
    Matrix tmp(_rows,_cols);
    for(unsigned long i=0; i<tmp.rows(); ++i)
        for(unsigned long j=0; j<tmp.cols(); ++j)
            tmp[i][j] = _coeffs[indexof(i,j)] / rhs;

    return tmp;
}

Matrix Matrix::operator/(int const& rhs) const
{
    Matrix tmp = *this;

    return tmp / (double)rhs;
}

Matrix& Matrix::operator/=(Matrix const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp / rhs;
    swap(tmp);

    return *this;
}

Matrix& Matrix::operator/=(double const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp / rhs;
    swap(tmp);

    return *this;
}

Matrix& Matrix::operator/=(int const& rhs)
{
    Matrix tmp = *this;
    tmp = tmp / rhs;
    swap(tmp);

    return *this;
}

Matrix Matrix::operator^(unsigned int const& rhs) const
{
    Matrix tmp = *this;
    for(unsigned int i=0; i<rhs; ++i)
        tmp*=tmp;

    return tmp;
}

#endif
