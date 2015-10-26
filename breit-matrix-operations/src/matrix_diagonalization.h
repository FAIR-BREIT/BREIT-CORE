/* 
 * File:   matrix_diagonalization.h
 * Author: winckler
 *
 * Created on July 15, 2015, 11:11 PM
 */

#ifndef MATRIX_DIAGONALIZATION_H
#define	MATRIX_DIAGONALIZATION_H

// std
#include <complex>


#if defined(__clang__)
_Pragma("clang diagnostic push") 
_Pragma("clang diagnostic ignored \"-Wshadow\"") 

// boost
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// bindings
#include "boost/numeric/bindings/lapack/syev.hpp"
#include "boost/numeric/bindings/lapack/gesvd.hpp"
#include "boost/numeric/bindings/lapack/gesdd.hpp"
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"
#include <boost/numeric/bindings/lapack/geev.hpp>
_Pragma("clang diagnostic pop")

#elif defined(__GNUC__) || defined(__GNUG__)    
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Wshadow\"")

// boost
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// bindings
#include "boost/numeric/bindings/lapack/syev.hpp"
#include "boost/numeric/bindings/lapack/gesvd.hpp"
#include "boost/numeric/bindings/lapack/gesdd.hpp"
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"
#include <boost/numeric/bindings/lapack/geev.hpp>

_Pragma("GCC diagnostic pop")
#endif    

// breit
#include "def.h"

namespace lapack = boost::numeric::bindings::lapack;
namespace breit
{
    
//namespace ublas  = boost::numeric::ublas;
//namespace lapack = boost::numeric::bindings::lapack;
    // general case of diagonalization
    template<typename T>
    inline int diagonalize_gen(
                            ublas::matrix<T, ublas::column_major>& A, 
                            ublas::vector<std::complex<T> >& eigen_values,
                            ublas::matrix<std::complex<T>, ublas::column_major>* eigen_vectors_inv,
                            ublas::matrix<std::complex<T>, ublas::column_major>* eigen_vectors
                          )
    {
        int i_err=lapack::geev(A, eigen_values, eigen_vectors_inv, eigen_vectors, lapack::optimal_workspace());
        return i_err;
    }
    
    // diagonalize symetric matrix : trivial case because eigenvalues are all real with orthogonal eigenmatrix
    template<typename M, typename V>
    inline int diagonalize_sym(M& eigen_vectors, V& eigen_values) 
    {
        int i_err = lapack::syev( 'V', 'U', eigen_vectors, eigen_values, lapack::minimal_workspace() );
        return i_err;
    }
    
} // breit namespace

#endif	/* MATRIX_DIAGONALIZATION_H */

