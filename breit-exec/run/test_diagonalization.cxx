/* 
 * File:   test_diagonalization.cxx
 * Author: winckler
 *
 * Created on July 15, 2015, 10:57 PM
 */

#include "equations_manager.h"
#include "breit_equations.h"
#include "solve_breit_equations.h"
#include "matrix_diagonalization.h"

#include <boost/numeric/ublas/storage.hpp>
using namespace breit;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack; 

typedef breit_equations<double>                                  equations_d;
typedef solve_breit_equations<double>                            solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d>    breit_manager;
typedef ublas::matrix<double,ublas::column_major>               matrix_d;
typedef ublas::matrix<std::complex<double>,ublas::column_major> matrix_c;
typedef ublas::vector<std::complex<double> >                    vector_c;
typedef ublas::vector<double >                                  vector_d;
/*
 * Assumptions : the matrix A of the system can be diagonalized in R and eigen values are non-zero.
 * (recall : assumption satisfied if A real symetric)
 * 
 * Then :
 * 
 * dF/dx = A * F <=> dF/dx = P * D * V      (with V=P^(-1) )
 *               <=> dY/dx = D * Y          (with Y = V * F)
 *               <=> solution of the form Y = ( c1 exp(-lambda1 * x) , ... , cN exp(-lambdaN * x) )
 * 
 * boundary for F(x=0) gives the unknown coefficients : (c1, ..., cN) = V * F(x=0)
 * 
 */
enum class diagonalisation
{
    unknown,
    in_R,
    in_C,
    no
};

int main(int argc, char** argv) 
{
    enum diagonalisation diag_case = diagonalisation::in_C;
    
    solve_breit_equations<double> solver;

    matrix_d A(3,3);
    matrix_d B(3,3);
    matrix_c P(3,3);
    matrix_c P_inv(3,3);
    
    vector_c eigen_values(3);
    vector_d snd_member(3);
    
    switch (diag_case)
    {
        case diagonalisation::in_R :
            LOG(INFO)<<"Define matrix diagonalizable in R";
            A(0,0)=1; A(0,1)=2;  A(0,2)=0;
            A(1,0)=0; A(1,1)=3;  A(2,2)=0;
            A(2,0)=2; A(2,1)=-4; A(2,2)=2;
            break;
            
        case diagonalisation::in_C :
            LOG(INFO)<<"Define matrix diagonalizable in C";
            //*
            A(0,0)=1;  A(0,1)=1;  A(0,2)=0;
            A(1,0)=-1; A(1,1)=2;  A(1,2)=1;
            A(2,0)=1;  A(2,1)=0;  A(2,2)=1;
            // */
            /*
            A(0,0)=58;  A(0,1)=9;  A(0,2)=2;
            A(1,0)=186; A(1,1)=383;  A(1,2)=96;
            A(2,0)=-912;  A(2,1)=-1551;  A(2,2)=-388;
            */
            break;
            
        case diagonalisation::no :
            LOG(INFO)<<"Define matrix not diagonalizable";
            A(0,0)=2; A(0,1)=0;  A(0,2)=0;
            A(1,0)=0; A(1,1)=3;  A(2,2)=1;
            A(2,0)=0; A(2,1)=0;  A(2,2)=3;
            break;
            
        case diagonalisation::unknown :
            break;
    }
    
            
    //int diag_gen_err=diagonalize_gen(A,eigen_values,&P_inv,&P);
    //if(diag_gen_err)
    //    LOG(DEBUG)<<"diagonalize_gen function returned error value "<<diag_gen_err;
    
    //solver.solve_homogeneous_system(A,snd_member);
    vector_d ini_cond(3);
    solver.solve(A,snd_member,ini_cond);
    
    
    
    
    // */
    
    return 0;
}


