/* 
 * File:   def.h
 * Author: winckler
 *
 * Created on August 7, 2015, 11:05 AM
 */

#ifndef DEF_H
#define	DEF_H
    


#if defined(__clang__)

_Pragma("clang diagnostic push") 
_Pragma("clang diagnostic ignored \"-Wshadow\"") 

// boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

_Pragma("clang diagnostic pop")
#elif defined(__GNUC__) || defined(__GNUG__)    
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Wshadow\"")

// boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

_Pragma("GCC diagnostic pop")
#endif

#include <map>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace ublas  = boost::numeric::ublas;


class breit_summary
{
  public :
    breit_summary() :    filename(), 
                        outfilename(),  
                        title(), 
                        F_index_map(), 
                        approximated_solutions(), 
                        equilibrium_solutions(), 
                        distance_to_equilibrium(),
                        analytical_solutions(),
                        max_fraction_index(0),
                        system_dim(0), 
                        reduced_system_dim(0) , 
                        offset(0)
    {}
    virtual ~breit_summary (){}


    std::string filename;
    std::string outfilename;
    std::string title;
    std::map<size_t,int> F_index_map;// matrix index -> real index

    std::map<size_t,double> approximated_solutions;// matrix index -> Fi values
    std::map<size_t,double> equilibrium_solutions;// matrix index -> Fi values
    std::map<size_t,double> distance_to_equilibrium;
    std::map<size_t,std::string> analytical_solutions;
    
    size_t max_fraction_index;

    std::size_t system_dim;
    std::size_t reduced_system_dim;
    std::size_t offset;

};

namespace breit
{
    
    /*
    template<>
    default_sink& operator << < ublas::matrix<double,ublas::column_major> >
                ( default_sink &sink, const ublas::matrix<double,ublas::column_major> &data ) 
    {
        if (sink.fOutput_file.is_open())
            sink.fOutput_file << data;
        return sink;
    }
    
    template<>
    default_sink& operator << < ublas::matrix<std::complex<double>,ublas::column_major> >
                ( default_sink &sink, const ublas::matrix<std::complex<double>,ublas::column_major> &data ) 
    {
        if (sink.fOutput_file.is_open())
            sink.fOutput_file << data;
        return sink;
    }
    */
    std::ostream& operator<<(std::ostream& os, const ublas::matrix<double,ublas::column_major> &mat)
    {
        for(std::size_t i(0); i<mat.size1(); i++)
        {
            os<<"Row "<<i+1<<" :  ";
            for(unsigned int j=0; j<mat.size2(); j++)
            {
                // P cells = a 
                os<< "( "<< std::setw(12) << mat(i,j)<<" )"<<std::setw(6);
            }
            os<<"\n";
        }
        
        return os;
    }
    
    std::ostream& operator << ( std::ostream& os, const ublas::matrix<std::complex<double>,ublas::column_major> &mat ) 
    {
        for(std::size_t i(0); i<mat.size1(); i++)
        {
            os<<"Row "<<i+1<<" :  ";
            for(unsigned int j=0; j<mat.size2(); j++)
            {
                // P cells = a + ib 
                os << "( "  << std::setw(12) << mat(i,j).real() 
                            <<" , "
                            << std::setw(12) << mat(i,j).imag()
                            <<" )"<<std::setw(6);
            }
            os<<"\n";
        }
        
        return os;
    }
    
    
    
    
    
    
    
    
    

    template<typename T>
    class breit_matrix : public ublas::matrix<std::complex<T>,ublas::column_major>
    {
    public:
        breit_matrix() : ublas::matrix<std::complex<T>,ublas::column_major>(), fName("matrix")
        {
            
        }
        
        breit_matrix(size_t dim1, size_t dim2, const std::string& name="matrix") : 
            ublas::matrix<std::complex<T>,ublas::column_major>(dim1,dim2), fName(name)
        {
            
        }
        
        virtual ~breit_matrix()
        {
            
        }
    private:
        std::string fName;
    };

    typedef po::variables_map                                    variables_map;
    //typedef ublas::vector<double>                                        vector_d;
    //typedef ublas::vector<std::complex<double> >                         vector_c;
    //typedef ublas::matrix<double,ublas::column_major>                    matrix_d;
    //typedef ublas::matrix<std::complex<double>,ublas::column_major>      matrix_c;
          
}


#endif	/* DEF_H */

