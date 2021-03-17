/* 
 * File:   breit_equations.h
 * Author: winckler
 *
 * Created on July 13, 2015, 4:51 PM
 */

#ifndef BREIT_EQUATIONS_H
#define	BREIT_EQUATIONS_H

// std
#include <limits>
#include <memory>

// boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

// breit
#include "generate_equations.h"
#include "options_manager.h"
#include "breit_user_interface.h"
#include "def.h"

namespace ublas = boost::numeric::ublas;

namespace breit
{
    
    // CRTP derived class
    template <typename T, typename U=breit_user_interface >
    class breit_equations : public generate_equations<breit_equations<T> >, public U
    {
        typedef T                                                        data_type;  // numerical data type (int, float, double, ...)
        typedef U                                                          ui_type;  // user interface policy (default=options_manager)
        typedef breit_equations<data_type,ui_type>                        self_type;  // this type 
        typedef generate_equations<breit_equations<T> >                  gener_type;  // generate equations policy
        typedef ublas::matrix<data_type,ublas::column_major>              matrix_d;  // boost matrix type
        typedef ublas::vector<double>                                     vector_d;
        using ui_type::parse_cfgFile;
        
        typedef po::options_description                        options_description;
        typedef po::variables_map                                    variables_map;
        typedef fs::path                                                      path;
        
        size_t fEqDim;
        size_t fCoefDim;
        ublas::range fCoef_range_i;
        ublas::range fCoef_range_j;
        ublas::range fSystem_range;
        std::map<std::pair<size_t,size_t>, data_type> fCoef_list;
        matrix_d fMat;
        vector_d f2nd_member;
        vector_d fF0;
        
        
        size_t fCoef_index_min;
        size_t fCoef_index_max;
        std::shared_ptr<breit_summary> fSummary;
        
    protected:
        using ui_type::fvarmap;// boost variable map or equivalent : need signature  "data_type val=fvarmap["key"].as<data_type>();"
        using ui_type::fConfig_file_path;
        using ui_type::fVarmap_input_file;
        using ui_type::N_Avogadro;
        std::map<size_t, data_type> fIni_cond_map;
    public:
        
        breit_equations();
        virtual ~breit_equations(){}
        
        matrix_d& output();
        vector_d& snd_member();
        virtual int parse(const int argc, char** argv, bool AllowUnregistered = false);
        
        // read input and get coefficients of the system
        int read_impl();
        // transform system of equations into matrix equations
        int generate_impl();
        // case dF/dx = AF + cte with dim(A) = N-1
        int dynamic_eq_system();
        //case dF/dx = MF = 0 with dim(M) = N
        int static_eq_system();
        // temp, compute a simple formula taken into account a capture and loss of a single electron (c.f. Betz)
        std::vector<double> get_1electron_approximation_solution();
        
        int init_summary(std::shared_ptr<breit_summary> const& summary) 
        {
            fSummary = summary;
            ui_type::init_summary(summary);
            return 0;
        }
        
        vector_d& initial_condition()
        {
            return fF0;
        }
        
    protected:
        /// ////////////////////////////////////////////////////////////////////////////////
        // Function below are helper functions to compute the matrix element of the final system
        // same as Kronecker-Delta symbol, to help finding the matrix elements of the system of equations
        data_type Fpq(size_t k, size_t q)
        {
            if(k==q)
                return data_type(1);
            else
                return data_type(0);
        }
        // provide element below diagonal of the matrix eq system
        data_type ionization_sum(size_t i, size_t q)
        {
            typedef std::pair<size_t,size_t> coef;
            data_type val=data_type();
            LOG(MAXDEBUG)<<"fCoef_index_min="<<fCoef_index_min;
            LOG(MAXDEBUG)<<"fEqDim="<<fEqDim;
            size_t maxrange = 0;
            if (i > 0) {
                maxrange = i-1;
            }

            for(size_t j(fCoef_index_min);j<=maxrange;j++)
            {
                LOG(MAXDEBUG)<<"i="<<i <<" j="<<j<<"  q="<<q;
                val+=fCoef_list.at(coef(j,i))*Fpq(j,q);
                LOG(MAXDEBUG)<<"val="<<val;
            }
            return val;
        }
        // provide element above diagonal of the matrix eq system
        data_type recombination_sum(size_t i, size_t q)
        {
            typedef std::pair<size_t,size_t> coef;
            data_type val=data_type();
            for(size_t s(i+1);s<=fCoef_index_min+fEqDim-1;s++)
                val+=fCoef_list.at(coef(s,i))*Fpq(s,q);
            return val;
        }
        // provide diagonal element of the matrix eq system
        data_type diagonal_sum(size_t i, size_t q)
        {
            typedef std::pair<size_t,size_t> coef;
            data_type val=data_type();
            for(size_t m(i+1);m<=fCoef_index_min+fEqDim-1;m++) {
                val+=fCoef_list.at(coef(i,m))*Fpq(i,q);
            }

            size_t maxrange = 0;
            if (i > 0) {
                maxrange = i-1;
            }

            for(size_t k(fCoef_index_min);k<=maxrange;k++) {
                val+=fCoef_list.at(coef(i,k))*Fpq(i,q);
            }
            return val;
        }
        data_type compute_matrix_element(size_t p, size_t q)
        {
            data_type val=data_type();
            LOG(MAXDEBUG)<<"ionization_sum";
            val=ionization_sum(p,q);
            LOG(MAXDEBUG)<<"recombination_sum";
            val+=recombination_sum(p,q);
            LOG(MAXDEBUG)<<"diagonal_sum";
            val-=diagonal_sum(p,q);
            return val;
        }
    };
}


#include "breit_equations.tpl"

#endif	/* BREIT_EQUATIONS_H */

