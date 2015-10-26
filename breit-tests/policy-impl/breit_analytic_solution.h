/* 
 * File:   breit_analytic_solution.h
 * Author: winckler
 *
 * Created on August 7, 2015, 10:41 AM
 */

#ifndef BREIT_ANALYTIC_SOLUTION_H
#define	BREIT_ANALYTIC_SOLUTION_H


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "def.h"

namespace breit
{
    
    
    template<typename T>
    class breit_analytic_solution
    {
        typedef T                                                              data_type;
        typedef ublas::vector<data_type>                                        vector_d;
        typedef ublas::vector<std::complex<data_type> >                         vector_c;
        typedef ublas::matrix<data_type,ublas::column_major>                    matrix_d;
        typedef ublas::matrix<std::complex<data_type>,ublas::column_major>      matrix_c;
        
        
        typedef std::map<size_t, std::complex<data_type> >                           eigen_value_map;
        typedef std::vector<std::tuple<size_t, size_t, std::complex<double> > > complex_eigen_values;
        
        double fPRECISON;
        std::shared_ptr<breit_summary> fSummary;
    protected:
        std::map<std::size_t, std::string> fGeneral_solution;
        double fUnit_convertor=1.;
        
    public:
        breit_analytic_solution() :  fPRECISON(1.e-15),
                                    fGeneral_solution(),
                                    fUnit_convertor(1.)
        {}
        virtual ~breit_analytic_solution(){}


        int init(const matrix_c& eigen_mat)
        {
            fPRECISON=1.e-15; // temporary, need a real error treatment
            for(size_t row; row<eigen_mat.size1(); row++)
            {
                fGeneral_solution[row]="";
            }
            
            return 0;
        }
        
        int init_summary(std::shared_ptr<breit_summary> const& summary) 
        {
            fSummary = summary;
            return 0;
        }
        int form_general_solution(const vector_d& particular_solution )
        {
            
            for(auto& p : fGeneral_solution)
            {
                p.second+="+";
                p.second+=to_string_scientific(particular_solution(p.first));
            }
            
            
            
            std::size_t last_key=0;
            std::size_t first_key=0;
            if( !fGeneral_solution.empty() )
            {
                first_key=fGeneral_solution.begin()->first;
                last_key=fGeneral_solution.rbegin()->first;
            }
            
            last_key+=1;// because store keys as 0, 1 ... N-2
            // 8 lvl system -> 0, 1, ... 6, and we print 1, ... 7
            // need to add the last eq F8= 1 - sum of the other states
            
            if(last_key!=0 && !fGeneral_solution.count(last_key))
            {
                std::string F_last("1. - (");
                for(const auto& p : fGeneral_solution)
                {
                    if(p.first != first_key)
                    F_last+=" + ";
                    F_last+=p.second;
                    
                }
                    F_last+=")";
                fGeneral_solution[last_key]=F_last;
            }
            
            // Copy final solution
            fSummary->analytical_solutions=fGeneral_solution;
            
            
            
            
            return 0;
        }
        
        
        int form_homogeneous_solution(  const matrix_c& eigen_mat, 
                                        const vector_d& unknown_coef, 
                                        eigen_value_map& eig_val_map, 
                                        complex_eigen_values& eig_val_c)
        {
            
            
            
            // CASE = complex eigenvectors
            // C_k exp(lambda_k x) * ev_k
            // with ev_k(i) = ai cos(omega_k x) - bi sin(omega_k x) 
            
            // and ev_k'(i) bar = ai sin(omega_k' x) + bi cos(omega_k' x)
            for(const auto& p : eig_val_c)
            {
                size_t index=0;
                size_t index_bar=0;
                std::complex<double> eigenvalue;
                std::tie(index,index_bar,eigenvalue) = p;
                
                double lambda=eigenvalue.real()*fUnit_convertor;
                double omega=eigenvalue.imag()*fUnit_convertor;
                
                std::string expLambdaX;
                std::string coswx;
                std::string sinwx;
                
                // to simplify string solution handle case where lambda and omega = +/- 1
                if(std::fabs(lambda-1.)<fPRECISON)
                {
                    expLambdaX="exp(x)";
                }
                else
                {
                    if(std::fabs(lambda+1.)<fPRECISON)
                        expLambdaX="exp(-x)";
                    else
                        expLambdaX="exp(" + to_string_scientific(lambda)+"*x)";
                }
                
                if(std::fabs(omega-1.)<fPRECISON)
                {
                    sinwx="sin(x)";
                    coswx="cos(x)";
                }
                else
                {
                    if(std::fabs(omega+1.)<fPRECISON)
                    {
                        sinwx="sin(-x)";
                        coswx="cos(x)";
                    }
                    else
                    {
                        sinwx="sin(" + to_string_scientific(omega) + "*x)";
                        coswx="cos(" + to_string_scientific(omega) + "*x)";
                    }
                }
                
                for(size_t row(0); row<eigen_mat.size1(); row++)
                {
                    // ///////////////////////////////
                    // contribution from eigenvector :
                    // coef of first eigenvector
                    double ai_val=eigen_mat(row,index).real();
                    double bi_val=eigen_mat(row,index).imag();
                    double C1_val=unknown_coef(index);
                    double C2_val=unknown_coef(index_bar);
                    
                    double C1xai_val=C1_val*ai_val;
                    double C1xbi_val=C1_val*bi_val;
                    
                    double C2xai_val=C2_val*ai_val;
                    double C2xbi_val=C2_val*bi_val;
                    
                    std::string ai;
                    std::string bi;
                    std::string C1;
                    std::string C2;
                    
                    std::string C1xai=to_string_scientific(C1xai_val);
                    std::string C1xbi=to_string_scientific(C1xbi_val);
                    std::string C2xai=to_string_scientific(C2xai_val);
                    std::string C2xbi=to_string_scientific(C2xbi_val);
                    
                    
                    
                   
                    ai=to_string_scientific(ai_val);
                    bi=to_string_scientific(bi_val);
                    
                    C1=to_string_scientific(C1_val);
                    C2=to_string_scientific(C2_val);
                    
                    // 
                    if(std::fabs(C1_val)>fPRECISON)                                                     //                  [if C1!=0]
                    {
                        if(!fGeneral_solution[row].empty())                             
                            fGeneral_solution[row]+="+";                                    // +                [if string not empty]

                        //fGeneral_solution[row]+="(";                                        // (                [if C1!=0]
                        //fGeneral_solution[row]+=C1;                                         // C1               [if C1!=0]
                        //fGeneral_solution[row]+=")";                                        // )                [if C1!=0]
                        //fGeneral_solution[row]+="*";                                        // *                [if C1!=0]
                        
                        if(std::fabs(lambda)>fPRECISON)
                            fGeneral_solution[row]+= expLambdaX + "*";                      // exp(lambda x) *  [if C1!=0] [if lambda!=0]
                        // component of first eigenvector
                        fGeneral_solution[row]+="(";                                        // (
                        if(std::fabs(ai_val)>fPRECISON)
                            fGeneral_solution[row]+="(" + C1xai + ")";                         // (ai)             [if C1!=0] [if ai!=0]
                        if(std::fabs(omega)>fPRECISON)
                        {
                            if(std::fabs(ai_val)>fPRECISON)
                                fGeneral_solution[row]+="*"+coswx;                          // * cos(omega x)   [if C1!=0] [if omega!=0] [if ai!=0]
                            if(std::fabs(bi_val)>fPRECISON)
                            {
                                fGeneral_solution[row]+="-1.*";                                // -                [if C1!=0] [if omega!=0] [if bi!=0]
                                fGeneral_solution[row]+="(" + C1xbi + ")";                     // (bi)             [if C1!=0] [if omega!=0] [if bi!=0]
                                fGeneral_solution[row]+="*"+sinwx;                          // * sin(omega x)   [if C1!=0] [if omega!=0] [if bi!=0]
                            }
                        }
                        fGeneral_solution[row]+=")";                                        // )
                    }
                    // ///////////////////////////////
                    // contribution from eigenvector' complex conjugate :
                    // coef of first eigenvector
                    
                    if(std::fabs(C2_val)>fPRECISON)                                                     //                  [if C2!=0]
                    {
                        if(!fGeneral_solution[row].empty()) 
                            fGeneral_solution[row]+="+";                                    // +                [if string not empty]

                        //fGeneral_solution[row]+="(";                                        // (                [if C2!=0]
                        //fGeneral_solution[row]+=C2;                                         // C2               [if C2!=0]
                        //fGeneral_solution[row]+=")";                                        // )                [if C2!=0]
                        //fGeneral_solution[row]+="*";                                        // *                [if C2!=0]
                        
                        if(std::fabs(lambda)>fPRECISON)
                            fGeneral_solution[row]+= expLambdaX + "*";                      // exp(lambda x) *  [if C2!=0] [if lambda!=0]
                        // component of first eigenvector
                        fGeneral_solution[row]+="(";                                        // (
                        if(std::fabs(bi_val)>fPRECISON)
                            fGeneral_solution[row]+="(" + C2xbi + ")";                         // (bi)             [if C2!=0] [if bi!=0]
                        if(std::fabs(omega)>fPRECISON)
                        {
                            if(std::fabs(bi_val)>fPRECISON)
                                fGeneral_solution[row]+="*"+coswx;                          // * cos(omega x)   [if C2!=0] [if omega!=0] [if bi!=0]
                            if(std::fabs(ai_val)>fPRECISON)
                            {
                                fGeneral_solution[row]+="+1.*";                                // -                [if C2!=0] [if omega!=0] [if ai!=0]
                                fGeneral_solution[row]+="(" + C2xai + ")";                     // (ai)             [if C2!=0] [if omega!=0] [if ai!=0]
                                fGeneral_solution[row]+="*"+sinwx;                          // * sin(omega x)   [if C2!=0] [if omega!=0] [if ai!=0]
                            }
                        }
                        fGeneral_solution[row]+=")";                                        // )
                    }
                    
                    
                }
            }
            
            // CASE = real eigenvectors
            // C_k exp(lambda_k x) * ev_k
            // with ev_k(i) = ai
            for(const auto& p : eig_val_map)
            {
                double lambda=p.second.real()*fUnit_convertor;
                size_t index=p.first;
                std::string expLambdaX;
                
                // to simplify string solution handle case where lambda and omega = +/- 1
                if(std::fabs(lambda-1.0)<fPRECISON)
                {
                    expLambdaX="exp(x)";
                }
                else
                {
                    if(std::fabs(lambda+1.0)<fPRECISON)
                        expLambdaX="exp(-x)";
                    else
                        expLambdaX="exp(" + to_string_scientific(lambda)+"*x)";
                }
                
                for(size_t row(0); row<eigen_mat.size1(); row++)
                {
                    double ai_val=eigen_mat(row,index).real();
                    double C1_val=unknown_coef(index);
                    double C1ai_val=C1_val*ai_val;
                    
                    
                    std::string ai;// ai == P(row=i,ev_index)
                    std::string C1;
                    std::string C1ai;
                    ai=to_string_scientific(ai_val);
                    C1=to_string_scientific(C1_val);
                    C1ai=to_string_scientific(C1ai_val);
                    
                    if(std::fabs(C1_val)>fPRECISON && std::fabs(ai_val)>fPRECISON)
                    {
                        if(!fGeneral_solution[row].empty())
                            fGeneral_solution[row]+="+";

                        
                        fGeneral_solution[row]+="(" + C1ai +")";                                             // C1
                        //fGeneral_solution[row]+="* ";                                           // * 
                        if(std::fabs(lambda)>fPRECISON)
                            fGeneral_solution[row]+= "*" + expLambdaX;                          // exp(lambda x) *   [if lambda!=0]
                        //fGeneral_solution[row]+=ai;                                             // ai
                    }
                }
            }
            
            
            /*for(const auto& p : fGeneral_solution)
            {
                LOG(INFO)<<"F"<< p.first+1 <<"(x) = "<< p.second;
            }*/
            
            
            return 0;
        }
        
    
        
        
        
        int form_homogeneous_solution_raw(  const matrix_c& eigen_mat, 
                                        const vector_d& unknown_coef, 
                                        eigen_value_map& eig_val_map, 
                                        complex_eigen_values& eig_val_c)
        {
            // CASE = complex eigenvectors
            // C_k exp(lambda_k x) * ev_k
            // with ev_k(i) = ai cos(omega_k x) - bi sin(omega_k x) 
            
            // and ev_k'(i) bar = ai sin(omega_k' x) + bi cos(omega_k' x)
            for(const auto& p : eig_val_c)
            {
                size_t index=0;
                size_t index_bar=0;
                std::complex<double> eigenvalue;
                std::tie(index,index_bar,eigenvalue) = p;
                
                double lambda=eigenvalue.real()*fUnit_convertor;
                double omega=eigenvalue.imag()*fUnit_convertor;
                
                std::string expLambdaX;
                std::string coswx;
                std::string sinwx;
                
                // to simplify string solution handle case where lambda and omega = +/- 1
                
                expLambdaX="exp(" + to_string_scientific(lambda)+"*x)";
                sinwx="sin(" + to_string_scientific(omega) + "*x)";
                coswx="cos(" + to_string_scientific(omega) + "*x)";
                        
                for(size_t row(0); row<eigen_mat.size1(); row++)
                {
                    // ///////////////////////////////
                    // contribution from eigenvector :
                    // coef of first eigenvector
                    double ai_val=eigen_mat(row,index).real();
                    double bi_val=eigen_mat(row,index).imag();
                    double C1_val=unknown_coef(index);
                    double C2_val=unknown_coef(index_bar);
                    
                    double C1xai_val=C1_val*ai_val;
                    double C1xbi_val=C1_val*bi_val;
                    
                    double C2xai_val=C2_val*ai_val;
                    double C2xbi_val=C2_val*bi_val;
                    
                    std::string ai;
                    std::string bi;
                    std::string C1;
                    std::string C2;
                    
                    std::string C1xai=to_string_scientific(C1xai_val);
                    std::string C1xbi=to_string_scientific(C1xbi_val);
                    std::string C2xai=to_string_scientific(C2xai_val);
                    std::string C2xbi=to_string_scientific(C2xbi_val);
                    
                    
                    
                   
                    ai=to_string_scientific(ai_val);
                    bi=to_string_scientific(bi_val);
                    
                    C1=to_string_scientific(C1_val);
                    C2=to_string_scientific(C2_val);
                    
                    // 
                    
                    if(!fGeneral_solution[row].empty())                             
                        fGeneral_solution[row]+="+";                                    // +                [if string not empty]

                    //fGeneral_solution[row]+="(";                                        // (                [if C1!=0]
                    //fGeneral_solution[row]+=C1;                                         // C1               [if C1!=0]
                    //fGeneral_solution[row]+=")";                                        // )                [if C1!=0]
                    //fGeneral_solution[row]+="*";                                        // *                [if C1!=0]

                    fGeneral_solution[row]+= expLambdaX + "*";                      // exp(lambda x) *  [if C1!=0] [if lambda!=0]
                    // component of first eigenvector
                    fGeneral_solution[row]+="(";                                        // (
                    fGeneral_solution[row]+="(" + C1xai + ")";                         // (ai)             [if C1!=0] [if ai!=0]
                    
                    fGeneral_solution[row]+="*"+coswx;                          // * cos(omega x)   [if C1!=0] [if omega!=0] [if ai!=0]
                    fGeneral_solution[row]+="-";                                // -                [if C1!=0] [if omega!=0] [if bi!=0]
                    fGeneral_solution[row]+="(" + C1xbi + ")";                     // (bi)             [if C1!=0] [if omega!=0] [if bi!=0]
                    fGeneral_solution[row]+="*"+sinwx;                          // * sin(omega x)   [if C1!=0] [if omega!=0] [if bi!=0]
                    fGeneral_solution[row]+=")";                                        // )
                    
                    // ///////////////////////////////
                    // contribution from eigenvector' complex conjugate :
                    // coef of first eigenvector
                    
                    if(!fGeneral_solution[row].empty()) 
                        fGeneral_solution[row]+="+";                                    // +                [if string not empty]

                    //fGeneral_solution[row]+="(";                                        // (                [if C2!=0]
                    //fGeneral_solution[row]+=C2;                                         // C2               [if C2!=0]
                    //fGeneral_solution[row]+=")";                                        // )                [if C2!=0]
                    //fGeneral_solution[row]+="*";                                        // *                [if C2!=0]

                    fGeneral_solution[row]+= expLambdaX + "*";                      // exp(lambda x) *  [if C2!=0] [if lambda!=0]
                    // component of first eigenvector
                    fGeneral_solution[row]+="(";                                        // (
                    fGeneral_solution[row]+="(" + C2xbi + ")";                         // (bi)             [if C2!=0] [if bi!=0]
                    
                    fGeneral_solution[row]+="*"+coswx;                          // * cos(omega x)   [if C2!=0] [if omega!=0] [if bi!=0]

                    fGeneral_solution[row]+="+";                                // -                [if C2!=0] [if omega!=0] [if ai!=0]
                    fGeneral_solution[row]+="(" + C2xai + ")";                     // (ai)             [if C2!=0] [if omega!=0] [if ai!=0]
                    fGeneral_solution[row]+="*"+sinwx;                          // * sin(omega x)   [if C2!=0] [if omega!=0] [if ai!=0]
                    fGeneral_solution[row]+=")";                                        // )
                    
                    
                }
            }
            
            // CASE = real eigenvectors
            // C_k exp(lambda_k x) * ev_k
            // with ev_k(i) = ai
            for(const auto& p : eig_val_map)
            {
                double lambda=p.second.real()*fUnit_convertor;
                size_t index=p.first;
                std::string expLambdaX;
                
                // to simplify string solution handle case where lambda and omega = +/- 1
                
                expLambdaX="exp(" + to_string_scientific(lambda)+"*x)";
                
                for(size_t row(0); row<eigen_mat.size1(); row++)
                {
                    double ai_val=eigen_mat(row,index).real();
                    double C1_val=unknown_coef(index);
                    double C1ai_val=C1_val*ai_val;
                    
                    
                    std::string ai;// ai == P(row=i,ev_index)
                    std::string C1;
                    std::string C1ai;
                    ai=to_string_scientific(ai_val);
                    C1=to_string_scientific(C1_val);
                    C1ai=to_string_scientific(C1ai_val);
                    
                    if(!fGeneral_solution[row].empty())
                        fGeneral_solution[row]+="+";


                    fGeneral_solution[row]+="(" + C1ai +")";                                             // C1
                    //fGeneral_solution[row]+="* ";                                           // * 
                    fGeneral_solution[row]+= "*" + expLambdaX;                          // exp(lambda x) *   [if lambda!=0]
                    //fGeneral_solution[row]+=ai;                                             // ai
                }
            }
            
            
            for(const auto& p : fGeneral_solution)
            {
                LOG(INFO)<<"F"<< p.first+1 <<"(x) = "<< p.second;
            }
            
            
            return 0;
        }
        
        
        
        
        
        
    };
}
#endif	/* BREIT_ANALYTIC_SOLUTION_H */

