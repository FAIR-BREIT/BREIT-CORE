/* 
 * File:   solve_breit_equations_RK.h
 * Author: winckler
 *
 * Created on August 11, 2015, 12:10 PM
 */

#ifndef SOLVE_BREIT_EQUATIONS_RK_H
#define	SOLVE_BREIT_EQUATIONS_RK_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "options_manager.h"
#include "matrix_inverse.hpp"
#include "storage_adaptors.hpp"
#include "matrix_diagonalization.h"

#include "def.h"
#include "breit_analytic_solution.h"

#include "TH1D.h"


#include <iostream>
#include <fstream>
#include <sstream>

#include <tuple>
#include <memory>
#include <algorithm>
#include <vector>


#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"


namespace breit
{
    
    
    
    
    template<typename T>
    class solve_breit_equations_RK 
    {
        public:
        
        enum class diagonalizable
        {
            unknown,
            in_R,
            in_C,
            no
        };
        
        private:
          typedef T                                                              data_type;
          typedef ublas::vector<data_type>                                        vector_d;
          typedef ublas::vector<std::complex<data_type> >                         vector_c;
          typedef ublas::matrix<data_type,ublas::column_major>                    matrix_d;
          typedef ublas::matrix<std::complex<data_type>,ublas::column_major>      matrix_c;
          typedef po::variables_map                                          variables_map;
          //  equation to solve : dF/dx = AF + g <=> dF/dx = P D P^1 F + g
          matrix_d fA;                       // A
          matrix_d fA_inv;                   // A^-1
          vector_d f2nd_member;              // g
          vector_d fF0;                      // Fi(x=0)
          vector_c fD;                       // D (stored in a vector, store only eigenvalues)
          matrix_c fEigen_mat;               // P
          matrix_c fEigen_mat_inv;           // P^-1
         
          vector_c fConstant_set;            // unknown coefficient from integration, to be determined with initial conditions
         
          vector_d fEquilibrium_solution;    // Fi at equilibrium, should be eq to part. sol.
          //vector_d fGeneral_solution;    // homogeneous solution
          //vector_d fGeneral_solution;        // general solution = homogeneous + particular solution
          enum diagonalizable diagonalisation_case;
          variables_map fvarmap;
          std::vector<double> fApproximated_solution;
          std::shared_ptr<breit_summary> fSummary;
    protected:
          std::map<std::size_t, std::shared_ptr<TH1D> > fGeneral_solution;
          //TCanvas* fCanvas;
          //TLegend* fLegend;
    public:
        
        solve_breit_equations_RK() : 
                                 fA(), 
                                 fA_inv(), 
                                 f2nd_member(),
                                 fF0(),
                                 fD(), 
                                 fEigen_mat(), 
                                 fEigen_mat_inv(),
                                 fEquilibrium_solution(),
                                 //fGeneral_solution(),
                                 diagonalisation_case(diagonalizable::unknown),
                                 fvarmap(),fGeneral_solution(), fApproximated_solution()
        {}
        virtual ~solve_breit_equations_RK(){}
        
        int init(const variables_map& vm)
        {
            fvarmap=vm;
            fs::path input=fvarmap["input-file"].template as<fs::path>();
            std::string filename=input.filename().string();
            std::string output=fvarmap["output-directory"].template as<fs::path>().string();
            output+="/Results_";
            output+=filename;
            LOG(INFO)<<"Print output to : "<<output;
            fSummary->filename=output;
            LOG(RESULTS)<<"Input file :"<<filename<<"\n";
            return 0;
        }
        int init_summary(std::shared_ptr<breit_summary> const& summary) 
        {
            fSummary = summary;
            return 0;
        }
        // main function 
        int solve(const matrix_d& mat, const vector_d& vec, const vector_d& initial_condition)
        {
            try
            {
                reset_to(mat);
                solve_dyneq_at_equilibrium(mat,vec);
                solve_dynamic_system(mat,vec);
            }
            catch(std::exception& e)
            {
                LOG(ERROR)<< "could not solve system. Reason : " << e.what();
                return 1;
            }
            return 0;
        }
        
        
        // temp
        int print_analytical_solution(const std::vector<double>& vec)
        {
            
            LOG(RESULTS)<<"ANALYTICAL SOLUTION\n";//<<std::endl;
            for(int i(0);i<vec.size()-1;i++)
                LOG(RESULTS)<<"F"<<i+1<<"="<<vec.at(i)<<"\n";//<<std::endl;
            
            LOG(RESULTS)<<"sum="<<vec.at(vec.size()-1)<<"\n";//<<std::endl;
            
            return 0;
        }
        
        
        int set_approximated_solution(const std::vector<double>& vec)
        {
            fApproximated_solution=vec;
            if(fApproximated_solution.size()<1)
                return 1;
            
            return 0;
            
        }
        // temp
        int print_approximated_solution()
        {
            if(fApproximated_solution.size()<1)
                return 1;
            
            LOG(RESULTS)<<" ";
            LOG(RESULTS)<<"##########################################################################";
            LOG(RESULTS)<<"#  EQUILIBRIUM CHARGE STATE DISTRIBUTION  (1-electron approximation)     #";
            LOG(RESULTS)<<"##########################################################################";
            LOG(RESULTS)<<" ";
            for(int i(0);i<fApproximated_solution.size()-1;i++)
                LOG(RESULTS)<<"F"<<i+1<<"="<<fApproximated_solution.at(i);
            LOG(RESULTS)<<"sum="<<fApproximated_solution.at(fApproximated_solution.size()-1);
            
            
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation 
        
        
        
        int solve_dynamic_system(const matrix_d& mat, const vector_d& vec)
        {
            vector_d F0(mat.size1());
            vector_d Fnp1(mat.size1());
            vector_d Fn(mat.size1());
            
            double dx=0.001;
            double xmin=0.0;
            double xmax=20.0;
            size_t it_number=size_t( (xmax-xmin)/dx );
            LOG(INFO)<<"Number of iterations = "<<it_number;
            
            //std::make_tuple
            //fLegend = new TLegend(0.1,0.7,0.2,0.9);
            //init
            size_t non_zero_coord=0;
            for(size_t i(0); i<F0.size(); i++)
            {
                std::string name = "F" + std::to_string(i+1);
                //fGeneral_solution[i] = new TH1D(name.c_str(),name.c_str(),it_number,xmin,xmax);
                fGeneral_solution[i] = std::make_shared<TH1D>(name.c_str(),name.c_str(),it_number,xmin,xmax);
                fGeneral_solution.at(i)->SetLineColor(i+1);
                //fLegend->AddEntry(fGeneral_solution.at(i), name.c_str());
                if(non_zero_coord!=i)
                    F0(i)=0.0;
                else
                    F0(i)=1.0;
            }
            
            
            Fn=F0;
            for(size_t i(0); i<it_number; i++)
            {
                vector_d temp(mat.size1());
                
                temp=prod(mat,Fn);
                Fnp1=temp*dx + vec*dx + Fn;
                //Fnp1=temp*dx + Fn;// homogeneous equation
                for(auto& p : fGeneral_solution)
                {
                    fGeneral_solution.at(p.first)->SetBinContent(i,Fn(p.first));
                }
                Fn=Fnp1;
            }
            //plot();
            return 0;
        }
        /*
        int plot()
        {
            LOG(DEBUG)<<"GUI start";
            fCanvas = new TCanvas("c1","Solutions",800,600);
            for(const auto& p : fGeneral_solution)
            {
                
                if(p.first!=0)
                    p.second->Draw("SAME");
                else
                {
                    //p.second->GetYaxis()->SetRangeUser(0., 2.);
                    p.second->Draw();
                }
            }
            fLegend->Draw();
            return 0;
        }
        //*/
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation at equilibrium method 2 (dim=N-1)
        // dynamic equations : case at equilibrium F = -A^1 g 
        // this solution should be a particular solution to the general system
        int solve_dyneq_at_equilibrium(const matrix_d& mat, const vector_d& vec)
        {
            LOG(MAXDEBUG)<<"running solve dynamic eq";
            LOG(RESULTS)<<"found a "<<mat.size1()+1<<" level system\n";
            
            LOG(RESULTS)<<"SOLUTION AT EQUILIBRIUM :\n";
            fA=mat;
            f2nd_member=vec;
            InvertMatrix<matrix_d>(fA,fA_inv);

            //std::cout << fA_inv << std::endl;
            // todo : need to get 
            // -index range
            double sum=0.0;
            double FN=1.0;
            vector_d neg_Fi = prod(fA_inv, f2nd_member);// dim N-1
            
            for(size_t i(0); i< neg_Fi.size(); i++)
            {
                fEquilibrium_solution(i)=-neg_Fi(i);
                FN+=neg_Fi(i);
                sum+=fEquilibrium_solution(i);
                LOG(INFO)<<"F"<<i+1<<"="<<fEquilibrium_solution(i);
                LOG(RESULTS)<<"F"<<i+1<<"="<<fEquilibrium_solution(i)<<"\n";
            }
            // add the last one (1-sum)
            fEquilibrium_solution(neg_Fi.size())=FN;
            sum+=FN;
            
            LOG(INFO)<<"F"<< neg_Fi.size()+1<<"="<<fEquilibrium_solution(neg_Fi.size());
            LOG(INFO)<<"sum = "<< sum;
            LOG(RESULTS)<<"F" << neg_Fi.size() + 1 << "="<<fEquilibrium_solution(neg_Fi.size())<<"\n";
            LOG(RESULTS)<<"sum = "<< sum<<"\n";
            return 0;
        }
        
        int reset_to(const matrix_d& mat)
        {
            // check first if input matrix is a square matrix
            LOG(MAXDEBUG)<<"input matrix dimension = "<<mat.size1()<<"x"<<mat.size2();
            if(mat.size1()!=mat.size2())
            {
                LOG(ERROR) << "input matrix is not a square matrix (dim1 = "
                           << mat.size1() 
                           << ", dim2 = " 
                           << mat.size2() 
                           << ").";
                return 1;
            }
            
            /// input to copy
            fA.clear();
            f2nd_member.clear();
            fF0.clear();
            
            fA.resize(mat.size1(),mat.size2());
            f2nd_member.resize(mat.size1());
            fF0.resize(mat.size1()+1);// to check
            
            /// intermediate matrix and vectors
            fA_inv.clear();
            fD.clear();
            fEigen_mat.clear();
            fEigen_mat_inv.clear();
            fConstant_set.clear();
            
            fA_inv.resize(mat.size1(),mat.size2());
            fD.resize(mat.size1());
            fEigen_mat.resize(mat.size1(),mat.size2());
            fEigen_mat_inv.resize(mat.size1(),mat.size2());
            fConstant_set.resize(mat.size1());
            
            /// vector solutions 
            // case we use the "dynamic equations" the dimensions of the vector solution are mat.size+1
            // otherwise they have same dimension
            
            fEquilibrium_solution.clear();
            //fGeneral_solution.clear();
            //fGeneral_solution.clear();
            
            fEquilibrium_solution.resize(mat.size1()+1);
            //fGeneral_solution.resize(mat.size1()+1);
            //fGeneral_solution.resize(mat.size1()+1);
            
            return 0;
        }
        
    

    };
    
    
}

#endif	/* SOLVE_BREIT_EQUATIONS_RK_H */

