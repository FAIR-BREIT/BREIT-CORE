/* 
 * File:   equations_manager.h
 * Author: winckler
 *
 * Created on July 12, 2015, 9:34 PM
 */

#ifndef EQUATIONS_MANAGER_H
#define	EQUATIONS_MANAGER_H
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "def.h"
#include "logger.h"
namespace breit
{

    
    struct no_gui
    {
        no_gui(){}
        virtual ~no_gui(){}
        
        template <typename... Args> int init(Args&... args){return 0;}
        template <typename... Args> int plot_fig(Args&... args){return 0;}
        template <typename... Args> int init_summary(Args&... args) {return 0;}
        template <typename... Args> int print_table(Args&... args){return 0;}
        template <typename... Args> int save_fig(Args&... args){return 0;}
        template <typename... Args> int compute_equilibrium_distance(){return 0;}
    };
    
    
    template<typename T, typename U, typename V, typename W=no_gui >
    class equations_manager : public U, public V, public W
    {
        typedef T                                    data_type;  // numerical data type (int, float, double, ...)
        typedef U                                      eq_type;  // equation type = generate equation policy
        typedef V                                solve_eq_type;  // solve eq. type = solve equtation policy
        typedef W                                     gui_type;
        typedef equations_manager<T,U,V,W>           self_type;  // this type
        
        std::shared_ptr<breit_summary> fSummary;
        using eq_type::fvarmap;
    public:
        
        equations_manager () :  eq_type(), solve_eq_type(), gui_type()
                                
        {
            fSummary = std::make_shared<breit_summary>();
            eq_type::init_summary(fSummary);
            solve_eq_type::init_summary(fSummary);
            gui_type::init_summary(fSummary);
        }
        
        virtual ~equations_manager (){}

        auto get_options() -> decltype(fvarmap)
        {
            return fvarmap;
        }
        


        int init()
        {
            if(eq_type::init())
                return 1;
            if( solve_eq_type::init(eq_type::fvarmap) )
                return 1;
            
            return 0;
        }

        int run()
        {
            
            solve_eq_type::set_approximated_solution(eq_type::get_1electron_approximation_solution());
            
            if(solve_eq_type::solve(eq_type::output(), 
                                    eq_type::snd_member(), 
                                    eq_type::initial_condition()
                                    ))
                return 1;
            
            //needed for printing table
            gui_type::init(eq_type::fVarmap_input_file,eq_type::fvarmap);
            gui_type::init(solve_eq_type::fGeneral_solution);
            
            
            return 0;
        }
        
        int save()
        {
            INIT_NEW_FILE(fSummary->outfilename,EQUAL,RESULTS);
            double Xmin=eq_type::fVarmap_input_file.at("thickness.minimum").template as<double>();
            double Xmax=eq_type::fVarmap_input_file.at("thickness.maximum").template as<double>();
            size_t Npoint=eq_type::fVarmap_input_file.at("thickness.point.number").template as<std::size_t>();
            
            double Ymin=eq_type::fVarmap_input_file.at("fraction.minimum").template as<double>();
            double Ymax=eq_type::fVarmap_input_file.at("fraction.maximum").template as<double>();
            
            std::string proj_symbol=eq_type::fVarmap_input_file.at("projectile.symbol").template as<std::string>();
            std::string proj_energy=eq_type::fVarmap_input_file.at("projectile.energy").template as<std::string>();
            std::string target_symbol=eq_type::fVarmap_input_file.at("target.symbol").template as<std::string>();
            int tmass=(int)eq_type::fVarmap_input_file.at("target.mass.number").template as<double>();
            std::string target_mass=std::to_string(tmass);
            std::string target_pressure=eq_type::fVarmap_input_file.at("target.pressure").template as<std::string>();
            std::string X_unit=eq_type::fVarmap_input_file.at("thickness.unit").template as<std::string>();
            
            std::stringstream ss;
            
            ss<<proj_symbol
                    <<" projectile at "
                    <<proj_energy
                    //<<" on ^{"
                    <<" on "
                    <<target_mass//<<"}"
                    <<target_symbol
                    <<" target"
                    ;
            ss<<" with "<< target_pressure <<" pressure.";
            
            
            LOG(RESULTS)<<" ";
            LOG(RESULTS)<<"##########################################################################";
            LOG(RESULTS)<<"#                          BREIT  -  RESULTS                              #";
            LOG(RESULTS)<<"##########################################################################";
            LOG(RESULTS)<<" ";
            
            LOG(RESULTS)<<"Title : "<<ss.str();
            LOG(RESULTS)<<" ";
            LOG(RESULTS)<<"Computed from input file : "<<fSummary->filename;
            LOG(RESULTS)<<" ";
            LOG(RESULTS)<<"Found a "<< fSummary->system_dim <<" level system\n";
            LOG(RESULTS)<<" ";
            
            
            bool print_equilibrium=eq_type::fvarmap["save-equilibrium"].template as<bool>();
            
            if(print_equilibrium)
            {
                
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<"#                EQUILIBRIUM CHARGE STATE DISTRIBUTION                   #";
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<" ";
                double sum=0;
                double mean=0;
                for(const auto& p : fSummary->equilibrium_solutions)
                {
                    double Fi=p.second;
                    double qi=fSummary->F_index_map.at(p.first);
                    mean+=qi*Fi;
                    sum+=Fi;
                    LOG(RESULTS)  << "F"
                                << std::to_string(fSummary->F_index_map.at(p.first))
                                << " = "
                                << Fi;
                }
                LOG(RESULTS)<<"sum = "<<sum;
                LOG(RESULTS)<<"<q> = "<<mean;
            }
            ////////////////////////////////////////////////////////////////////////////////////////
            
            bool print_approximation=eq_type::fvarmap["save-approximation"].template as<bool>();
            
            if(print_approximation)
            {
                double sum=0;
                double mean=0;
                LOG(RESULTS)<<" ";
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<"#  EQUILIBRIUM CHARGE STATE DISTRIBUTION  (1-electron approximation)     #";
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<" ";
                for(const auto& p : fSummary->approximated_solutions)
                {
                    double Fi=p.second;
                    double qi=fSummary->F_index_map.at(p.first);
                    mean+=qi*Fi;
                    sum+=Fi;
                    LOG(RESULTS)    << "F"
                                    << std::to_string(fSummary->F_index_map.at(p.first))
                                    << " = "
                                    << Fi;
                }
                LOG(RESULTS)<<"sum = "<<sum;
                LOG(RESULTS)<<"<q> = "<<mean;
            }
                
            
            ////////////////////////////////////////////////////////////////////////////////////////
            
            bool save_table=eq_type::fvarmap["save-table"].template as<bool>();
            bool save_analytic=eq_type::fvarmap["save-analytic"].template as<bool>();
            bool save_fig_ne=eq_type::fvarmap["save-fig-ne"].template as<bool>();
            
            
            if(save_table || save_analytic)
            {
                LOG(RESULTS)<<" ";
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<"#             NON-EQUILIBRIUM CHARGE STATE DISTRIBUTION                  #";
                LOG(RESULTS)<<"##########################################################################";
                LOG(RESULTS)<<" ";
            }
                
            
            
            if(save_analytic)
            {
                LOG(RESULTS)<<"##################################";
                LOG(RESULTS)<<"#computed analytical formulae :";
                for(const auto& p : fSummary->analytical_solutions)
                {
                    LOG(RESULTS)   <<" ";
                    LOG(RESULTS)   << "F"
                                << std::to_string(fSummary->F_index_map.at(p.first))
                                << "(x) = "
                                << p.second;
                    LOG(RESULTS)   <<" ";
                }
            }
            


            bool print_distance=eq_type::fvarmap["save-distance-to-eq"].template as<bool>();
            if(print_distance)
            {
                gui_type::compute_equilibrium_distance();
                LOG(RESULTS)<<"##################################";
                LOG(RESULTS)<<"#estimated distance to equilibrium :";
                for(const auto& p : fSummary->distance_to_equilibrium)
                {
                    if(p.second > 0.0)
                        LOG(RESULTS)<< "X(F"
                                    << std::to_string(fSummary->F_index_map.at(p.first))
                                    << ") = "
                                    << p.second 
                                    << " "
                                    << X_unit;
                    else
                        LOG(RESULTS)<< "X(F"
                                    << std::to_string(fSummary->F_index_map.at(p.first))
                                    << ") = "
                                    << Xmax
                                    << " "
                                    << X_unit
                                    << " (warning : estimation found at range limit."
                                    << " Increase thickness range or relative deviation epsilon if less than 0.001)";
                }
            }
            
            if(save_table)
            {
                LOG(RESULTS)<<" ";
                LOG(RESULTS)<<"##################################";
                LOG(RESULTS)<<"#TABLE :";
                LOG(RESULTS)<<"X unit : "<<X_unit;
                LOG(RESULTS)<<"X range : "<<Xmin<<" - "<<Xmax;
                LOG(RESULTS)<<"Point number : "<<Npoint;
                gui_type::print_table();
            }
            LOG(INFO)<<"- saving output to : "<<fSummary->outfilename;
            
            
            if(save_fig_ne)
            {
                /*
                fs::path input=eq_type::fvarmap["input-file"].template as<fs::path>();
                std::string filename=input.stem().string();
                std::string output=eq_type::fvarmap["output-directory"].template as<fs::path>().string();
                output+="/Bear-results-figure-ne-";
                output+=filename;
                output+=".pdf";
                LOG(INFO)<<"save figure to : "<<output;
                 //*/
                //gui_type::plot();
                //gui_type::save_fig(output);
            }
            
            
            return 0;
        }
        
        int plot()
        {
            gui_type::plot_fig();
            bool save_figure = false;//eq_type::fvarmap["save-fig-ne"].template as<bool>();
            
            
            return 0;
        }
        
        
        
    };
}
#endif	/* EQUATIONS_MANAGER_H */

