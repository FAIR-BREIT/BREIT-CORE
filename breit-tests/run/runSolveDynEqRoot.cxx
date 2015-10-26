/* 
 * File:   runSolveDynEqRoot.cxx
 * Author: winckler
 *
 * Created on August 10, 2015, 8:07 PM
 */

#include "equations_manager.h"
#include "breit_equations.h"
#include "solve_breit_equations.h"
#include "breit_user_interface.h"
#include "breit_gui_root.h"
#include "logger.h"
#include "TApplication.h"

#include "def.h"

using namespace breit;

typedef breit_equations<double> equations_d;
typedef solve_breit_equations<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d,breit_gui_root> breit_manager;


int main(int argc, char** argv) 
{
    try
    {
        init_log_console(breit::severity_level::INFO,log_op::operation::GREATER_EQ_THAN);
        /// /////////////////////////////////////////////////////
        // CREATE EQUATION MANAGER
        LOG(STATE)<<"start BREIT : Ballance Equations for Atomic Reactions";
        breit_manager man;
        man.use_cfgFile();
        
        /// /////////////////////////////////////////////////////
        // PARSE OPTIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"parsing command line and input files ...";
        LOG(INFO)<<" ";
        if(man.parse(argc, argv,true))
            return 1;
        
        bool plot = man.get_varMap()["plot"].as<bool>();
        bool save = man.get_varMap()["save"].as<bool>();
        bool save_fig_ne = man.get_varMap()["save-fig-ne"].as<bool>();
        bool save_fig_e = man.get_varMap()["save-fig-e"].as<bool>();
        bool save_root_ne = man.get_varMap()["save-root-ne"].as<bool>();
        bool save_root_e = man.get_varMap()["save-root-e"].as<bool>();
        /// /////////////////////////////////////////////////////
        // INIT EQUATIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"initializing ...";
        if(man.init())
            return 1;
        
        /// /////////////////////////////////////////////////////
        // RUN SOLVE EQUATIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"running ...";
        if(man.run())
            return 1;
        
        
        
        /// /////////////////////////////////////////////////////
        // SAVE
        if(save)
        {
            LOG(INFO)<<" ";
            LOG(STATE)<<"saving to file ...";
            if(man.save()) 
                return 1;
        }
            
        /// /////////////////////////////////////////////////////
        // PLOT
        
        if(plot)
        {
            TApplication app("App", nullptr, 0);
            LOG(INFO)<<" ";
            LOG(STATE)<<"plotting ...";
            if(man.plot()) 
                return 1;
            // run event loop
            app.Run();
        }
        else
        {
            if(save_fig_ne || save_fig_e || save_root_ne || save_root_e)
            {
                
                if(man.plot()) 
                    return 1;
            }
        }
        
    }
    catch(std::exception& e)
    {
        LOG(ERROR) << e.what();
        return 1;
    }
    
    LOG(INFO)<<"Execution successful!";
    return 0;
}
