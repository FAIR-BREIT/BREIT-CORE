/* 
 * File:   runSolveDynEqRoot_RK.cxx
 * Author: winckler
 *
 * Created on August 11, 2015, 2:10 PM
 */

#include "equations_manager.h"
#include "breit_equations.h"
#include "solve_breit_equations_RK.h"
#include "breit_user_interface.h"
#include "breit_gui_root.h"

#include "TApplication.h"

using namespace breit;

typedef breit_equations<double> equations_d;
typedef solve_breit_equations_RK<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d,breit_gui_root> breit_manager;


int main(int argc, char** argv) 
{
    try
    {
        TApplication app("App", nullptr, 0);
        breit_manager man;
        man.use_cfgFile();
        
        LOG(INFO)<<"parsing ...";
        
        if(man.parse(argc, argv,true))
            return 1;
        
        LOG(INFO)<<"initializing ...";
        if(man.init())
            return 1;
        
        LOG(INFO)<<"running ...";
        if(man.run())
            return 1;
        
        LOG(INFO)<<"saving ...";
        if(man.save()) 
            return 1;
        
        
        
        LOG(INFO)<<"plotting ...";
        if(man.plot()) 
            return 1;
        
        app.Run();
        
    }
    catch(std::exception& e)
    {
        LOG(ERROR) << e.what();
        return 1;
    }
    
    LOG(INFO)<<"Execution successful!";
    return 0;
}
