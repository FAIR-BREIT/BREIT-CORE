/* 
 * File:   runSolveSystemAtEquilibrium.cxx
 * Author: winckler
 *
 * Created on July 14, 2015, 11:18 AM
 */

#include "equations_manager.h"
#include "breit_equations.h"
#include "solve_breit_equations.h"
#include "breit_user_interface.h"

using namespace breit;

typedef breit_equations<double> equations_d;
typedef solve_breit_equations<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d> breit_manager;
int main(int argc, char** argv) 
{
    try
    {
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
        
        
    }
    catch(std::exception& e)
    {
        LOG(ERROR) << e.what();
        return 1;
    }
    
    LOG(INFO)<<"Execution successful!";
    return 0;
}
