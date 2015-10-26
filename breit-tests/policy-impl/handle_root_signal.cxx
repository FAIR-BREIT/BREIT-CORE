/* 
 * File:   handle_root_signal.cxx
 * Author: winckler
 * 
 * Created on September 13, 2015, 4:34 PM
 */

#include "handle_root_signal.h"

handle_root_signal::handle_root_signal() : fCanvas(NULL)
{
}



handle_root_signal::~handle_root_signal() 
{
}

int handle_root_signal::set_canvas(TCanvas* canvas)
{
    fCanvas=canvas;
    fCanvas->Connect("Closed()","handle_root_signal",this,"exit_root()");
}

int handle_root_signal::set_canvas2(TCanvas* canvas)
{
    fCanvas2=canvas;
    fCanvas2->Connect("Closed()","handle_root_signal",this,"exit_root()");
}

int handle_root_signal::exit_root()
{
    gApplication->Terminate(0);
    return 0;
}



ClassImp(handle_root_signal)