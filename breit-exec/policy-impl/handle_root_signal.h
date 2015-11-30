/* 
 * File:   handle_root_signal.h
 * Author: winckler
 *
 * Created on September 13, 2015, 4:34 PM
 */

#ifndef HANDLE_ROOT_SIGNAL_H
#define	HANDLE_ROOT_SIGNAL_H

#include "TROOT.h"
#include "TPad.h"
#include "TH1.h"
#include "TCanvas.h"
#include "RQ_OBJECT.h"
#include "TApplication.h"


class handle_root_signal 
{
    RQ_OBJECT("handle_root_signal")
public:
    handle_root_signal();
    virtual ~handle_root_signal();
    
    int set_canvas(TCanvas* canvas);
    int set_canvas2(TCanvas* canvas);
    int exit_root();
private:

    TCanvas* fCanvas;
    TCanvas* fCanvas2;
    ClassDef(handle_root_signal, 0);
};

#endif	/* HANDLE_ROOT_SIGNAL_H */

