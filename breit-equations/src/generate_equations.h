/* 
 * File:   GenerateEquations.h
 * Author: winckler
 *
 * Created on July 12, 2015, 9:01 PM
 */

#ifndef GENERATEEQUATIONS_H
#define	GENERATEEQUATIONS_H

#include <string>
#include "logger.h"
namespace breit
{

    //  CRTP base class
    template <typename TDerived >
    class generate_equations
    {
    public:
        generate_equations() 
        {}

        virtual ~generate_equations()
        {}

        int read()
        {
            return static_cast<TDerived*>(this)->read_impl();
        }

        int generate()
        {
            return static_cast<TDerived*>(this)->generate_impl();
        }

        int init()
        {
            try
            {
                if(read())
                    return 1;
                if(generate())
                    return 1;
            }
            catch(std::exception& e)
            {
                LOG(ERROR) << e.what();
                return 1;
            }
            return 0;
        }
    };
}
#endif	/* GENERATEEQUATIONS_H */

