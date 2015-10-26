/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

/* 
 * File:   options_manager.h
 * Author: winckler
 *
 * Created on July 13, 2015, 10:30 PM
 */
// refactored from FairRoot::FairMQ project

#ifndef OPTIONS_MANAGER_H
#define	OPTIONS_MANAGER_H

// std
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <tuple>

// breit
#include "logger.h"
#include "def.h"
namespace breit
{
    template<class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
    {
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, "  "));
        return os;
    }

    //namespace po = boost::program_options;
    //namespace fs = boost::filesystem;
    class options_manager
    {
        typedef po::options_description                         options_description;
        typedef po::variables_map                                     variables_map;
        typedef po::variable_value                                   variable_value;
        typedef fs::path                                                       path;
    public:
        enum class source {cmd, cfg, file, env};
        options_manager();
        virtual ~options_manager();

        //  add options_description
        int addTo_cmdLine(const options_description& optdesc, bool visible = true);
        int addTo_cfgFile(const options_description& optdesc, bool visible = true);
        int addTo_env(const options_description& optdesc);

        int add_to();
        void use_cfgFile(const std::string& filename = "");

        // get value corresponding to the key
        template<typename T>
        T get_value(const std::string& key) const
        {
            T val = T();
            try
            {
                if (fvarmap.count(key))
                {
                    val = fvarmap[key].as<T>();
                }
            }
            catch(std::exception& e)
            {
                LOG(ERROR) << "Exception thrown for the key '" << key << "'";
                LOG(ERROR) << e.what();
                this->print_help();
            }

            return val;
        }

        // convert value to string that corresponds to the key
        std::string get_stringVal(const std::string& key);

        const variables_map& get_varMap() const {return fvarmap;}

        // boost prog options parsers
        int parse_cmdLine(const int argc, char** argv, const options_description& desc, variables_map& varmap, bool AllowUnregistered = false);
        int parse_cmdLine(const int argc, char** argv, const options_description& desc, bool AllowUnregistered = false);

        int parse_cfgFile(const std::string& filename, const options_description& desc, variables_map& varmap, bool AllowUnregistered = false);
        int parse_cfgFile(const std::string& filename, const options_description& desc, bool AllowUnregistered = false);
        int parse_cfgFile(std::ifstream& ifs, const options_description& desc, variables_map& varmap, bool AllowUnregistered = false);
        int parse_cfgFile(std::ifstream& ifs, const options_description& desc, bool AllowUnregistered = false);

        int parse_environment(const std::function<std::string(std::string)>&);

        virtual int parse(const int argc, char** argv, bool AllowUnregistered = false) = 0 ;

        virtual int print_options();
        int print_help() const;
        int init_summary(std::shared_ptr<breit_summary> const& summary) 
        {
            fSummary = summary;
            return 0;
        }

        
        
    protected:
        // options container
        variables_map fvarmap;

        // basic description categories
        options_description fGenericDesc;
        options_description fConfigDesc;
        options_description fEnvironmentDesc;
        options_description fHiddenDesc;

        // Description of cmd line and simple configuration file (configuration file like txt, but not like xml json ini)
        options_description fCmdline_options;
        options_description fConfig_file_options;

        // Description which is printed in help command line
        options_description fVisible_options;

        std::string fVerboseLvl;
        bool fUse_cfgFile;
        path fConfig_file_path;
        
        std::map<std::string,bool> fVisible_key_map;
        
        
        virtual int notifySwitch_options();

        // update_varMap() and replace() --> helper functions to modify the value of variable map after calling po::store
        template<typename T>
        void update_varMap(const std::string& key, const T& val)
        {
            replace(fvarmap, key, val);
        }

        template<typename T>
        void replace(std::map<std::string, variable_value>& vm, const std::string& opt, const T& val)
        {
            vm[opt].value() = boost::any(val);
        }

        int register_parsedOptions_to_print(const options_description& descriptions, bool visible=true)
        {
            LOG(MAXDEBUG)<<"register keys ";
            auto vec_opt_shptr = descriptions.options();
            for(const auto& p : vec_opt_shptr)
            {
                std::string key = p->canonical_display_name();
                if(!fVisible_key_map.count(key))
                {
                    LOG(MAXDEBUG)<<"register key "<< key <<" from print_option() at level>=INFO";
                    fVisible_key_map[key] = visible;
                }
                else
                {
                    LOG(MAXDEBUG)<<"key "<< key <<" already registered to print_option() at level>=INFO";
                    
                }
                
            }
            
            return 0;
        }
        
        
        
    private:
        std::shared_ptr<breit_summary> fSummary;
        // /////////////////////////////////////////////
        // Methods below are helper functions used in the print_options method
        typedef std::tuple<std::string, std::string,std::string, std::string> VarValInfo_t;
        typedef std::map<std::string, VarValInfo_t > MapVarValInfo_t;

        VarValInfo_t get_varInfo(const variable_value& var_val);

        template<typename T>
        std::string variable_toString(const variable_value& var_val)
        {
            auto& value = var_val.value();
            std::ostringstream ostr;
            
            if (auto q = boost::any_cast<bool>(&value))
            {
                std::string val_str;
                if(*q)
                    val_str="true";
                else
                    val_str="false";
                return val_str;
            }
            
            if (auto q = boost::any_cast<T>(&value))
            {
                ostr << *q;
            }
            std::string val_str = ostr.str();
            return val_str;
        }

        static void max(int &val, const int &comp)
        {
            if (comp > val)
            {
                val = comp;
            }
        }
    };

}

#endif	/* OPTIONS_MANAGER_H */

