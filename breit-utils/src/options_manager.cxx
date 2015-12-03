/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

/*
 * File:   options_manager.cxx
 * Author: winckler
 * 
 * Created on July 13, 2015, 10:30 PM
 */
// refactored from FairRoot::FairMQ project

#include "options_manager.h"

namespace breit
{

    /// //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor / destructor
    options_manager::options_manager() : 
                            fGenericDesc("Generic options description"), 
                            fConfigDesc("Configuration options description"), 
                            fHiddenDesc("Hidden options description"), 
                            fEnvironmentDesc("Environment Variables"),
                            fCmdline_options("Command line options"), 
                            fConfig_file_options("Configuration file options"), 
                            fVisible_options("User options"),
                            fVerboseLvl("RESULTS"), fUse_cfgFile(false), fConfig_file_path(), fVisible_key_map()
    {
        // //////////////////////////////////
        // define generic options
        fGenericDesc.add_options()
            ("help,h", "produce help")
            ("version,v", "print version")
            ("verbose", po::value<std::string>(&fVerboseLvl)->default_value("INFO"), "Verbosity level : \n"
                                                                        "  0=DEBUG \n"
                                                                        "  1=INFO \n"
                                                                        "  2=WARN \n"
                                                                        "  3=ERROR \n"
                                                                        "  4=STATE \n"
                                                                        "  5=NOLOG"
                )
            ("log-color", po::value<bool>()->default_value(true), "logger color format : true or false")
            ;

    }

    options_manager::~options_manager() 
    {
    }


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Add option descriptions

    int options_manager::addTo_cmdLine(const options_description& optdesc, bool visible)
    {
        fCmdline_options.add(optdesc);
        if(visible)
            fVisible_options.add(optdesc);
        register_parsedOptions_to_print(optdesc,visible);

        return 0;
    }

    int options_manager::addTo_cfgFile(const options_description& optdesc, bool visible)
    {
        //if use_cfgFile() not yet called, then enable it with required file name to be provided by command line
        if(!fUse_cfgFile)
            use_cfgFile();

        fConfig_file_options.add(optdesc);
        if(visible)
            fVisible_options.add(optdesc);
        register_parsedOptions_to_print(optdesc,visible);
        return 0;
    }

    int options_manager::addTo_env(const options_description& optdesc)
    {
        fEnvironmentDesc.add(optdesc);
        return 0;
    }


    void options_manager::use_cfgFile(const std::string& filename)
    {
            fUse_cfgFile = true;
            if (filename.empty())
            {
                fConfigDesc.add_options()
                    ("config,c", po::value<path>(&fConfig_file_path)->required(), "Path to configuration file (required argument)");
                addTo_cmdLine(fConfigDesc);
            }
            else
            {
                fConfig_file_path = filename;
            }
    }

    /// //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Parser

    int options_manager::parse_cmdLine(const int argc, char** argv, const options_description& desc, variables_map& varmap, bool AllowUnregistered)
    {
        try
        {
            // //////////////////////////////////
            // get options from cmd line and store in variable map
            // here we use command_line_parser instead of parse_command_line, to allow unregistered and positional options
            if(AllowUnregistered)
            {
                po::command_line_parser parser{argc, argv};
                parser.options(desc).allow_unregistered();
                po::parsed_options parsedOptions = parser.run();
                po::store(parsedOptions,varmap);
            }
            else
                po::store(po::parse_command_line(argc, argv, desc), varmap);

            // //////////////////////////////////
            // call the virtual notifySwitch_options method to handle switch options like e.g. "--help" or "--version"
            // return 1 if switch options found in varmap
            if(notifySwitch_options())
                return 1;

            po::notify(varmap);
        }
        catch(std::exception& e)
        {
            LOG(ERROR) << e.what();
            return 1;
        }
        return 0;
    }


    int options_manager::parse_cmdLine(const int argc, char** argv, const options_description& desc, bool AllowUnregistered)
    {
        return parse_cmdLine(argc,argv,desc,fvarmap,AllowUnregistered);
    }



    int options_manager::parse_cfgFile(std::ifstream& ifs, const options_description& desc, variables_map& varmap, bool AllowUnregistered)
    {
        try
        {
            if (!ifs)
            {
                LOG(ERROR) << "can not open configuration file \n";
                return -1;
            }
            else
            {
                po:store(parse_config_file(ifs, desc, AllowUnregistered), varmap);
                po::notify(varmap);
            }
        }
        catch(std::exception& e)
        {
            LOG(ERROR) << e.what();
            return 1;
        }
        return 0;
    }

    int options_manager::parse_cfgFile(const std::string& filename, const options_description& desc, variables_map& varmap, bool AllowUnregistered)
    {
        try
        {
            if (!fs::exists(filename))
            {
                LOG(ERROR)<<"file '"<< filename <<"' not found";
                return 1;
            }

            std::ifstream ifs(filename.c_str());
            if (!ifs)
            {
                LOG(ERROR) << "can not open file: " << filename <<"'";
                return 1;
            }
            else
            {
                po:store(parse_config_file(ifs, desc, AllowUnregistered), varmap);
                po::notify(varmap);
            }
        }
        catch(std::exception& e)
        {
            LOG(ERROR) << e.what();
            return 1;
        }
        
        return 0;
    }


    int options_manager::parse_cfgFile(const std::string& filename, const options_description& desc, bool AllowUnregistered)
    {
        return parse_cfgFile(filename,desc,fvarmap,AllowUnregistered);
    }

    int options_manager::parse_cfgFile(std::ifstream& ifs, const options_description& desc, bool AllowUnregistered)
    {
        return parse_cfgFile(ifs,desc,fvarmap,AllowUnregistered);
    }


    int options_manager::parse_environment(const std::function<std::string(std::string)> & environment_mapper)
    {
        po::store(po::parse_environment(fEnvironmentDesc, environment_mapper), fvarmap);
        po::notify(fvarmap);

        return 0;
    }

    // Given a key, convert the variable value to string
    std::string options_manager::get_stringVal(const std::string& key)
        {
            std::string val_str;
            try
            {
                if ( fvarmap.count(key) )
                {
                    auto& value = fvarmap[key].value();


                    // string albeit useless here
                    if(auto q = boost::any_cast< std::string >(&value ))
                    {
                        val_str = variable_toString< std::string >(fvarmap[key]);
                        return val_str;
                    }

                    // vector<string>
                    if(auto q = boost::any_cast< std::vector<std::string> >(&value ))
                    {
                        val_str = variable_toString< std::vector<std::string> >(fvarmap[key]);
                        return val_str;
                    }

                    // int
                    if(auto q = boost::any_cast< int >(&value ))
                    {
                        val_str = variable_toString< int >(fvarmap[key]);
                        return val_str;
                    }

                    // vector<int>
                    if(auto q = boost::any_cast< std::vector<int> >(&value ))
                    {
                        val_str = variable_toString< std::vector<int> >(fvarmap[key]);
                        return val_str;
                    }

                    // size_t
                    if(auto q = boost::any_cast< size_t >(&value ))
                    {
                        val_str = variable_toString< size_t >(fvarmap[key]);
                        return val_str;
                    }

                    // vector<size_t>
                    if(auto q = boost::any_cast< std::vector<size_t> >(&value ))
                    {
                        val_str = variable_toString< std::vector<size_t> >(fvarmap[key]);
                        return val_str;
                    }
                    
                    // float
                    if(auto q = boost::any_cast< float >(&value ))
                    {
                        val_str = variable_toString< float >(fvarmap[key]);
                        return val_str;
                    }

                    // vector float
                    if(auto q = boost::any_cast< std::vector<float> >(&value ))
                    {
                        val_str = variable_toString< std::vector<float> >(fvarmap[key]);
                        return val_str;
                    }

                    // double
                    if(auto q = boost::any_cast< double >(&value ))
                    {
                        val_str = variable_toString< double >(fvarmap[key]);
                        return val_str;
                    }

                    // vector double
                    if(auto q = boost::any_cast< std::vector<double> >(&value ))
                    {
                        val_str = variable_toString< std::vector<double> >(fvarmap[key]);
                        return val_str;
                    }
                    
                    
                    if(auto q = boost::any_cast< bool >(&value ))
                    {
                        val_str = variable_toString< bool >(fvarmap[key]);
                        
                        
                        return val_str;
                    }


                }
            }
            catch(std::exception& e)
            {
                LOG(ERROR) << "Exception thrown for the key '" << key << "'";
                LOG(ERROR) << e.what();
            }

            return val_str;
        }


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Print/notify options

    int options_manager::print_help()  const
    {
        std::cout << fVisible_options << "\n";
        return 0;
    }

    int options_manager::print_options()
    {
        // //////////////////////////////////
        // Method to overload.
        // -> loop over variable map and print its content
        // -> In this example the following types are supported:
        // std::string, int, float, double, path
        // std::vector<std::string>, std::vector<int>, std::vector<float>, std::vector<double>


        MapVarValInfo_t mapinfo;

        // get string length for formatting and convert varmap values into string
        int maxlen_1st = 0;
        int maxlen_2nd = 0;
        int maxlen_TypeInfo = 0;
        int maxlen_default =0;
        int maxlen_empty = 0;
        int total_len=0;
        for (const auto& m : fvarmap)
        {
            max(maxlen_1st, m.first.length());

            VarValInfo_t valinfo=get_varInfo(m.second);
            mapinfo[m.first]=valinfo;
            std::string val_str;
            std::string typeInfo_str;
            std::string default_str;
            std::string empty_str;
            std::tie(val_str,typeInfo_str,default_str,empty_str)=valinfo;

            max(maxlen_2nd, val_str.length());
            max(maxlen_TypeInfo, typeInfo_str.length());
            max(maxlen_default, default_str.length());
            max(maxlen_empty, empty_str.length());

        }

        // TODO : limit the value length field in a better way
        if(maxlen_2nd>100) 
            maxlen_2nd=100;
        total_len=maxlen_1st+maxlen_2nd+maxlen_TypeInfo+maxlen_default+maxlen_empty;


        //maxlen_2nd=200;

        // formatting and printing

        LOG(INFO)<<std::setfill ('*') << std::setw (total_len+3)<<"*";// +3 because of string " = "
        std::string print_optionsTitle="     Program options found     ";

        int leftSpace_len=0;
        int rightSpace_len=0;
        int leftTitle_shift_len=0;
        int rightTitle_shift_len=0;

        leftTitle_shift_len=print_optionsTitle.length()/2;

        if ((print_optionsTitle.length()) % 2)
            rightTitle_shift_len=leftTitle_shift_len+1;
        else
            rightTitle_shift_len=leftTitle_shift_len;

        leftSpace_len=(total_len+3)/2-leftTitle_shift_len;
        if ((total_len+3) % 2) 
            rightSpace_len=(total_len+3)/2-rightTitle_shift_len+1;
        else
            rightSpace_len=(total_len+3)/2-rightTitle_shift_len;


        LOG(INFO) <<std::setfill ('*') << std::setw(leftSpace_len) <<"*"
                    <<std::setw(print_optionsTitle.length()) 
                << std::setfill(' ')<< print_optionsTitle 
                    <<std::setfill ('*') << std::setw(rightSpace_len) <<"*";

        LOG(INFO) <<std::setfill ('*') << std::setw (total_len+3)<<"*";

        for (const auto& p : mapinfo)
        {
            std::string key_str;
            std::string val_str;
            std::string typeInfo_str;
            std::string default_str;
            std::string empty_str;
            key_str=p.first;
            std::tie(val_str,typeInfo_str,default_str,empty_str)=p.second;
            
            std::ostringstream os;
            os  << std::setfill (' ')
                << std::setw(maxlen_1st)
                << std::left 
                << p.first << " = " 
                << std::setw(maxlen_2nd) 
                << val_str 
                << std::setw(maxlen_TypeInfo) 
                << typeInfo_str 
                << std::setw(maxlen_default)
                << default_str 
                << std::setw(maxlen_empty)
                << empty_str;
            
            if(fVisible_key_map.count(key_str))
            {
                if(fVisible_key_map.at(key_str))
                    LOG(INFO) << os.str();
                else
                    LOG(DEBUG) << os.str();
            }
            else
            {
                LOG(WARN) << os.str();// for the unregistered option keys
            }
            
        }
        LOG(INFO)<<std::setfill ('*') << std::setw (total_len+3)<<"*";// +3 for " = "
        return 0;
    }



    int options_manager::notifySwitch_options()
    {
        // //////////////////////////////////
        // Method to overload.
        if ( fvarmap.count("help") )
        {
            std::cout << "***** FAIR Program Options ***** \n" << fVisible_options;
            return 1;
        }

        if (fvarmap.count("version")) 
        {
            std::cout << "alpha version 0.0\n";
            return 1;
        }

        return 0;
    }

    options_manager::VarValInfo_t options_manager::get_varInfo(const variable_value& var_val)
        {
            // tuple<val_str, type_info_str, default_str, empty>
            auto& value = var_val.value();
            std::string defaulted_val;
            std::string empty_val;

            if(var_val.empty())
                empty_val="  [empty]";
            else
                if(var_val.defaulted())
                    defaulted_val="  [default value]";
                else
                    defaulted_val="  [provided value]";

            empty_val+=" *";
            // string
            if(auto q = boost::any_cast< std::string >(&value))
            {
                std::string val_str = *q;
                return std::make_tuple(val_str,std::string("  [Type=string]"),defaulted_val,empty_val);
            }

            // vector<string>
            if(auto q = boost::any_cast< std::vector<std::string> >(&value))
            {
                std::string val_str = variable_toString< std::vector<std::string> >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=vector<string>]"),defaulted_val,empty_val);
            }

            // int
            if(auto q = boost::any_cast< int >(&value))
            {
                std::string val_str =  variable_toString< int >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=int]"),defaulted_val,empty_val);
            }

            // vector<int>
            if(auto q = boost::any_cast< std::vector<int> >(&value))
            {
                std::string val_str = variable_toString< std::vector<int> >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=vector<int>]"),defaulted_val,empty_val);
            }

            // float
            if(auto q = boost::any_cast< float >(&value))
            {
                std::string val_str = variable_toString< float >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=float]"),defaulted_val,empty_val);
            }

            // vector<float>
            if(auto q = boost::any_cast< std::vector<float> >(&value))
            {
                std::string val_str = variable_toString< std::vector<float> >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=vector<float>]"),defaulted_val,empty_val);
            }

            // double
            if(auto q = boost::any_cast< double >(&value))
            {
                std::string val_str = variable_toString< double >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=double]"),defaulted_val,empty_val);
            }

            // vector<double>
            if(auto q = boost::any_cast< std::vector<double> >(&value))
            {
                std::string val_str = variable_toString< std::vector<double> >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=vector<double>]"),defaulted_val,empty_val);
            }

            // path
            if(auto q = boost::any_cast<path>(&value))
            {
                std::string val_str = (*q).string();
                //std::string val_str = (*q).filename().generic_string();
                return std::make_tuple(val_str,std::string("  [Type=path]"),defaulted_val,empty_val);
            }
            
            if(auto q = boost::any_cast<std::size_t>(&value))
            {
                std::string val_str = variable_toString< std::size_t >(var_val);
                //std::string val_str = (*q).filename().generic_string();
                return std::make_tuple(val_str,std::string("  [Type=size_t]"),defaulted_val,empty_val);
            }

            if(auto q = boost::any_cast< bool >(&value))
            {
                std::string val_str = variable_toString< bool >(var_val);
                return std::make_tuple(val_str,std::string("  [Type=bool]"),defaulted_val,empty_val);
            }
            
            // if we get here, the type is not supported return unknown info
            return std::make_tuple(std::string("Unknown value"), std::string("  [Type=Unknown]"),defaulted_val,empty_val);
        }


}