/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
/* 
 * File:   logger_def.h
 * Author: winckler
 *
 * Created on September 2, 2015, 11:58 AM
 */

#ifndef LOGGER_DEF_H
#define	LOGGER_DEF_H

#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/log/utility/manipulators/to_log.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace breit
{
    enum severity_level
    {
        MAXDEBUG,
        DEBUG,
        RESULTS,
        INFO,
        WARN,
        ERROR,
        STATE,
        NOLOG
    };        


    static const std::array<std::string, 8> g_LogSeverityLevelString
    {
        { 
            "[MAXDEBUG]",
            "[DEBUG]   ", 
            "[RESULTS] ",
            "[INFO]    ", 
            "[WARN]    ", 
            "[ERROR]   ", 
            "[STATE]   ", 
            "[NOLOG]   " 
        }
    };



    namespace color 
    {
        enum code 
        {
            FG_BLACK    = 30,
            FG_RED      = 31,
            FG_GREEN    = 32,
            FG_YELLOW   = 33,
            FG_BLUE     = 34,
            FG_MAGENTA  = 35,
            FG_CYAN     = 36,
            FG_WHITE    = 37,
            FG_DEFAULT  = 39,
            BG_RED      = 41,
            BG_GREEN    = 42,
            BG_BLUE     = 44,
            BG_DEFAULT  = 49
        };
    }
}

// helper function to format in color console output
template <breit::color::code color> 
inline std::string write_in(const std::string& text_in_bold)
{
    std::ostringstream os;
    os << "\033[01;" << color << "m" << text_in_bold << "\033[0m";

    return os.str();
}

// typedef
typedef breit::severity_level custom_severity_level;
#define SEVERITY_THRESHOLD custom_severity_level::MAXDEBUG
#define SEVERITY_ERROR custom_severity_level::ERROR
#define SEVERITY_NOLOG custom_severity_level::NOLOG
// tags used for log console or file formatting
struct tag_console;
struct tag_file;


// overload operator for console output
inline boost::log::formatting_ostream& operator<<
        (
            boost::log::formatting_ostream& strm,
            boost::log::to_log_manip< custom_severity_level, tag_console > const& manip
        )
{
    custom_severity_level level = manip.get();
    std::size_t idx=static_cast< std::size_t >(level);
    if ( idx < breit::g_LogSeverityLevelString.size() )
    {
        //strm <<" idx = "<<idx <<" ";
        switch (level)
        {
            case custom_severity_level::MAXDEBUG :
                strm << write_in<breit::color::FG_BLUE>(breit::g_LogSeverityLevelString.at(idx));
                break;
            
            case custom_severity_level::DEBUG :
                strm << write_in<breit::color::FG_BLUE>(breit::g_LogSeverityLevelString.at(idx));
                break;

            case custom_severity_level::RESULTS :
                strm << write_in<breit::color::FG_MAGENTA>(breit::g_LogSeverityLevelString.at(idx));
                break;

            case custom_severity_level::INFO :
                strm << write_in<breit::color::FG_GREEN>(breit::g_LogSeverityLevelString.at(idx));
                break;    
                
            case custom_severity_level::WARN :
                strm << write_in<breit::color::FG_YELLOW>(breit::g_LogSeverityLevelString.at(idx));
                break;

            case custom_severity_level::STATE :
                strm << write_in<breit::color::FG_MAGENTA>(breit::g_LogSeverityLevelString.at(idx));
                break;    
                
            case custom_severity_level::ERROR :
                strm << write_in<breit::color::FG_RED>(breit::g_LogSeverityLevelString.at(idx));
                break;    
                
            case custom_severity_level::NOLOG :
                strm << write_in<breit::color::FG_DEFAULT>(breit::g_LogSeverityLevelString.at(idx));
                break;
            
            default:
                break;
        }
        
    }
    else
    {
        strm    << write_in<breit::color::FG_RED>("Unknown log level ")
                << "(int level = "<<static_cast< int >(level)
                <<")";
    }
    return strm;
}

// overload operator for file output
inline boost::log::formatting_ostream& operator<<
        (
            boost::log::formatting_ostream& strm,
            boost::log::to_log_manip< custom_severity_level, tag_file > const& manip
        )
{
    custom_severity_level level = manip.get();
    std::size_t idx=static_cast< std::size_t >(level);
    if ( idx < breit::g_LogSeverityLevelString.size() )
        strm << breit::g_LogSeverityLevelString.at(idx);
    else
    {
        strm    << write_in<breit::color::FG_RED>("Unknown log level ")
                << "(int level = "<<static_cast< int >(level)
                <<")";
    }
    return strm;
}




// other helper functions

namespace breit
{

    
    template <typename T>
    std::string to_string_scientific(const T a_value)
    {
        std::ostringstream out;
        out << std::scientific << a_value;
        return out.str();
    }
    
    
    
    template<typename charT, typename traits = std::char_traits<charT> >
    class bstream_center_helper 
    {
        std::basic_string<charT, traits> str_;
    public:
        bstream_center_helper(std::basic_string<charT, traits> str) : str_(str) {}
        template<typename a, typename b>
        friend std::basic_ostream<a, b>& operator<<(std::basic_ostream<a, b>& s, const bstream_center_helper<a, b>& c);
    };

    template<typename charT, typename traits = std::char_traits<charT> >
    inline bstream_center_helper<charT, traits> bstream_centered(std::basic_string<charT, traits> str) {
        return bstream_center_helper<charT, traits>(str);
    }

    // redeclare for std::string directly so we can support anything that implicitly converts to std::string
    inline bstream_center_helper<std::string::value_type, std::string::traits_type> bstream_centered(const std::string& str) {
        return bstream_center_helper<std::string::value_type, std::string::traits_type>(str);
    }

    template<typename charT, typename traits>
    inline std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& s, const bstream_center_helper<charT, traits>& c) {
        std::streamsize w = s.width();
        if (w > c.str_.length()) {
            std::streamsize left = (w + c.str_.length()) / 2;
            s.width(left);
            s << c.str_;
            s.width(w - left);
            s << "";
        } else {
            s << c.str_;
        }
        return s;
    }
}

#endif	/* LOGGER_DEF_H */

