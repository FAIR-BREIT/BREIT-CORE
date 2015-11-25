/* 
 * File:   breit_equations.tpl
 * Author: winckler
 *
 * Created on August 4, 2015, 10:01 AM
 */

namespace breit
{
    /// ////////////////////////////////////////////////////////////////////////////////
    template <typename T, typename U >
    breit_equations<T,U>::breit_equations():  generate_equations<self_type>(), 
                                ui_type(), 
                                fEqDim(0), 
                                fCoefDim(9),
                                fCoef_range_i(),
                                fCoef_range_j(),
                                fCoef_list(),
                                fMat(),
                                f2nd_member(),
                                fF0(), fCoef_index_min(0), fCoef_index_max(0),fIni_cond_map()
            {
            }


    
    /// ////////////////////////////////////////////////////////////////////////////////
    template <typename T, typename U >
    typename breit_equations<T,U>::matrix_d& breit_equations<T,U>::output()
    {
        return fMat;
    }
    /// ////////////////////////////////////////////////////////////////////////////////
    template <typename T, typename U >
    typename breit_equations<T,U>::vector_d& breit_equations<T,U>::snd_member()
    {
        return f2nd_member;
    }
    /// ////////////////////////////////////////////////////////////////////////////////
    template <typename T, typename U >
    int breit_equations<T,U>::parse(const int argc, char** argv, bool AllowUnregistered)
    {
        if(ui_type::parse(argc,argv,AllowUnregistered))
            return 1;

        return 0;
    }
    /// ////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    /// ////////////////////////////////////////////////////////////////////////////////
    // read input and get coefficients of the system
    template <typename T, typename U >
    int breit_equations<T,U>::read_impl()
    {
        LOG(DEBUG)<<"in read_impl() ...";
        ublas::range input_coef_range_i(  fvarmap.at("coef.index.i.min").template as<size_t>() , 
                                          fvarmap.at("coef.index.i.max").template as<size_t>()+1
                                       );
        ublas::range input_coef_range_j(  fvarmap.at("coef.index.j.min").template as<size_t>() , 
                                          fvarmap.at("coef.index.j.max").template as<size_t>()+1
                                       );
        
        

        LOG(DEBUG)<<"init coefficients description ...";
        /// use prog options to parse input data file
        // define options_description of the coefs
        options_description input_file_desc("input file description");
        options_description coef_desc("cross-sections description");
        options_description header_desc("cross-sections header description");
        options_description ic_desc("Initial conditions");
        
        ui_type::init_initial_condition_descriptions(
                                            ic_desc,
                                            fvarmap.at("coef.index.i.min").template as<size_t>() , 
                                            fvarmap.at("coef.index.i.max").template as<size_t>()+1
                                       );
        ui_type::init_coef_descriptions(coef_desc);
        ui_type::init_input_header_descriptions(header_desc);
        input_file_desc.add(header_desc).add(coef_desc).add(ic_desc);
        
        
        
        // parse and fill vm
        variables_map vm;

        // But before parsing must handle UTF8 BOM case...
        path filename=fvarmap["input-file"].template as<path>();
        path filename_out=filename.parent_path();
        filename_out/=filename.stem();
        filename_out+=".breit.txt";

        std::ifstream infile(filename.string());
        

        bool UTF8=false;
        if (infile.good())
        {
            std::string line;
            getline(infile, line);
            if (line.compare( 0, 3, "\xEF\xBB\xBF" ) == 0)
            {
                line.erase( 0, 3 );
                UTF8=true;
            }
            
            if(UTF8)
            {
                std::ofstream outfile(filename_out.string());
                while(getline(infile, line))
                {
                    outfile<<line<<std::endl;
                }        
                outfile.close();
            }
        }

        infile.close();
        
        path file_to_parse;
        if(UTF8)
            file_to_parse=filename_out;
        else
            file_to_parse=filename;

        LOG(MAXDEBUG)<<"parse data file "<<file_to_parse.string()<<" ...";

        if(ui_type::parse_cfgFile(file_to_parse.string(),input_file_desc,vm,false))
            return 1;

         double scale_factor=ui_type::get_scale_factor(vm);
        
        // get vm data into the fCoef_list container 
        fCoef_list.clear();

        size_t index_i_min = std::numeric_limits<size_t>::max(); 
        size_t index_i_max = std::numeric_limits<size_t>::min();
        size_t index_j_min = std::numeric_limits<size_t>::max(); 
        size_t index_j_max = std::numeric_limits<size_t>::min();

        // go over user-provided indices and if the 
        // matrix element Qij is found as non defaulted, 
        // store the indices and value in fCoef_list map
        LOG(DEBUG)<<"searching for coefficients ...";
        
        
        
        for(const auto& i : input_coef_range_i)
            for(const auto& j : input_coef_range_j)
            {
                std::string key=ui_type::form_coef_key(i,j);
                //LOG(DEBUG)<<"Q."<<i<<"."<<j<<" = "<<key;
                if(vm.count(key))
                {
                    //LOG(DEBUG)<<"counted";

                    if(!vm.at(key).defaulted())
                    {
                        LOG(DEBUG)<<"found cross-section coefficient : "<< key <<" = "<< vm.at(key).as<data_type>();

                        std::pair<size_t,size_t> coef_key(i,j);
                        data_type coef_val=vm.at(key).as<data_type>()*scale_factor;
                        fCoef_list.insert( std::make_pair(coef_key, coef_val) );

                        // to resize matrix properly :
                        // get the max/min indices of the coef.
                        // this is necessary to avoid having rows or column full of zeros
                        if(i<index_i_min)
                            index_i_min=i;

                        if(i>index_i_max)
                            index_i_max=i;

                        if(j<index_j_min)
                            index_j_min=j;

                        if(j>index_j_max)
                            index_j_max=j;

                        LOG(MAXDEBUG)<<" i="<< i << " min="<<index_i_min <<" max="<<index_i_max;
                        LOG(MAXDEBUG)<<" j="<< j << " min="<<index_j_min <<" max="<<index_j_max;
                    }
                }
            }

        fCoef_index_min=index_i_min;
        fCoef_index_max=index_i_max;
        // get new range objects (index_i_max + 1 to include last index)
        ublas::range coef_range_i(index_i_min,index_i_max+1);
        ublas::range coef_range_j(index_j_min,index_j_max+1);
        // and update the corresponding private members
        fCoef_range_i=coef_range_i;
        fCoef_range_j=coef_range_j;
        
        // do the same for initial conditions
        index_i_min = std::numeric_limits<size_t>::max(); 
        index_i_max = std::numeric_limits<size_t>::min();
        index_j_min = std::numeric_limits<size_t>::max(); 
        index_j_max = std::numeric_limits<size_t>::min();

        
        
        for(const auto& i : input_coef_range_i)
        {
            std::string key=ui_type::form_init_cond_key(i);
            if(vm.count(key))
            {
                if(!vm.at(key).defaulted())
                {
                    data_type val=vm.at(key).as<data_type>();
                    LOG(DEBUG)<<"found initial conditions : "<< key <<" == "<< val;
                    
                    data_type ic_val=vm.at(key).as<data_type>();
                    fIni_cond_map.insert(std::make_pair(i, ic_val) );
                    if(i<index_i_min)
                            index_i_min=i;

                    if(i>index_i_max)
                        index_i_max=i;
                }
            }
        }
        
        ublas::range init_cond_range(index_i_min,index_i_max+1);
        fEqDim=fvarmap["eq-dim"].template as<size_t>();



        LOG(MAXDEBUG)   <<" generating equations parameters : dim_i = " << fCoef_range_i.size()
                        <<" offset_i = " << fCoef_range_i.start() << " last index_i = " 
                        << fCoef_range_i.size()+fCoef_range_i.start()-1;

        for(const auto& i : fCoef_range_i)
            LOG(MAXDEBUG)<<" range-i="<<i;

        LOG(MAXDEBUG)  <<" generating equations parameters : dim_j = " << fCoef_range_j.size()
                    <<" offset_j = " << fCoef_range_j.start() << " last index_j = " 
                    << fCoef_range_j.size()+fCoef_range_j.start()-1;

        for(const auto& j : fCoef_range_j)
            LOG(MAXDEBUG)<<" range-j="<<j;

        // perform some checks. 
        LOG(DEBUG)<<"checking for index ranges ...";
        
        if(init_cond_range!=fCoef_range_i)
        {
            LOG(ERROR)<<"provided index range of initial conditions does not match the provided index range of coefficient";
            //LOG(ERROR)<<"check the given cross-section coefficients in file "<<filename.string();
            return 1;
        }
        
        if(fCoef_range_i!=fCoef_range_j)
        {
            LOG(ERROR)<<"i and j range index are different and should be the same.";
            LOG(ERROR)<<"check the given cross-section coefficients in file "<<filename.string();
            return 1;
        }
        /*else
        {
            
            if(fEqDim>fCoef_range_i.size())
            {
                LOG(ERROR)  <<"the found matrix dimension (dim="<< fCoef_range_i.size() 
                            <<") is smaller than the user provided one ("<<fEqDim<<").";
                
                LOG(ERROR)  <<"check the provided cross-section coefficients "
                            <<filename.string()
                            <<" or the configuration options ";
                //LOG(WARNING)<<"check the cross-section coefficient indices or the user provided dimension";
                fEqDim=fCoef_range_i.size();
                return 1;
            }
            
            // temporary : force user dim to found dim
            
        }*/

        fEqDim=fCoef_range_i.size();

        // if all checks are fine, fill missing coefs in new reduced range with zeros :
        for(const auto& i : fCoef_range_i)
            for(const auto& j : fCoef_range_j)
            {
                std::pair<size_t,size_t> coef_key(i,j);
                if(!fCoef_list.count(coef_key))
                {
                    fCoef_list.insert( std::make_pair(coef_key, 0) );
                }
            }

        size_t dim=fCoef_range_i.size();
        size_t offset=fCoef_range_i.start();        

        if (ui_type::fUse_cfgFile)
        {
            if (fs::exists(ui_type::fConfig_file_path))
            {
                if (ui_type::parse_cfgFile(ui_type::fConfig_file_path.string(), ic_desc, ui_type::fvarmap, true))
                    return 1;
            }
            else
            {
                LOG(ERROR)<<"config file '"<< fConfig_file_path <<"' not found";
                return 1;
            }

        }
        
        
        vector_d F0(dim); 
        for(const auto& p : fIni_cond_map)
        {
            F0(p.first-offset)=p.second;
        }
        fF0=F0;
        
        
        
        return 0;
    }// end read_impl
    
    
    
    /// ////////////////////////////////////////////////////////////////////////////////
    // transform system of equations into matrix equations
    template <typename T, typename U >
    int breit_equations<T,U>::generate_impl()
    {
        //fEqDim=fvarmap["eq-dim"].template as<int>();
        LOG(DEBUG)<<"generating equations";
        // temporary if(staticeq) :
        bool staticeq=false;

        if(staticeq)
        {
            if(static_eq_system())
                return 1;
        }
        else
        {
            if(dynamic_eq_system())
                return 1;
        }

        std::string verbose=fvarmap["verbose"].template as<std::string>();
        LOG(DEBUG) << "printing generated matrix to process : ";
        LOG(DEBUG) << fMat << std::endl;
        LOG(DEBUG) << "and second member vector : ";
        LOG(DEBUG) << f2nd_member << std::endl;

        return 0;
    }



    /// ////////////////////////////////////////////////////////////////////////////////
    //case dF/dx = MF = 0 with dim(M) = N
    template <typename T, typename U >
    int breit_equations<T,U>::static_eq_system()
    {
        size_t dim=fCoef_range_i.size();
        size_t offset=fCoef_range_i.start();
        size_t last_index=offset+dim-1;

        matrix_d mat(dim,dim);
        for(const auto& p : fCoef_range_i)
            for(const auto& q : fCoef_range_j)
                if(last_index != p)
                    mat(p-offset,q-offset)=compute_matrix_element(p,q);
                else
                    mat(dim-1,q-offset)=1.0;

        // clear and copy
        fMat.clear();
        fMat=mat;

        return 0;
    }


    /// ////////////////////////////////////////////////////////////////////////////////
    // case dF/dx = AF + cte with dim(A) = N-1
    // note :   we could eventually use this method as well for dF/dx = 0, 
    //          and works in dim = N-1 since condition 
    //          FN=1-Sum(k=1 to N-1) Fk is satisfied.
    template <typename T, typename U >
    int breit_equations<T,U>::dynamic_eq_system()
    {

        size_t dim=fCoef_range_i.size();
        size_t offset=fCoef_range_i.start();
        size_t last_index=offset+dim-1;
        fEqDim=dim;// todo : clean previous assignment
        fSummary->system_dim=dim;
        fSummary->reduced_system_dim=dim-1;
        fSummary->offset=offset;
        LOG(DEBUG)   <<"in dynamic_eq_system() equations parameters : dim = " << dim
                        <<" offset = " << offset << " last index = " << last_index;
        ublas::range reduced_range(offset,last_index);// last point excluded : dim=N-1
        if(dim>1)
        {
            // store and map indexes
            for(const auto& p : fCoef_range_i)
                if(p-offset>=0)
                    fSummary->F_index_map[p-offset] = static_cast<int>(p);
                else
                {
                    LOG(ERROR)<<"index difference with the offset is negative";
                    return 1;
                }
            LOG(DEBUG)<<"size of map index "<< fSummary->F_index_map.size();
            for(const auto& p : fSummary->F_index_map)
                LOG(DEBUG)<<" index " << p.first << " -> "<<p.second;

            // system is reduced by one dimension due to condition FN=1-Sum(k=1 to N-1) Fk
            matrix_d mat(dim-1,dim-1);
            for(const auto& p : reduced_range)
                for(const auto& q : reduced_range)
                {
                    LOG(MAXDEBUG)<<" coef("<<p<<","<<q <<")"<< " offset="<<offset;
                    LOG(MAXDEBUG)<<" mat("<<p-offset<<","<<q-offset <<")"<< " offset="<<offset;
                    
                    mat(p-offset,q-offset)=compute_matrix_element(p,q)-compute_matrix_element(p,last_index);
                    LOG(MAXDEBUG)<<"mat("<<p<<","<<q <<") = "<< mat(p-offset,q-offset);
                }

            // clear and copy
            fMat.clear();
            fMat=mat;
            vector_d Cte(reduced_range.size());
            for(const auto& p : reduced_range)
            {
                Cte(p-offset)=compute_matrix_element(p,last_index);
                LOG(MAXDEBUG)<<"p-offset="<<p-offset << " <-> v("<< p <<")="<<compute_matrix_element(p,last_index);
            }
            //std::cout<<Cte<<std::endl;
            // clear and copy
            f2nd_member.clear();
            f2nd_member=Cte;
        }

        return 0;
    }
    
    
    /// ////////////////////////////////////////////////////////////////////////////////
    // temporary
    template <typename T, typename U >
    std::vector<double> breit_equations<T,U>::get_1electron_approximation_solution()
    {
        LOG(MAXDEBUG)<<"1-electron approximation solution- matrix dim = "<<fMat.size1();
        
        size_t syst_dim=fMat.size1()+1;
        size_t coef_index_min=fCoef_index_min-1;
        
        //F1
        typedef std::pair<size_t,size_t> coef;
        std::vector<double> ana_sol;
        //data_type coef_val=vm.at(coef(14,15)).as<data_type>();
        //*
        double denominator=1+fCoef_list.at(coef(coef_index_min+syst_dim-1,coef_index_min+syst_dim))/fCoef_list.at(coef(coef_index_min+syst_dim,coef_index_min+syst_dim-1));
        //double denominator=1+fQ[14][15]/fQ[15][14];
        //std::cout<<"denominator (1) = "<<denominator<<std::endl;

        for(int i(coef_index_min+syst_dim-2);i>coef_index_min;i--)
        {
            //denominator*=fQ[i][i+1]/fQ[i+1][i];
            denominator*=fCoef_list.at(coef(i,i+1))/fCoef_list.at(coef(i+1,i));
            denominator+=1.0;
            //std::cout<<"denominator ("<<i<<") = "<<denominator<<std::endl;
        }


        //std::cout<<"denominator="<<denominator<<std::endl;
        double F1=1/denominator;
        ana_sol.push_back(F1);
        //std::cout<<"denominator="<<denominator<<std::endl;
        LOG(MAXDEBUG)<<"F1="<<F1;
        double Fi=F1;
        double Fip1=0.0;
        double sum=F1;
        for(int i(coef_index_min+1);i<coef_index_min+syst_dim;i++)
        {
            LOG(MAXDEBUG)<<"i="<<i;
            //Fip1=Fi*(fQ[i][i+1]/fQ[i+1][i]);
            Fip1=Fi*(fCoef_list.at(coef(i,i+1))/fCoef_list.at(coef(i+1,i)));
            LOG(MAXDEBUG)<<"F"<<i+1<<"="<<Fip1;
            ana_sol.push_back(Fip1);
            Fi=Fip1;
            sum+=Fip1;

        }
        ana_sol.push_back(sum);
        LOG(MAXDEBUG)<<"sum="<<sum;
        
        
        //*/
        return ana_sol;
    }

}//end namespace breit