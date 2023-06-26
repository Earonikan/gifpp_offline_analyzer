#include "InputArguments.hpp"

InputArguments::InputArguments(int argv, char *argc[])
    {
        if (argv == 1)
        {
            std::cout << "You need to use arguments!" << std::endl;
            exit(0);
        }
        std::string arg;
        ParameterParser parser1(parameters.argKeys);
        ParameterParser parser2(parameters.runTypes);
        for (int i = 1; i < argv; ++i)
        {
            arg = argc[i];
            switch (parser1.ParseThis(arg))
            {
                case 1:
                {
                    parameters.runkey = argc[i+1];
                    if (parser2.ParseThis(parameters.runkey) > 0)
			        {
                        if (parameters.verbose > 0) std::cout << "Using runkey: " << parameters.runkey << std::endl;
                    }
                    else
                    {
                        std::cout << "Wrong runkey provided!!!" << std::endl;
                        exit(0);
                    }
                    break;
                }

                case 2:
                {
                    parameters.runlist_filename = argc[i+1]; break;
                }

                case 3: 
                {
                    parameters.runstart = atoi(argc[i+1]);
                    parameters.runstop = parameters.runstart;
		            break;
                }

                case 4:
                {
                    parameters.runstart = atoi(argc[i+1]);
                    parameters.runstop = atoi(argc[i+2]);
	    	        break;
                }

                case 5:
                {
                    parameters.stdout_flag = true; break;
                }
		
		        case 6:
                {
                    parameters.verbose = atoi(argc[i+1]); break;
                }

                case 7:
                {
                    parameters.root_tree = true; break;
                }

            } 
        }
    }
