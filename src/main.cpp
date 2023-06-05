#include "stdafx.hpp"
// #include "Fitter.hpp"
#include "Analyzer.hpp"


int main(int argv, char *argc[]) // 
{
    std::string arg, key, runlist_filename = "/home/qfl/online/midas_digi/daq/runs_history.txt";
    int runstart = 0, runstop = 0;

    for (int i = 1; i < argv-1; ++i)
    {
        arg = argc[i];
        if (arg == "--type")
        {
            key = argc[i+1];
            if ((key == "BVSCAN") || (key == "LED") || (key == "ALL")) std::cout << "Using key: " << key << std::endl;
            else
            {
                std::cout << "Wrong key parsed!!!" << std::endl;
                exit(0);
            }
        }
        if (arg == "--filename") runlist_filename = argc[i+1];
        
    }

    Analyzer analyzer(key, runlist_filename, runstart, runstop);





    // const char *ext = ".root";
    // const char *ext = "ch1.png";
    // TSystemDirectory dir(".", "/home/qfl/online/midas_digi/daq");
    // TList *files = dir.GetListOfFiles();
    
    // int num1, num2;
    // if (files) {

    //     TSystemFile *file;
    //     TString fname;
    //     TIter next(files);

    //     while ((file=(TSystemFile*)next())) {
    //         fname = file->GetName();
    //         if (!file->IsDirectory() && fname.EndsWith(ext)) {
    //             // std::cout << fname.Data() << std::endl;
    //             // sscanf(fname.Data(), "BVscan_run%d_run%d_ch1.png", &num1, &num2);
    //             // if ((num1>1576))
    //             // {
    //                 // std::cout << "Processing BV scan, started from run " << num1+1 << " finished with run " << num2;// << std::endl;
    //                 // auto start = std::chrono::steady_clock::now();
    //                 // fitBVscan(num1+1, num2);
    //                 // auto end = std::chrono::steady_clock::now();
    //                 // std::cout << ", elapsed time in seconds: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " sec." << std::endl;
    //             // }
    //             sscanf(fname.Data(), "output0000%d.root", &num1);
    //             if ((num1 >= 4204)) //&&(num1 < 3531))
    //             {
    //                 std::cout << "Processing run " << num1 << std::endl;
    //                 // fitGain(num1);


    //             }
    //         }
    //     }

    //     delete file;
    // }

    // delete files;

    return 0;
}

