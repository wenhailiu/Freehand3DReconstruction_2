#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

#include "FreehandReconstruction.h"

//Arguments parser:
class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    /// @author iain
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }
    /// @author iain
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
                != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

void PrintUsage(std::ostream& os){
    os << "Usage: " << std::endl;
    os << std::setw(15) << "-h:    Print usages. " << std::endl;
    os << std::setw(15) << "-c:    Load configuration file, which is necessary for image reconstruction. " << std::endl;
    os << std::setw(15) << "-m:    Data export mode: (default) I for integer (with hole-filling); F for float (witout hole-filling). " << std::endl;
}

int main(int argc, char* argv[]){

    InputParser input(argc, argv);

    if(input.cmdOptionExists("-h") || argc == 1){
        PrintUsage(std::cout);
        return 0;
    }

    std::string fileName;
    if(input.cmdOptionExists("-c")){
        fileName = input.getCmdOption("-c");

        if(fileName.empty()){
            PrintUsage(std::cout);
            return 0;
        }
    }
    else{
        PrintUsage(std::cout);
        return 0;
    }

    std::string Mode; 
    if(input.cmdOptionExists("-m")){
        Mode = input.getCmdOption("-m");

        if(Mode.empty()){
            std::cout << "Default mode: I is used. " << std::endl; 
            Mode = "I"; 
        }
    }
    else{
        std::cout << "Default mode: I is used. " << std::endl; 
        Mode = "I"; 
    }

    char mode_char = Mode[0]; 

    FreehandReconstruction reconstruction_handle(fileName);
    reconstruction_handle.ExtractImageToVolumeMatrices();
    reconstruction_handle.MallocHostSpace();
    reconstruction_handle.LaunchGPU_Reconstruction();
    reconstruction_handle.SaveReconstructedVolume(mode_char);
}