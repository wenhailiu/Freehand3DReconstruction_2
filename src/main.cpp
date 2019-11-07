#include <iostream>

#include "FreehandReconstruction.h"

int main(int argc, char* argv[]){
    FreehandReconstruction reconstruction_handle("/home/wenhai/vsc_workspace/Freehand3DReconstruction_2/init_data/Parameters.yaml");
    reconstruction_handle.ExtractImageToVolumeMatrices();
    reconstruction_handle.MallocHostSpace();
    reconstruction_handle.LaunchGPU_Reconstruction();
    reconstruction_handle.SaveReconstructedVolume();
}