#ifndef FREEHANDRECON
#define FREEHANDRECON

#include "../utilities/ultrasound_frames.h"
#include "../utilities/volumetric.h"

class FreehandReconstruction: private ImageBase::Ultrasound_Container, private ImageBase::Volume_Container{

public:
    FreehandReconstruction(std::string Initialization_parameter_path);
    FreehandReconstruction(std::string Initialization_parameter_path, const float StartAngle, const float EndAngle);
    bool ExtractImageToVolumeMatrices(int smoothPatch = 5);
    bool MallocHostSpace();
    bool LaunchGPU_Reconstruction();
    bool SaveReconstructedVolume();

    ~FreehandReconstruction();

private:

    //Extracted Matrices:
    int MatchedMatrixNumber;
    std::vector<std::vector<float>> ImageToVolumeMatrices;
    std::vector<std::vector<float>> TiltedImageMatrices;

    //scan parameters:
    float angleBegin; 
    float angleEnd;
    float Motor_Radius;
    float TrimThreshold;

    //Flags:
    bool RECON_ERROR; //should be always false;
    bool RawPtrSet;

    //yaml handler:
    YAML::Node Reconstructor_ParamHandle;

    //raw data pointers:
    float *ImageToVolume_hostPtr;

};

#endif