#ifndef US_FRAMES
#define US_FRAMES

#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"

#include "img_base.h"

namespace ImageBase
{

struct us_parameters_structure{
    ImageBase::Dim_pxl dim_pxl;
    ImageBase::PixelSize_mm p_size_mm;
    ImageBase::SpatialOrigin origin_mm;
};

class Ultrasound_Container{

public:
    Ultrasound_Container(std::string Input_Parameter_Path);

protected:

    //Member functions:
    //contents setup:
    bool setup_OnetimeReconstruction();
    //TODO:
    bool setup_ProgressiveReconstruction();
    bool add_One_FrameWithMatrix(std::vector<uint8_t> InputFrame, std::vector<std::vector<float>> InputProbeToWorld);


    //Image parameters:
    us_parameters_structure us_parameters;

    //Frames info:
    int NumFrames;
    std::string frame_data_path;
    std::vector<uint8_t> Frames_Raw_Data;

    //Transform data:
    std::string ProbeTracker_list_path;
    int ProbeTrackingNumber;
    std::vector<float> ImageToProbeMatrix;
    std::string ProbeCalibrationDate;

    //Tracking data container:
    std::vector<std::vector<float>> ProbeToWorldMatrices; //used for onetime_reconstruction
    std::vector<std::vector<float>> ProbeToWorldMatrix; //used for import matrices stream, both for ontime and progressively reconstruct. 
    // std::vector<std::vector<float>> ReferenceToWorldMatrices;

private:
    //yaml parser handle:
    YAML::Node parameters_file_handle;

    //Flags: true, indicates ERROR happended; 
    bool US_ERROR; //should be always false;
};

} 


#endif