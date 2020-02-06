#ifndef VOLUMETRIC
#define VOLUMETRIC

#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"

#include "img_base.h"

namespace ImageBase
{

struct volume_parameters_structure{
    ImageBase::Dim_vxl dim_vxl;
    ImageBase::VoxelSize_mm v_size_mm; 
    ImageBase::SpatialOrigin origin_mm;
};

class Volume_Container{

public:
    Volume_Container(std::string Input_Parameter_Path); 
    bool IsVolumeDataValid();
    bool IsVolumeDataReconstructed();
    ~Volume_Container();

protected:

    //Member functions:
    void SetVolumeDataValid();
    void SetVolumeDataReconstructed();
    
    //contents initialization:
    bool setup_VolumeToReconstruction(); 
    bool setup_ImportReconstruction(); 

    //Volumetric image patameter:
    volume_parameters_structure vol_parameters;

    //Output path:
    std::string Output_Path;

    //Tracking data container:
    std::string ReferenceTracker_list_path;
    int ReferenceTrackingNumber;
    std::vector<std::vector<float>> ReferenceToWorldMatrices;
    std::vector<std::vector<float>> ReferenceToWorldMatrix;

    //volume data:
    std::vector<float> Volume_data;
    std::vector<float> Weights_data;
    std::vector<uint8_t> VolumeToSave; 

    // Img_Plane XY_Plane
    std::vector<float> XY_Plane;
    std::vector<float> XZ_Plane;
    std::vector<float> YZ_Plane;

private:
    std::string Input_Parameter_Path_internal;
    //yaml parser handle:
    YAML::Node parameters_file_handle;

    //Flags:
    bool VOL_ERROR; //should be always false;
    bool Data_Valid;
    bool Data_Reconstructed;
};

}

#endif