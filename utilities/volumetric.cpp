#include <iostream>
#include <fstream>

#include "utilities/volumetric.h"
#include "utilities/utilities.h"

bool Vol_yaml_parameters_parser(const YAML::Node& parser, ImageBase::volume_parameters_structure& vol_parameters){
    if(parser["ReconstructionVolume"]["Dimension"]){
        vol_parameters.dim_vxl.x = parser["ReconstructionVolume"]["Dimension"]["x"].as<int>();
        vol_parameters.dim_vxl.y = parser["ReconstructionVolume"]["Dimension"]["y"].as<int>();
        vol_parameters.dim_vxl.z = parser["ReconstructionVolume"]["Dimension"]["z"].as<int>();
    }
    else{
        std::cout << "[ERROR]: missing entry of Dimension" << std::endl;
        return true;
    }

    if(parser["ReconstructionVolume"]["VoxelSize"]){
        vol_parameters.v_size_mm.x = parser["ReconstructionVolume"]["VoxelSize"]["x"].as<float>();
        vol_parameters.v_size_mm.y = parser["ReconstructionVolume"]["VoxelSize"]["y"].as<float>();
        vol_parameters.v_size_mm.z = parser["ReconstructionVolume"]["VoxelSize"]["z"].as<float>();
    }
    else{
        std::cout << "[ERROR]: missing entry of VoxelSize" << std::endl;
        return true;
    }

    if(parser["ReconstructionVolume"]["Origin"]){
        vol_parameters.origin_mm.x = parser["ReconstructionVolume"]["Origin"]["x"].as<float>();
        vol_parameters.origin_mm.y = parser["ReconstructionVolume"]["Origin"]["y"].as<float>();
        vol_parameters.origin_mm.z = parser["ReconstructionVolume"]["Origin"]["z"].as<float>();
    }
    else{
        std::cout << "[ERROR]: missing entry of Origin" << std::endl;
        return true;
    }

    return false;
}

ImageBase::Volume_Container::Volume_Container(std::string Input_Parameter_Path){
    Input_Parameter_Path_internal = Input_Parameter_Path;

    VOL_ERROR = false;
    parameters_file_handle = YAML::LoadFile(Input_Parameter_Path);

    if(!parameters_file_handle["ReconstructionVolume"]){
        // looking for ReconstructionVolume entry in parameter file
        std::cout << "[ERROR]: please check if ReconstructionVolume entry exists in file: " << Input_Parameter_Path << std::endl;
        VOL_ERROR = true;
        return;
    }

    //parse basic parameters:
    if(Vol_yaml_parameters_parser(parameters_file_handle, vol_parameters)){
        VOL_ERROR = true;
        return;
    }

    //Output Path:
    if(parameters_file_handle["ReconstructionVolume"]["OutputPath"]){
        Output_Path = parameters_file_handle["ReconstructionVolume"]["OutputPath"].as<std::string>();
    }
    else{
        std::cout << "[ERROR]: missing entry of OutputPath" << std::endl;
        VOL_ERROR = true;
        return;
    }

    //parse transform info:
    if(parameters_file_handle["ReconstructionVolume"]["Tramsforms"]){
        ReferenceTracker_list_path = parameters_file_handle["ReconstructionVolume"]["Tramsforms"]["TrackerListsPath"].as<std::string>();
        ReferenceTrackingNumber = parameters_file_handle["ReconstructionVolume"]["Tramsforms"]["TrackingDataNumber"].as<int>();
    }
    else{
        std::cout << "[ERROR]: missing entry of Tramsforms" << std::endl;
        VOL_ERROR = true;
        return;
    }
} 

bool ImageBase::Volume_Container::setup_VolumeToReconstruction(){
    if(ReferenceTrackingNumber > 1){
        //Initialize volumetric data:
        //TODO

        //Import Reference tracking data:
        ReferenceToWorldMatrices.resize( 
            ReferenceTrackingNumber, 
            std::vector<float>(16, 0.0f)
        );
        YAML::Node TrackingData_handle = YAML::LoadFile(ReferenceTracker_list_path);
        if(TrackingData_handle.size() != ReferenceTrackingNumber){
            std::cout << "[ERROR]: Tracking data is insufficient." << std::endl;
            VOL_ERROR = true;
            return VOL_ERROR;
        }
        else{
            for(int it_TrackingData = 0; it_TrackingData < TrackingData_handle.size(); ++it_TrackingData){
                std::vector<std::vector<float>> GotMatrix = TrackingData_handle["Frame#" + std::to_string(it_TrackingData)]["ReferenceToWorldMatrix"].as<std::vector<std::vector<float>>>();
                ReferenceToWorldMatrices[it_TrackingData] = Utilities::MatrixVector2D_Deserializer(GotMatrix);
            }
            std::cout << "[STATE]: " << ReferenceToWorldMatrices.size() << " ReferenceToWorld tracking data imported. " << std::endl;
        }

        Data_Reconstructed = false;
        Data_Valid = false;
        return VOL_ERROR;
    }
    else{
        std::cout << "[ERROR]: Trackingdata for Reference is invalid, please use ImportReconstruction" << std::endl;
        VOL_ERROR = true;
        return VOL_ERROR;
    }
}

void ImageBase::Volume_Container::SetVolumeDataValid(){
    Data_Valid = true;
    // parameters_file_handle = YAML::LoadFile(Input_Parameter_Path_internal);
    parameters_file_handle["ReconstructionVolume"]["Dimension"]["x"] = vol_parameters.dim_vxl.x;
    parameters_file_handle["ReconstructionVolume"]["Dimension"]["y"] = vol_parameters.dim_vxl.y;
    parameters_file_handle["ReconstructionVolume"]["Dimension"]["z"] = vol_parameters.dim_vxl.z;

    parameters_file_handle["ReconstructionVolume"]["Origin"]["x"] = vol_parameters.origin_mm.x;
    parameters_file_handle["ReconstructionVolume"]["Origin"]["y"] = vol_parameters.origin_mm.y;
    parameters_file_handle["ReconstructionVolume"]["Origin"]["z"] = vol_parameters.origin_mm.z;
}

void ImageBase::Volume_Container::SetVolumeDataReconstructed(){
    Data_Reconstructed = true;
    // parameters_file_handle = YAML::LoadFile(Input_Parameter_Path_internal);
    parameters_file_handle["ReconstructionVolume"]["IsReconstructed"] = true;
}

bool ImageBase::Volume_Container::IsVolumeDataValid(){
    return Data_Valid;
}

bool ImageBase::Volume_Container::IsVolumeDataReconstructed(){
    return Data_Reconstructed;
}

ImageBase::Volume_Container::~Volume_Container(){
    std::ofstream fout(Input_Parameter_Path_internal); 
    fout << parameters_file_handle;
}