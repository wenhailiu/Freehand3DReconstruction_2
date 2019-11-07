#include <iostream>
#include <string>

#include "ultrasound_frames.h"
#include "utilities/utilities.h"

bool US_yaml_parameters_parser(const YAML::Node& parser, ImageBase::us_parameters_structure& us_parameters){
    if(parser["Ultrasound"]["Dimension"]){
        us_parameters.dim_pxl.x = parser["Ultrasound"]["Dimension"]["x"].as<int>();
        us_parameters.dim_pxl.y = parser["Ultrasound"]["Dimension"]["y"].as<int>();
    }
    else{
        std::cout << "[ERROR]: missing entry of Dimension" << std::endl;
        return true;
    }

    if(parser["Ultrasound"]["PixelSize"]){
        us_parameters.p_size_mm.x = parser["Ultrasound"]["PixelSize"]["x"].as<float>();
        us_parameters.p_size_mm.y = parser["Ultrasound"]["PixelSize"]["y"].as<float>();
        us_parameters.p_size_mm.z = parser["Ultrasound"]["PixelSize"]["z"].as<float>();
    }
    else{
        std::cout << "[ERROR]: missing entry of PixelSize" << std::endl;
        return true;
    }

    if(parser["Ultrasound"]["Origin"]){
        us_parameters.origin_mm.x = parser["Ultrasound"]["Origin"]["x"].as<float>();
        us_parameters.origin_mm.y = parser["Ultrasound"]["Origin"]["y"].as<float>();
        us_parameters.origin_mm.z = parser["Ultrasound"]["Origin"]["z"].as<float>();
    }
    else{
        std::cout << "[ERROR]: missing entry of Origin" << std::endl;
        return true;
    }

    return false;
}

ImageBase::Ultrasound_Container::Ultrasound_Container(std::string Input_Parameter_Path){
    US_ERROR = false;

    parameters_file_handle = YAML::LoadFile(Input_Parameter_Path);

    if(!parameters_file_handle["Ultrasound"]){
        // looking for Ultrasound entry in parameter file
        std::cout << "[ERROR]: please check if Ultrasound entry exists in file: " << Input_Parameter_Path << std::endl;
        US_ERROR = true;
        return;
    }

    //parse basic parameters:
    if(US_yaml_parameters_parser(parameters_file_handle, us_parameters)){
        US_ERROR = true;
        return;
    }

    //parse Frame info:
    if(parameters_file_handle["Ultrasound"]["FramesInfo"]){
        frame_data_path = parameters_file_handle["Ultrasound"]["FramesInfo"]["Path"].as<std::string>();
        NumFrames = parameters_file_handle["Ultrasound"]["FramesInfo"]["FrameNumber"].as<int>();
    }
    else{
        std::cout << "[ERROR]: missing entry of FramesInfo" << std::endl;
        US_ERROR = true;
        return;
    }

    //parse transform info:
    if(parameters_file_handle["Ultrasound"]["Tramsforms"]){
        ProbeTracker_list_path = parameters_file_handle["Ultrasound"]["Tramsforms"]["TrackerListsPath"].as<std::string>();
        ProbeTrackingNumber = parameters_file_handle["Ultrasound"]["Tramsforms"]["TrackingDataNumber"].as<int>();
        std::vector<std::vector<float>> ImageToProbeMatrix_2Dvector = parameters_file_handle["Ultrasound"]["Tramsforms"]["ImageToProbe"].as<std::vector<std::vector<float>>>();
        ImageToProbeMatrix = Utilities::MatrixVector2D_Deserializer(ImageToProbeMatrix_2Dvector);
        ProbeCalibrationDate = parameters_file_handle["Ultrasound"]["Tramsforms"]["ProbeCalibrationDate"].as<std::string>();

        if(ProbeTrackingNumber != NumFrames){
            std::cout << "[ERROR]: Tracking data does not match with FrameNumber." << std::endl;
            US_ERROR = true;
            return;
        }
    }
    else{
        std::cout << "[ERROR]: missing entry of Tramsforms" << std::endl;
        US_ERROR = true;
        return;
    }

}

bool ImageBase::Ultrasound_Container::setup_OnetimeReconstruction(){
    if(NumFrames > 1 && ProbeTrackingNumber > 1){
        //Initialize Frame raw data:
        Frames_Raw_Data.resize( 
            us_parameters.dim_pxl.x * us_parameters.dim_pxl.y * 
            NumFrames, 
            0);
        Utilities::readFromBin(Frames_Raw_Data.data(), us_parameters.dim_pxl.x * us_parameters.dim_pxl.y * NumFrames, frame_data_path);
        std::cout << "[STATE]: " << us_parameters.dim_pxl.x << " * " << us_parameters.dim_pxl.y << " * " << NumFrames << " bytes FrameData initialized. " << std::endl;

        //Import Probe tracking data:
        ProbeToWorldMatrices.resize( 
            ProbeTrackingNumber, 
            std::vector<float>(16, 0.0f)
        );
        YAML::Node TrackingData_handle = YAML::LoadFile(ProbeTracker_list_path);
        if(TrackingData_handle.size() != ProbeTrackingNumber){
            std::cout << "[ERROR]: Tracking data is insufficient." << std::endl;
            US_ERROR = true;
            return US_ERROR;
        }
        else{
            for(int it_TrackingData = 0; it_TrackingData < TrackingData_handle.size(); ++it_TrackingData){
                std::vector<std::vector<float>> GotMatrix = TrackingData_handle["Frame#" + std::to_string(it_TrackingData)]["ProbeToWorldMatrix"].as<std::vector<std::vector<float>>>();
                ProbeToWorldMatrices[it_TrackingData] = Utilities::MatrixVector2D_Deserializer(GotMatrix);
            }
            std::cout << "[STATE]: " << ProbeToWorldMatrices.size() << " ProbeToWorld tracking data imported. " << std::endl;
        }
        return US_ERROR;
    }
    else{
        std::cout << "[ERROR]: NumFrames is invalid for OnetimeReconstruction, please use ProgressiveReconstruction." << std::endl;
        US_ERROR = true;
        return US_ERROR;
    }
}