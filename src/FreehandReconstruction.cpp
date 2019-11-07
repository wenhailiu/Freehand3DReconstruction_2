#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>

#include "FreehandReconstruction.h"
#include "utilities/utilities.h"

#include "freehand_reconstruction_GPU.cuh"

uint8_t sliceMinus(uint8_t x, uint8_t y){
    return (uint8_t)std::abs(int(x) - int(y));
}

FreehandReconstruction::FreehandReconstruction(std::string Initialization_parameter_path): 
ImageBase::Volume_Container(Initialization_parameter_path), 
ImageBase::Ultrasound_Container(Initialization_parameter_path)
{   
    Reconstructor_ParamHandle = YAML::LoadFile(Initialization_parameter_path);

    RECON_ERROR = false;
    RawPtrSet = false;
    ImageToVolume_hostPtr = NULL;

    if(!Reconstructor_ParamHandle["ScanInfo"]){
        std::cout << "[ERROR]: Entry: ScanInfo was not found. " << std::endl;
        RECON_ERROR = true;
        return;
    }

    //Interpret as 2D freehand Scan:
    //No need to extract valid slices from the scan.
    if(Reconstructor_ParamHandle["ScanInfo"]["Type"].as<std::string>() == "2D"){

        //Reset angle begin and end points:
        angleBegin = 0.0f;
        angleEnd = 0.0f;
        Motor_Radius = 0.0f;
        TrimThreshold = 0.0f;

        if(ImageBase::Ultrasound_Container::setup_OnetimeReconstruction()){
            std::cout << "[ERROR]: Error occurred while constructing Ultrasound container. " << std::endl;
            RECON_ERROR = true;
            return;
        }

        if(ImageBase::Volume_Container::setup_VolumeToReconstruction()){
            std::cout << "[ERROR]: Error occurred while constructing Volumetrix container. " << std::endl;
            RECON_ERROR = true;
            return;
        }

        if(ReferenceToWorldMatrices.size() != ProbeToWorldMatrices.size()){
            std::cout << "[ERROR]: Error occurred: Reference and Probe TrackingData not match. " << std::endl;
            RECON_ERROR = true;
            return;
        }
        else{
            TiltedImageMatrices.resize(ReferenceToWorldMatrices.size(), std::vector<float>(16, 0.0f));
            Utilities::MatrixS4X4Identity(TiltedImageMatrices);
        }
    }

    //Interpret as 3D onetime Scan:
    //Need to extract valid slices from the scan. 
    else if(Reconstructor_ParamHandle["ScanInfo"]["Type"].as<std::string>() == "3D"){

        //Set angle range:
        angleBegin = Reconstructor_ParamHandle["ScanInfo"]["AngleBegin"].as<float>();
        angleEnd = Reconstructor_ParamHandle["ScanInfo"]["AngleEnd"].as<float>();
        Motor_Radius = Reconstructor_ParamHandle["ScanInfo"]["MotorRadius"].as<float>();
        TrimThreshold = Reconstructor_ParamHandle["ScanInfo"]["TrimThreshold"].as<float>();

        if(ImageBase::Ultrasound_Container::setup_OnetimeReconstruction()){
            std::cout << "[ERROR]: Error occurred while constructing Ultrasound container. " << std::endl;
            RECON_ERROR = true;
            return;
        }

        if(ImageBase::Volume_Container::setup_VolumeToReconstruction()){
            std::cout << "[ERROR]: Error occurred while constructing Volumetrix container. " << std::endl;
            RECON_ERROR = true;
            return;
        }

        if(ReferenceToWorldMatrices.size() != ProbeToWorldMatrices.size()){
            std::cout << "[ERROR]: Error occurred: Reference and Probe TrackingData not match. " << std::endl;
            RECON_ERROR = true;
            return;
        }
        else{

            //Step 1: Calculate inter-slice standard deviation. 
            std::vector<float> std_slices;
            {
                std::vector<uint8_t> CurrentSlice(us_parameters.dim_pxl.x * us_parameters.dim_pxl.y, 0);
                for(int slice_it = 0; slice_it < (NumFrames - 1); ++slice_it){
                    std::transform( 
                        Frames_Raw_Data.begin() + (slice_it + 1) * us_parameters.dim_pxl.x * us_parameters.dim_pxl.y, Frames_Raw_Data.begin() + (slice_it + 1) * us_parameters.dim_pxl.x * us_parameters.dim_pxl.y + us_parameters.dim_pxl.x * us_parameters.dim_pxl.y, 
                        Frames_Raw_Data.begin() + slice_it * us_parameters.dim_pxl.x * us_parameters.dim_pxl.y, CurrentSlice.begin(), 
                        sliceMinus
                    );

                    float mean_slice = std::accumulate(CurrentSlice.begin(), CurrentSlice.end(), 0) / (float)CurrentSlice.size();
                    float slice_var = 0.0f;
                    for(int i = 0; i < CurrentSlice.size(); ++i){
                        slice_var += static_cast<float>(((float)CurrentSlice[i] - mean_slice) * ((float)CurrentSlice[i] - mean_slice));
                    }
                    slice_var /= (CurrentSlice.size());

                    std_slices.push_back(std::sqrt(slice_var));
                }
            }
            
            //Step 2: Determine valid slices range. 
            int ValidSlice_BegIndex = 0, ValidSlice_EndIndex = 0, TrimmedSize = 0;
            {
                ValidSlice_BegIndex = std::lower_bound(std_slices.begin(), std_slices.end(), TrimThreshold) - std_slices.begin();
                
                std::reverse(std_slices.begin(), std_slices.end());
                ValidSlice_EndIndex = std_slices.size() - (std::lower_bound(std_slices.begin(), std_slices.end(), TrimThreshold) - std_slices.begin());
                
                TrimmedSize = ValidSlice_EndIndex - ValidSlice_BegIndex + 1;

                std::cout << "[STATE]: Valid slices range: " << ValidSlice_BegIndex << " to " << ValidSlice_EndIndex << std::endl;
            }
            
            //Step 3: Trim: Frames_Raw_Data, ReferenceToTrackerMatrices, ProbeToTrackerMatrices.
            {
                ProbeToWorldMatrices.erase(ProbeToWorldMatrices.begin() + ValidSlice_EndIndex + 1, ProbeToWorldMatrices.end());
                ProbeToWorldMatrices.erase(ProbeToWorldMatrices.begin(), ProbeToWorldMatrices.begin() + ValidSlice_BegIndex);
                std::cout << "[STATE]: ProbeToWorldMatrices trimmed, new size: " << ProbeToWorldMatrices.size() << std::endl;

                ReferenceToWorldMatrices.erase(ReferenceToWorldMatrices.begin() + ValidSlice_EndIndex + 1, ReferenceToWorldMatrices.end());
                ReferenceToWorldMatrices.erase(ReferenceToWorldMatrices.begin(), ReferenceToWorldMatrices.begin() + ValidSlice_BegIndex);
                std::cout << "[STATE]: ReferenceToWorldMatrices trimmed, new size: " << ReferenceToWorldMatrices.size() << std::endl;

                Frames_Raw_Data.erase( 
                    Frames_Raw_Data.begin() + (ValidSlice_EndIndex + 1) * us_parameters.dim_pxl.x * us_parameters.dim_pxl.y, 
                    Frames_Raw_Data.end()
                );
                Frames_Raw_Data.erase( 
                    Frames_Raw_Data.begin(), 
                    Frames_Raw_Data.begin() + (ValidSlice_BegIndex) * us_parameters.dim_pxl.x * us_parameters.dim_pxl.y
                );
                std::cout << "[STATE]: Frames_Raw_Data trimmed, new size: " << Frames_Raw_Data.size() / us_parameters.dim_pxl.x / us_parameters.dim_pxl.y << std::endl;

                NumFrames = TrimmedSize;
                ProbeTrackingNumber = TrimmedSize;
                ReferenceTrackingNumber = TrimmedSize;
            }

            //Step 4: Deduce the relative transform matrices:
            {
                TiltedImageMatrices.resize(TrimmedSize, std::vector<float>(16, 0.0f));
                Utilities::MatrixS4X4Identity(TiltedImageMatrices);

                std::vector<float> AnglesVector(TrimmedSize, 0.0f);
                for(int i = 0; i < AnglesVector.size(); ++i){
                    AnglesVector[i] = angleBegin + i * ((angleEnd - angleBegin) / (TrimmedSize - 1));
                }

                //Initialize Matrix:
                for(int i = 0; i < AnglesVector.size(); ++i){
                    float cosAnagle = std::cos(AnglesVector[i] * M_PI / 180.0);
                    float sinAnagle = std::sin(AnglesVector[i] * M_PI / 180.0);
                    //Translation Y:
                    TiltedImageMatrices[i][7] = Motor_Radius * cosAnagle - Motor_Radius;
                    TiltedImageMatrices[i][11] = Motor_Radius * sinAnagle;

                    //Rotation matrix:
                    TiltedImageMatrices[i][5] = cosAnagle;
                    TiltedImageMatrices[i][6] = -sinAnagle;
                    TiltedImageMatrices[i][9] = sinAnagle;
                    TiltedImageMatrices[i][10] = cosAnagle;

                    Utilities::MatrixS4X4Invert(TiltedImageMatrices[i].data(), TiltedImageMatrices[i].data());
                }
                // std::cout << std::endl;
            }
        }
    }
    else{
        std::cout << "[ERROR]: ScanInfo Type should be: 2D or 3D. " << std::endl;
        RECON_ERROR = true;
        return;
    }
}

FreehandReconstruction::FreehandReconstruction(std::string Initialization_parameter_path, const float StartAngle, const float EndAngle): 
ImageBase::Volume_Container(Initialization_parameter_path), 
ImageBase::Ultrasound_Container(Initialization_parameter_path)
{
    RECON_ERROR = false;
    RawPtrSet = false;
    ImageToVolume_hostPtr = NULL;

    if(ReferenceToWorldMatrices.size() != ProbeToWorldMatrices.size()){
        std::cout << "[ERROR]: Error occurred: Reference and Probe TrackingData not match. " << std::endl;
        RECON_ERROR = true;
        return;
    }
    else{
        TiltedImageMatrices.resize(ReferenceToWorldMatrices.size(), std::vector<float>(16, 0.0f));
        Utilities::MatrixS4X4Identity(TiltedImageMatrices);

        
    }
}

bool FreehandReconstruction::ExtractImageToVolumeMatrices(int smoothPatch){
    //Check if the ImageToProbeMatrix is normalized. 
    //If not, then normalized with current Ultrasound pixel size:
    std::cout << "[STATE]: Image To Probe Matrix: " << std::endl;
    Utilities::MatrixPrint(ImageBase::Ultrasound_Container::ImageToProbeMatrix);
    { //Normalization Check: 
        float norm_x = std::sqrt(
            ImageToProbeMatrix[0] * ImageToProbeMatrix[0] + 
            ImageToProbeMatrix[1] * ImageToProbeMatrix[1] + 
            ImageToProbeMatrix[2] * ImageToProbeMatrix[2]
        );

        float norm_y = std::sqrt(
            ImageToProbeMatrix[4] * ImageToProbeMatrix[4] + 
            ImageToProbeMatrix[5] * ImageToProbeMatrix[5] + 
            ImageToProbeMatrix[6] * ImageToProbeMatrix[6]
        );

        float norm_z = std::sqrt(
            ImageToProbeMatrix[8] * ImageToProbeMatrix[8] + 
            ImageToProbeMatrix[9] * ImageToProbeMatrix[9] + 
            ImageToProbeMatrix[10] * ImageToProbeMatrix[10]
        );

        std::cout << "[STATE]: ImageToProbeMatrix has NormX: " << norm_x << ", NormY: " << norm_y << ", NormZ: " << norm_z << std::endl;

        std::cout << "[STATE]: Adding current pixel size X: " << us_parameters.p_size_mm.x << ", Y: " << us_parameters.p_size_mm.y << ", Z: " << us_parameters.p_size_mm.z << std::endl;
        ImageToProbeMatrix[0] /= norm_x;
        // ImageToProbeMatrix[0] *= us_parameters.p_size_mm.x;
        ImageToProbeMatrix[1] /= norm_x;
        // ImageToProbeMatrix[1] *= us_parameters.p_size_mm.x;
        ImageToProbeMatrix[2] /= norm_x;
        // ImageToProbeMatrix[2] *= us_parameters.p_size_mm.x;
        ImageToProbeMatrix[3] /= norm_x;
        ImageToProbeMatrix[3] *= us_parameters.p_size_mm.x;

        ImageToProbeMatrix[4] /= norm_y;
        // ImageToProbeMatrix[4] *= us_parameters.p_size_mm.y;
        ImageToProbeMatrix[5] /= norm_y;
        // ImageToProbeMatrix[5] *= us_parameters.p_size_mm.y;
        ImageToProbeMatrix[6] /= norm_y;
        // ImageToProbeMatrix[6] *= us_parameters.p_size_mm.y;

        ImageToProbeMatrix[8] /= norm_z;
        // ImageToProbeMatrix[8] *= us_parameters.p_size_mm.z;
        ImageToProbeMatrix[9] /= norm_z;
        // ImageToProbeMatrix[9] *= us_parameters.p_size_mm.z;
        ImageToProbeMatrix[10] /= norm_z;
        // ImageToProbeMatrix[10] *= us_parameters.p_size_mm.z;
        
        //Print new matrix:
        std::cout << "[STATE]: Image To Probe Matrix with current PixelSize: " << std::endl;
        Utilities::MatrixPrint(ImageBase::Ultrasound_Container::ImageToProbeMatrix);
    }

    //Go through matrices, calculate ImageToVolumeMatrix:
    //Make correction to Volume Dimension:
    {
        if(ReferenceToWorldMatrices.size() == ProbeToWorldMatrices.size()){
            MatchedMatrixNumber = ReferenceToWorldMatrices.size();
            ImageToVolumeMatrices.resize(MatchedMatrixNumber, std::vector<float>(16, 0.0f));
            
            float point1[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            float point1_out[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            float point2[4] = {(float)us_parameters.dim_pxl.x, 0.0f, 0.0f, 1.0f};
            float point2_out[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            float point3[4] = {0.0f, (float)us_parameters.dim_pxl.y, 0.0f, 1.0f};
            float point3_out[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            float point4[4] = {(float)us_parameters.dim_pxl.x, (float)us_parameters.dim_pxl.y, 0.0f, 1.0f};
            float point4_out[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            std::vector<float> vPoint_X, vPoint_Y, vPoint_Z;

            for(int Matrix_it = 0; Matrix_it < MatchedMatrixNumber; ++Matrix_it){
                float TiledImageToProbe[16] = {0.0f};
                Utilities::MatrixMatrixS4X4Multiply(
                    ImageToProbeMatrix.data(), 
                    TiltedImageMatrices[Matrix_it].data(), 
                    TiledImageToProbe
                );

                TiledImageToProbe[0] *= us_parameters.p_size_mm.x;
                TiledImageToProbe[1] *= us_parameters.p_size_mm.x;
                TiledImageToProbe[2] *= us_parameters.p_size_mm.x;

                TiledImageToProbe[4] *= us_parameters.p_size_mm.y;
                TiledImageToProbe[5] *= us_parameters.p_size_mm.y;
                TiledImageToProbe[6] *= us_parameters.p_size_mm.y;

                TiledImageToProbe[8] *= us_parameters.p_size_mm.z;
                TiledImageToProbe[9] *= us_parameters.p_size_mm.z;
                TiledImageToProbe[10] *= us_parameters.p_size_mm.z;

                float ImageToWorld[16] = {0.0f};
                Utilities::MatrixMatrixS4X4Multiply( 
                    ProbeToWorldMatrices[Matrix_it].data(), 
                    TiledImageToProbe, 
                    ImageToWorld);
                
                float WorldToReference[16] = {0.0f};
                Utilities::MatrixS4X4Invert(ReferenceToWorldMatrices[Matrix_it].data(), WorldToReference);

                std::vector<float> ImageToVolume(16, 0.0f);
                Utilities::MatrixMatrixS4X4Multiply( 
                    WorldToReference, 
                    ImageToWorld, 
                    ImageToVolume.data()
                );

                ImageToVolumeMatrices[Matrix_it] = ImageToVolume;
                Utilities::MatrixPointS4X4Multiply(ImageToVolume.data(), point1, point1_out);
                vPoint_X.push_back(point1_out[0]);
                vPoint_Y.push_back(point1_out[1]);
                vPoint_Z.push_back(point1_out[2]);

                ImageToVolumeMatrices[Matrix_it] = ImageToVolume;
                Utilities::MatrixPointS4X4Multiply(ImageToVolume.data(), point2, point2_out);
                vPoint_X.push_back(point2_out[0]);
                vPoint_Y.push_back(point2_out[1]);
                vPoint_Z.push_back(point2_out[2]);

                ImageToVolumeMatrices[Matrix_it] = ImageToVolume;
                Utilities::MatrixPointS4X4Multiply(ImageToVolume.data(), point3, point3_out);
                vPoint_X.push_back(point3_out[0]);
                vPoint_Y.push_back(point3_out[1]);
                vPoint_Z.push_back(point3_out[2]);

                ImageToVolumeMatrices[Matrix_it] = ImageToVolume;
                Utilities::MatrixPointS4X4Multiply(ImageToVolume.data(), point4, point4_out);
                vPoint_X.push_back(point4_out[0]);
                vPoint_Y.push_back(point4_out[1]);
                vPoint_Z.push_back(point4_out[2]);
            }

            float max_position_X_Volume = *std::max_element(vPoint_X.begin(), vPoint_X.end());
            float min_position_X_Volume = *std::min_element(vPoint_X.begin(), vPoint_X.end());

            float max_position_Y_Volume = *std::max_element(vPoint_Y.begin(), vPoint_Y.end());
            float min_position_Y_Volume = *std::min_element(vPoint_Y.begin(), vPoint_Y.end());

            float max_position_Z_Volume = *std::max_element(vPoint_Z.begin(), vPoint_Z.end());
            float min_position_Z_Volume = *std::min_element(vPoint_Z.begin(), vPoint_Z.end());

            std::cout << "[STATE]: " << "Calculated Volume Origin X: " << min_position_X_Volume << "mm to " << max_position_X_Volume << 
            ", Y: " << min_position_Y_Volume << "mm to " << max_position_Y_Volume << 
            ", Z: " << min_position_Z_Volume << "mm to " << max_position_Z_Volume << std::endl;

            vol_parameters.origin_mm.x = min_position_X_Volume;
            vol_parameters.origin_mm.y = min_position_Y_Volume;
            vol_parameters.origin_mm.z = min_position_Z_Volume;

            vol_parameters.dim_vxl.x = int(std::ceil((max_position_X_Volume - min_position_X_Volume) / vol_parameters.v_size_mm.x) + 1);
            vol_parameters.dim_vxl.y = int(std::ceil((max_position_Y_Volume - min_position_Y_Volume) / vol_parameters.v_size_mm.y) + 1);
            vol_parameters.dim_vxl.z = int(std::ceil((max_position_Z_Volume - min_position_Z_Volume) / vol_parameters.v_size_mm.z) + 1);

            std::cout << "[STATE]: " << "Calculated Volume Dimension X: " << vol_parameters.dim_vxl.x << ", Y: " << vol_parameters.dim_vxl.y << ", Z: " << vol_parameters.dim_vxl.z << std::endl;

            ImageBase::Volume_Container::Volume_data.resize(vol_parameters.dim_vxl.x * vol_parameters.dim_vxl.y * vol_parameters.dim_vxl.z, 0.0f);
            ImageBase::Volume_Container::Weights_data.resize(vol_parameters.dim_vxl.x * vol_parameters.dim_vxl.y * vol_parameters.dim_vxl.z, 0.0f);

            SetVolumeDataValid();
        }
        else{
            std::cout << "[ERROR]: Number of Reference and Probe Tracking matrices Not Match. " << std::endl;
            RECON_ERROR = true;
            return RECON_ERROR;
        }
    }
    
    //TODO: Smoothen translations, rotations...
}

bool FreehandReconstruction::MallocHostSpace(){
    //Allocate Space for transform matrix: 
    ImageToVolume_hostPtr = new float[16 * NumFrames];

    for(int Matrix_it = 0; Matrix_it < NumFrames; ++Matrix_it){
        for(int ele_it = 0; ele_it < 16; ++ele_it){
            ImageToVolume_hostPtr[Matrix_it * 16 + ele_it] = ImageToVolumeMatrices[Matrix_it][ele_it];
        }
    }
    
    RawPtrSet = true;

}

bool FreehandReconstruction::LaunchGPU_Reconstruction(){
    
    uint8_t *UltrasoundFrames_hostPtr = ImageBase::Ultrasound_Container::Frames_Raw_Data.data();
    float *Volume_hostPtr = ImageBase::Volume_Container::Volume_data.data();
    float *weightings_hostPtr = ImageBase::Volume_Container::Weights_data.data();

    if(RawPtrSet && IsVolumeDataValid()){
        GPU_Setups(us_parameters, vol_parameters, NumFrames, ImageToVolume_hostPtr, Volume_hostPtr, weightings_hostPtr, UltrasoundFrames_hostPtr);
        SetVolumeDataReconstructed();
    }
    else{
        std::cout << "[ERROR]: Allocate space for raw pointers first. " << std::endl;
        RECON_ERROR = true;
        return RECON_ERROR;
    }   
}

bool FreehandReconstruction::SaveReconstructedVolume(){
    if(IsVolumeDataReconstructed()){
        Utilities::writeToBin( 
            ImageBase::Volume_Container::Volume_data.data(), 
            vol_parameters.dim_vxl.x * vol_parameters.dim_vxl.y * vol_parameters.dim_vxl.z, 
            ImageBase::Volume_Container::Output_Path
        );

        std::cout << "[STATE]: Reconstruction finished! Result volume saved: " << ImageBase::Volume_Container::Output_Path << std::endl;
    }
    else{
        std::cout << "[ERROR]: Reconstruction Failed. " << std::endl;
        RECON_ERROR = true;
        return RECON_ERROR;
    }
}

FreehandReconstruction::~FreehandReconstruction(){
    if(RawPtrSet && ImageToVolume_hostPtr != NULL){
        delete[] ImageToVolume_hostPtr;
        std::cout << "[STATE]: Clear all host raw memory pointers" << std::endl;
    }
    else{
        std::cout << "[ERROR]: Pointers are NULL, and were released somewhere else. " << std::endl;
    }
}