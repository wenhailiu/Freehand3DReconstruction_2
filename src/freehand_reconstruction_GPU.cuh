#ifndef CUDA_FREERECON
#define CUDA_FREERECON

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "utilities/img_base.h"
#include "utilities/volumetric.h"
#include "utilities/ultrasound_frames.h"

void GPU_Setups(const ImageBase::us_parameters_structure US_Params, const ImageBase::volume_parameters_structure Vol_Params, const int NumFrames, const float* TotalMatrices, float* Recon_Volume, float* Weighting_Volume, uint8_t *VolumeTosave, uint8_t* US_Frames);

#endif