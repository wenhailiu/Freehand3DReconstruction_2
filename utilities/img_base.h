#ifndef IMAGE_BASE
#define IMAGE_BASE

namespace ImageBase{
    struct Dim_vxl{
        int x;
        int y;
        int z;
    };

    struct Dim_pxl{
        int x;
        int y;
    };

    struct VoxelSize_mm{
        float x;
        float y;
        float z;
    };

    struct PixelSize_mm{
        float x;
        float y;
        float z;
    };

    struct SpatialOrigin{
        float x;
        float y;
        float z;
    };
}

#endif