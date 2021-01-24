
/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#include <opencv2\opencv.hpp>
#include "moments.hpp"

template <typename T>
void projections(Volume<T> volume, cv::Mat& x_y, cv::Mat &y_z, cv::Mat& x_z, cv::Mat& xz_y, cv::Mat& zx_y, int order)
{
    int w = volume.w, h = volume.h, d = volume.d;
    int slice = w * h;

    if(order < 0 || order > 4) return;

    const T* py = volume.voxels;

    for(int z = 0; z < d; z++) {
        int* p_y_z = y_z.ptr<int>(z);
        int* p_x_z = x_z.ptr<int>(z);

        const T* p = py;
        for(int y = 0; y < h; y++) {
            int sumx = 0;

            int* p_x_y = x_y.ptr<int>(y);

            if(order == 4) {
                int* p_xz_y = &xz_y.ptr<int>(y)[z];
                int* p_zx_y = &zx_y.ptr<int>(y)[d-z-1];
                for(int x = 0; x < w; x++) {
                    p_x_y[x] += p[x];
                    p_x_z[x] += p[x];
                    p_xz_y[x] += p[x];
                    p_zx_y[x] += p[x];
                    sumx += p[x];
                }

                p_y_z[y] += sumx;
            }
            else if(order == 3) {
                int* p_xz_y = &xz_y.ptr<int>(y)[z];
                for(int x = 0; x < w; x++) {
                    p_x_y[x] += p[x];
                    p_x_z[x] += p[x];
                    p_xz_y[x] += p[x];
                    sumx += p[x];
                }

                p_y_z[y] += sumx;
            }
            else if(order == 2) {
                for(int x = 0; x < w; x++) {
                    p_x_y[x] += p[x];
                    p_x_z[x] += p[x];
                    sumx += p[x];
                }

                p_y_z[y] += sumx;                
            }
            else if(order == 1) {
                for(int x = 0; x < w; x++) {
                    p_x_y[x] += p[x];
                    sumx += p[x];
                }

                p_y_z[y] += sumx;                
            }

            p += w;
        }

        py += slice;
    }
}
 
template void projections<char>(Volume<char>, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, int);
template void projections<uchar>(Volume<uchar>, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, int);
template void projections<short>(Volume<short>, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, int);
template void projections<ushort>(Volume<ushort>, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, int);

