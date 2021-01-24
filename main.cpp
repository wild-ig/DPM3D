/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#include <iostream>
#include <vector>
#include <numeric>
#include <Windows.h>
#include <opencv2\opencv.hpp>
#include <opencv2\highgui.hpp>
#include <opencv2\core\core.hpp>

#include "moments.hpp"

#include "NRRD/nrrd_image.hxx"


using namespace std;
using namespace cv;

// Default dummy test "image"
uchar dummy[] = {
0,	2,	3,	1,
43,	54,	61,	1,
4,	52,	25,	4,
11,	33,	42,	7,
1,	0,	2,	1,
0,	4,	6,	2,
86,	108,	122,	2,
8,	104,	50,	8,
22,	66,	84,	14,
2,	0,	4,	2,
0,	4,	6,	2,
94,	118,	134,	2,
8,	114,	55,	8,
24,	72,	92,	15,
2,	0,	4,	2,
0,	4,	7,	2,
112,	141,	160,	2,
9,	136,	66,	9,
28,	86,	110,	18,
2,	0,	4,	2,
0,	2,	3,	1,
56,	70,	80,	1,
4,	68,	33,	4,
14,	43,	55,	9,
1,	0,	2,	1,
0,	1,	2,	0,
43,	53,	61,	0,
3,	52,	25,	3,
10,	33,	42,	6,
0,	0,	1,	0 };

// Order of moments to generate - can be overirdden by cmd line
int order = 3;

// Display an image
void ShowOff(const Mat& image, char *name) {

    double maxValue, minValue;
    cv::namedWindow(name);
    Mat convert;
    minMaxIdx(image, &minValue, &maxValue);
    image.convertTo(convert, CV_16UC1, 65536.0/(maxValue - minValue), -minValue/65536.0);
    cv::imshow(name, convert);

}

void SaveOff(const Mat& image, char *name) {
    int width = image.size().width;
    int height = image.size().height;
    NRRD::Image<ushort> volume;
    
    volume.set(width, height, 256);

    for(int k = 0; k < 256; k++) {
        for(int j = 0; j < height; j++) {
            const int *p = image.ptr<int>(j);
            for(int i = 0; i < width; i++) {

                if(p[i]/256 > k)
                    volume.pixel(i, j, k) = static_cast<ushort>(p[i]);
            }
        }
    }

    volume.save(name);
}

// Print the moments
void DumpMoments(Moments3D m0, Moments3D m1, Moments3D m2, Moments3D m3) {
    cout << "m000: " << m0.m000 << " m000: " << m1.m000 << " m000: " << m2.m000 << " m000: " << m3.m000 << endl;
    cout << "m100: " << m0.m100 << " m100: " << m1.m100 << " m100: " << m2.m100 << " m100: " << m3.m100 << endl;
    cout << "m010: " << m0.m010 << " m010: " << m1.m010 << " m010: " << m2.m010 << " m010: " << m3.m010 << endl;
    cout << "m001: " << m0.m001 << " m001: " << m1.m001 << " m001: " << m2.m001 << " m001: " << m3.m001 << endl;
    cout << "m200: " << m0.m200 << " m200: " << m1.m200 << " m200: " << m2.m200 << " m200: " << m3.m200 << endl;
    cout << "m020: " << m0.m020 << " m020: " << m1.m020 << " m020: " << m2.m020 << " m020: " << m3.m020 << endl;
    cout << "m002: " << m0.m002 << " m002: " << m1.m002 << " m002: " << m2.m002 << " m002: " << m3.m002 << endl;
    cout << "m110: " << m0.m110 << " m110: " << m1.m110 << " m110: " << m2.m110 << " m110: " << m3.m110 << endl;
    cout << "m101: " << m0.m101 << " m101: " << m1.m101 << " m101: " << m2.m101 << " m101: " << m3.m101 << endl;
    cout << "m011: " << m0.m011 << " m011: " << m1.m011 << " m011: " << m2.m011 << " m011: " << m3.m011 << endl;
    cout << "m300: " << m0.m300 << " m300: " << m1.m300 << " m300: " << m2.m300 << " m300: " << m3.m300 << endl;
    cout << "m030: " << m0.m030 << " m030: " << m1.m030 << " m030: " << m2.m030 << " m030: " << m3.m030 << endl;
    cout << "m003: " << m0.m003 << " m003: " << m1.m003 << " m003: " << m2.m003 << " m003: " << m3.m003 << endl;
    cout << "m210: " << m0.m210 << " m210: " << m1.m210 << " m210: " << m2.m210 << " m210: " << m3.m210 << endl;
    cout << "m201: " << m0.m201 << " m201: " << m1.m201 << " m201: " << m2.m201 << " m201: " << m3.m201 << endl;
    cout << "m120: " << m0.m120 << " m120: " << m1.m120 << " m120: " << m2.m120 << " m120: " << m3.m120 << endl;
    cout << "m021: " << m0.m021 << " m021: " << m1.m021 << " m021: " << m2.m021 << " m021: " << m3.m021 << endl;
    cout << "m102: " << m0.m102 << " m102: " << m1.m102 << " m102: " << m2.m102 << " m102: " << m3.m102 << endl;
    cout << "m012: " << m0.m012 << " m012: " << m1.m012 << " m012: " << m2.m012 << " m012: " << m3.m012 << endl;
    cout << "m111: " << m0.m111 << " m111: " << m1.m111 << " m111: " << m2.m111 << " m111: " << m3.m111 << endl;

    if(order > 3) {
        cout << "m400: " << m0.m400 << " m400: " << m1.m400 << " m400: " << m2.m400 << " m400: " << m3.m400 << endl;
        cout << "m040: " << m0.m040 << " m040: " << m1.m040 << " m040: " << m2.m040 << " m040: " << m3.m040 << endl;
        cout << "m004: " << m0.m004 << " m004: " << m1.m004 << " m004: " << m2.m004 << " m004: " << m3.m004 << endl;
        cout << "m310: " << m0.m310 << " m310: " << m1.m310 << " m310: " << m2.m310 << " m310: " << m3.m310 << endl;
        cout << "m301: " << m0.m301 << " m301: " << m1.m301 << " m301: " << m2.m301 << " m301: " << m3.m301 << endl;
        cout << "m130: " << m0.m130 << " m130: " << m1.m130 << " m130: " << m2.m130 << " m130: " << m3.m130 << endl;
        cout << "m031: " << m0.m031 << " m031: " << m1.m031 << " m031: " << m2.m031 << " m031: " << m3.m031 << endl;
        cout << "m103: " << m0.m103 << " m103: " << m1.m103 << " m103: " << m2.m103 << " m103: " << m3.m103 << endl;
        cout << "m013: " << m0.m013 << " m013: " << m1.m013 << " m013: " << m2.m013 << " m013: " << m3.m013 << endl;
        cout << "m211: " << m0.m211 << " m211: " << m1.m211 << " m211: " << m2.m211 << " m211: " << m3.m211 << endl;
        cout << "m121: " << m0.m121 << " m121: " << m1.m121 << " m121: " << m2.m121 << " m121: " << m3.m121 << endl;
        cout << "m112: " << m0.m112 << " m112: " << m1.m112 << " m112: " << m2.m112 << " m112: " << m3.m112 << endl;
        cout << "m220: " << m0.m220 << " m220: " << m1.m220 << " m220: " << m2.m220 << " m220: " << m3.m220 << endl;
        cout << "m202: " << m0.m202 << " m202: " << m1.m202 << " m202: " << m2.m202 << " m202: " << m3.m202 << endl;
        cout << "m022: " << m0.m022 << " m022: " << m1.m022 << " m022: " << m2.m022 << " m022: " << m3.m022 << endl;
    }

    cout << "x: " << m1.m100 / m1.m000 << endl;
    cout << "y: " << m1.m010 / m1.m000 << endl;
    cout << "z: " << m1.m001 / m1.m000 << endl;
}

template <typename T>
void ProcessImage(NRRD::Image<T> &img, bool show = true) {
    Volume<T> volume;
    
    // Supports N-dimensional images
    int dim = img.dimension();
    std::cout << "Number of dimensions: " << dim << std::endl;

    std::cout << "Spacing:";
    for (int i=0; i < dim;i++) {
        std::cout << " " << img.spacing(i);
    }
    std::cout << std::endl;

    std::cout << "Sizes:";
    for (int i=0; i < dim;i++) {
        std::cout << " " << img.size(i);
    }
    std::cout << std::endl;

    if(dim != 3) return;

    // Create our volume representation
    volume.voxels = img;
    volume.w = img.size(0);
    volume.h = img.size(1);
    volume.d = img.size(2);

    // maximum height/width for our processing
    int maxH = std::max(volume.h, volume.d);
    int maxW = std::max(volume.w, volume.h);
    Size s(maxW + volume.d, maxH);

    pre_compute_power_arrays(s);

    // for projections
    Mat x_y, y_z, x_z, xz_y, zx_y;

    if(show) {
        x_y = Mat(volume.h, volume.w, CV_32SC1, Scalar(0));
        y_z = Mat(volume.d, volume.h, CV_32SC1, Scalar(0));
        x_z = Mat(volume.d, volume.w, CV_32SC1, Scalar(0));
        xz_y = Mat(volume.h, volume.w + volume.d, CV_32SC1, Scalar(0));
        zx_y = Mat(volume.h, volume.w + volume.d, CV_32SC1, Scalar(0));
    }

    Moments3D m0, m1, m2, m3;
    double t_ocv = 100.0, t_dpm = 100.0, t_naive = 100.0, t_simd = 100.0;

    // If image type is uchar/char then optimise for 16/32 bit processing
    if(sizeof(T) == 1) {
        std::cout << "Runs:";
        for (int i=0; i < dim;i++) {
            std::cout << " " << (int)ceil(img.size(i) / 256.0);
        }
        std::cout << std::endl;
    }

    for(int i = 0; i < 10; i++)
    {
        CreateOperationalImages<T, int>(volume.cuboid);

        double t0 = (double)getTickCount();
        m3 = voxel_moments<T, int>(volume, order);
        t0 = ((double)getTickCount()-t0)/getTickFrequency();
        if(t0 < t_dpm) t_dpm = t0;

        if(sizeof(T) == 1) {
            Cuboid cuboid(volume.w, std::min(256, volume.h), std::min(256, volume.d));

            CreateOperationalImages<T, ushort>(cuboid);

            double tt = (double)getTickCount();
            m0 = voxel_moments<T, ushort>(volume, order);
            tt = ((double)getTickCount()-tt)/getTickFrequency();
            if(tt < t_simd) t_simd = tt;
        }

    }
        double t1 = (double)getTickCount();
        m1 = opencv_3dmoments<T>(volume, order);
        t1 = ((double)getTickCount()-t1)/getTickFrequency();
        if(t1 < t_ocv) t_ocv = t1;

        // if(i % 20 == 0) {
            double t2 = (double)getTickCount();
            m2 = straight_3dmoments<T>(volume, order);
            t2 = ((double)getTickCount()-t2)/getTickFrequency();
            if(t2 < t_naive) t_naive = t2;
        // }

    if(show) {
        projections<T>(volume, x_y, y_z, x_z, xz_y, zx_y, order);
        cout << "DPM_SIMD         DPM      OpenCV     Straight" << endl;
        DumpMoments(m0, m3, m1, m2);
    }

    cout << setprecision(6);
    if(sizeof(T) == 1) {
        cout << "SIMD DPM 3D : " << t_simd * 1000.0 << endl;
    }
    cout << "DPM 3D Total: " << t_dpm * 1000.0 << endl;
    cout << "OCV 3D Total: " << t_ocv * 1000.0 << endl;
    cout << "Naiive Total: " << t_naive * 1000.0 << endl;

    if(show) {
        ShowOff(x_y, "x_y");
        ShowOff(y_z, "y_z");
        ShowOff(x_z, "x_z");
        ShowOff(xz_y, "xz_y");
        if(order > 3)
            ShowOff(zx_y, "zx_y");
    }

    waitKey();
}

int main(int argc, char *argv[]) {

    std::string file;
    std::string nrrd_type;
    NRRD::Image<char> img_char;
    NRRD::Image<uchar> img_uchar(4, 5, 6, dummy);
    NRRD::Image<short> img_short;
    NRRD::Image<ushort> img_ushort;

    if(!SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS))
        cout << "Could not set process priority";

    if(!SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL))
        cout << "Could not set thread priority";

    if(argc == 1)
    {
        nrrd_type = "local";
    }
    else if(argc >= 5 && std::strcmp(argv[1], "-g") == 0)
    {
        nrrd_type = "local";
        img_uchar.set(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
        
        if(argc > 5) {
            order = atoi(argv[5]);
            if(order < 1 || order > 4) order = 3;
        }        
    }
    else {
        file = argv[1];
        nrrd_type = NRRD::getDataType(file);

        if(argc > 2) {
            order = atoi(argv[2]);
            if(order < 1 || order > 4) order = 3;
        }
    }

    if(nrrd_type == "") {
        std::cerr << "Failed to read file.nrrd.\n";
        return 1;

    }
    else if(nrrd_type != "local") {
        std::cout << "File: " << file << std::endl;
        std::cout << "Type: " << nrrd_type << std::endl;
    }

    // Loading an image
    if(nrrd_type == "char") img_char.load(file);
    if(nrrd_type == "unsigned char") img_uchar.load(file);
    if(nrrd_type == "short") img_short.load(file);
    if(nrrd_type == "unsigned short") img_ushort.load(file);

    if (nrrd_type != "local" && !img_char && !img_uchar && !img_short && !img_ushort) {
        std::cerr << "Failed to read file.nrrd.\n";
        return 1;
    }

    if(nrrd_type == "local")
        ProcessImage<uchar>(img_uchar, false);
    else if(nrrd_type == "char")
        ProcessImage<char>(img_char);
    else if(nrrd_type == "unsigned char") 
        ProcessImage<uchar>(img_uchar);
    else if(nrrd_type == "unsigned short") 
        ProcessImage<ushort>(img_ushort);
    else if(nrrd_type == "short") 
        ProcessImage<short>(img_short);

    return 0;
}
