/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#ifndef MOMEMTS
#define MOMENTS

#include <opencv2\opencv.hpp>

struct Moments3D {
    double m000 = 0.0;
    double m100 = 0.0, m200 = 0.0, m300 = 0.0, m400 = 0.0;
    double m010 = 0.0, m020 = 0.0, m030 = 0.0, m040 = 0.0;
    double m001 = 0.0, m002 = 0.0, m003 = 0.0, m004 = 0.0;
    double m110 = 0.0, m101 = 0.0, m011 = 0.0, m111 = 0.0;
    double m210 = 0.0, m201 = 0.0, m102 = 0.0, m012 = 0.0, m021 = 0.0, m120 = 0.0;
    double m220 = 0.0, m202 = 0.0, m022 = 0.0;
    double m112 = 0.0, m121 = 0.0, m211 = 0.0;
    double m310 = 0.0, m301 = 0.0, m103 = 0.0, m013 = 0.0, m031 = 0.0, m130 = 0.0;

    Moments3D& operator=(Moments3D& m) {
        m000 = m.m000;
        m100 = m.m100;
        m200 = m.m200;
        m300 = m.m300;
        m400 = m.m400;
        m010 = m.m010;
        m020 = m.m020;
        m030 = m.m030;
        m040 = m.m040;
        m001 = m.m001;
        m002 = m.m002;
        m003 = m.m003;
        m004 = m.m004;
        m110 = m.m110;
        m101 = m.m101;
        m011 = m.m011;
        m111 = m.m111;
        m210 = m.m210;
        m201 = m.m201;
        m102 = m.m102;
        m012 = m.m012;
        m021 = m.m021;
        m120 = m.m120;
        m220 = m.m220;
        m202 = m.m220;
        m022 = m.m022;
        m112 = m.m112;
        m121 = m.m121;
        m211 = m.m211;
        m310 = m.m310;
        m301 = m.m301;
        m103 = m.m103;
        m013 = m.m013;
        m031 = m.m031;
        m130 = m.m130;

        return *this;
    }

    Moments3D& operator+=(Moments3D& m) {
        m000 += m.m000;
        m100 += m.m100;
        m200 += m.m200;
        m300 += m.m300;
        m400 += m.m400;
        m010 += m.m010;
        m020 += m.m020;
        m030 += m.m030;
        m040 += m.m040;
        m001 += m.m001;
        m002 += m.m002;
        m003 += m.m003;
        m004 += m.m004;
        m110 += m.m110;
        m101 += m.m101;
        m011 += m.m011;
        m111 += m.m111;
        m210 += m.m210;
        m201 += m.m201;
        m102 += m.m102;
        m012 += m.m012;
        m021 += m.m021;
        m120 += m.m120;
        m220 += m.m220;
        m202 += m.m220;
        m022 += m.m022;
        m112 += m.m112;
        m121 += m.m121;
        m211 += m.m211;
        m310 += m.m310;
        m301 += m.m301;
        m103 += m.m103;
        m013 += m.m013;
        m031 += m.m031;
        m130 += m.m130;

        return *this;
    }
};

struct Moments4th {
    double m00 = 0.0, m01 = 0.0, m10 = 0.0, m11 = 0.0, m20 = 0.0, 
            m02 = 0.0, m30 = 0.0, m12 = 0.0, m21 = 0.0, m03 = 0.0,
            m04 = 0.0, m40 = 0.0, m13 = 0.0, m31 = 0.0, m22 = 0.0;
};


struct Origin {
    int x = 0, y = 0, z = 0;

    Origin() {}
    Origin(int x_origin, int y_origin, int z_origin) :x(x_origin), y(y_origin), z(z_origin) {}
    
    Origin& operator=(Origin& o) {
        x = o.x;
        y = o.y;
        z = o.z;

        return *this;
    }
};

struct Cuboid {
    int w=0, h=0, d=0;

    Cuboid() {}
    Cuboid(int width, int height, int depth) :w(width), h(height), d(depth) {}

    Cuboid& operator=(Cuboid& c) {
        w = c.w;
        h = c.h;
        d = c.d;

        return *this;
    }
};

// Volume of interest
struct VOI {
    union {
        Origin origin;
        struct {
            int x, y, z;
        };
    };
    union {
        Cuboid cuboid;
        struct {
            int w, h, d;
        };
    };
    
    VOI() {}
    VOI(Origin o, Cuboid c) :origin(o), cuboid(c) {}

    VOI& operator=(VOI& voi) {
        origin = voi.origin;
        cuboid = voi.cuboid;

        return *this;
    }
};

template <typename T> struct Volume {
    T* voxels;
    union {
        Cuboid cuboid;
        struct {
            int w, h, d;
        };
    };

    Volume() {}
};

void pre_compute_power_arrays(const cv::Size s);
template <typename T, typename U>
void CreateOperationalImages(Cuboid c);
void ResetOperationalImages();

// Compute the moments of an 3D ROI
template <typename T, typename U = int>
Moments3D voxel_moments(Volume<T> volume, VOI v, int order = 3);
template <typename T, typename U = int>
Moments3D voxel_moments(Volume<T> volume, int order=3);

// Compute the 2D moments of a 2D image
template <typename U>
Moments4th drt_moments(const cv::Mat& image, int order = 3);

// Generate the projections of a 3D image
template <typename T>
void projections(Volume<T> volume, cv::Mat& x_y, cv::Mat& y_z, cv::Mat& x_z, cv::Mat& xz_y, cv::Mat& zx_y, int order = 3);

// Compute the naive 3D moments
template <typename T>
Moments3D straight_3dmoments(Volume<T> volume, int order = 3);

// Compute open CV version of the 3d moments
template <typename T>
Moments3D opencv_3dmoments(Volume<T> volume, int order = 3);

#endif // MOMENTS
