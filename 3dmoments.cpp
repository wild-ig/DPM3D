/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#include <vector>
#include <numeric>
#include <opencv2\opencv.hpp>
#include "moments.hpp"

using namespace cv;

// Processing/operational matrices
Mat XY, YZ, XZ, XZY, X_ZY;

void ResetOperationalImages() {
    XY = 0;
    YZ = 0;
    XZ = 0;
    XZY = 0;
    X_ZY = 0;
}

template <typename T, typename U>
void CreateOperationalImages(Cuboid c) {
    int cv_type = CV_32SC1;
    if(sizeof(T) == 1 && sizeof(U) == 2)
    {
        cv_type = CV_16UC1;
        c.h = std::min(c.h, 256);
        c.d = std::min(c.d, 256);
    }

    XY = Mat(c.h, c.w, cv_type, Scalar(0));
    // Always set to integer
    YZ = Mat(c.d, c.h, CV_32SC1, Scalar(0));
    XZ = Mat(c.d, c.w, cv_type, Scalar(0));
    XZY = Mat(c.h, c.w + c.d, cv_type, Scalar(0));
    X_ZY = Mat(c.h, c.w + c.d, cv_type, Scalar(0));
}

template void CreateOperationalImages<uchar, ushort>(Cuboid);
template void CreateOperationalImages<uchar, int>(Cuboid);
template void CreateOperationalImages<char, ushort>(Cuboid);
template void CreateOperationalImages<char, int>(Cuboid);
template void CreateOperationalImages<short, int>(Cuboid);
template void CreateOperationalImages<ushort, int>(Cuboid);
template void CreateOperationalImages<short, ushort>(Cuboid);
template void CreateOperationalImages<ushort, ushort>(Cuboid);


// Adjust for origin
Moments3D origin_compesation(Moments3D m, Origin origin, int order) {
    Moments3D mo;
    int x = origin.x, y = origin.y, z = origin.z;

    double x2 = x * x, x3 = x2 * x, x4 = x3 * x;
    double y2 = y * y, y3 = y2 * y, y4 = y3 * y;
    double z2 = z * z, z3 = z2 * z, z4 = z3 * z;
    double xy = x * y, xz = x * z, yz = y * z;
    double x2y = x2 * y, x2z = x2 * z;
    double xy2 = y2 * x, xz2 = z2 * x;
    double y2z = y2 * z, yz2 = z2 * y;
    double x2y2 = x2 * y2, x2z2 = x2 * z2, y2z2 = y2 * z2;
    double xyz = x * y * z;

    mo.m000 = m.m000;
    mo.m100 = m.m100 + m.m000*x;
    mo.m010 = m.m010 + m.m000*y;
    mo.m001 = m.m001 + m.m000*z;

    if(order > 1) {
        mo.m200 = m.m200 + 2*m.m100*x + m.m000*x2;
        mo.m020 = m.m020 + 2*m.m010*y + m.m000*y2;
        mo.m002 = m.m002 + 2*m.m001*z + m.m000*z2;

        mo.m110 = m.m110 + m.m100*y + m.m010*x + m.m000*xy;
        mo.m011 = m.m011 + m.m010*z + m.m001*y + m.m000*yz;
        mo.m101 = m.m101 + m.m100*z + m.m001*x + m.m000*xz;
    }

    if(order > 2) {
        // (x+xbar)^3
        // (x^2 + 2x.xbar + xbar^2)(x+xbar)
        // m300 + 3M200.xbar + 3M100.xbar^2 + xbar^3
        mo.m300 = m.m300 + 3*m.m200*x + 3*m.m100*x2 + m.m000*x3;
        mo.m030 = m.m030 + 3*m.m020*y + 3*m.m010*y2 + m.m000*y3;
        mo.m003 = m.m003 + 3*m.m002*z + 3*m.m001*z2 + m.m000*z3;

        // (x+xbar)^2(y+ybar)
        // (x^2 + 2x.xbar + xbar^2)(y+ybar)
        // m210 = M210 + M200.ybar + 2M110.xbar + 2M100.xbar.ybar + M010.xbar^2 +xbar^2.ybar
        mo.m210 = m.m210 + m.m200*y + 2*m.m110*x + 2*m.m100*xy + m.m010*x2 + m.m000*x2y;
        mo.m201 = m.m201 + m.m200*z + 2*m.m101*x + 2*m.m100*xz + m.m001*x2 + m.m000*x2z;
        mo.m120 = m.m120 + m.m020*x + 2*m.m110*y + 2*m.m010*xy + m.m100*y2 + m.m000*xy2;
        mo.m021 = m.m021 + m.m020*z + 2*m.m011*y + 2*m.m010*yz + m.m001*y2 + m.m000*y2z;
        mo.m102 = m.m102 + m.m002*x + 2*m.m101*z + 2*m.m001*xz + m.m100*z2 + m.m000*xz2;
        mo.m012 = m.m012 + m.m002*y + 2*m.m011*z + 2*m.m001*yz + m.m010*z2 + m.m000*yz2;
    
        // (x+xbar)(y+ybar)(x+zbar)
        // (xy + x.ybar + y.xbar + xbar.ybar)(z+zbar)
        // m111 = M111 + M101.ybar + M011.xbar + M001.xbar.ybar + M110.zbar + M100.ybar.zbar + M010.xbar.zbar + xbar.ybar.zbar
        mo.m111 = m.m111 + m.m101*y + m.m011*x + m.m110*z + m.m001*xy + m.m100*yz + m.m010*xz + xyz*m.m000;
    }

    if(order > 3) {
        // (x+xbar)^4
        // (x^2 + 2x.xbar + xbar^2)(x^2 + 2x.xbar + xbar^2)
        // x^4 + 4x^3.xbar + 6x^2.xbar^2 + 4x.xbar^3 + xbar^4
        // m400 = M400 + 4M300.x + 6M200.x.x + 4.M100.x^3 + M000.x^4
        mo.m400 = m.m400 + 4*m.m300*x + 6*m.m200*x2 + 4*m.m100*x3 + m.m000*x4;
        mo.m040 = m.m040 + 4*m.m030*y + 6*m.m020*y2 + 4*m.m010*y3 + m.m000*y4;
        mo.m004 = m.m004 + 4*m.m003*z + 6*m.m002*z2 + 4*m.m001*z3 + m.m000*z4;

        // (x+xbar)^2(y+ybar)(z+zbar)
        // (x^2 + 2x.xbar + xbar^2)(y+ybar)(z+zbar)
        // (x^2.y + x^2.ybar + 2xy.xbar + 2x.xbar.ybar + y.xbar^2 + xbar^2.ybar)(z+zbar)
        // x^2.y.z +x^2.z.ybar + 2xyz.xbar + 2xz.xbar.ybar + yz.xbar^2 + z.xbar^2.ybar
        //      +x^2.y.zbar + x^2.ybar.zbar + 2xy.xbar.zbar + 2x.xbar.ybar.zbar + y.xbar^2.zbar + xbar^2.ybar.zbar
        // m211 = M211 + M201.ybar + 2M111.xbar + 2M101.xbar.ybar + M011.xbar^2 + M001xbar^2.ybar
        //      + M210.zbar + M200.ybar.zbar + 2M110.xbar.zbar + 2M100.xbar.ybar.zbar + M010.xbar^2.zbar + M000.xbar^2.ybar.zbar
        mo.m211 = m.m211 + m.m201*y + 2*m.m111*x + 2*m.m101*xy + m.m011*x2 + m.m001*x2y +
                    m.m210*z + m.m200*yz + 2*m.m110*xz + 2*m.m100*xyz + m.m010*x2z + m.m000*x*xyz;
        mo.m121 = m.m121 + m.m021*x + 2*m.m111*y + 2*m.m011*xy + m.m101*y2 + m.m001*xy2 +
                    m.m120*z + m.m020*xz + 2*m.m110*yz + 2*m.m010*xyz + m.m100*y2z + m.m000*xyz*y;
        mo.m112 = m.m112 + m.m012*x + 2*m.m111*z + 2*m.m011*xz + m.m110*z2 + m.m010*xz2 +
                    m.m102*y + m.m002*xy + 2*m.m101*yz + 2*m.m001*xyz + m.m100*yz2 + m.m000*xyz*z;

        // (x+xb)^2(y+yb)^2
        // (x^2 + 2x.xb + xb^2)(y^2 + 2y.yb + yb^2)
        // x^2y^2 + 2x^2y.yb + x^2yb^2 + 2xy^2.xb + 4xy.xb.yb + 2x.xb.yb^2 + y^2.xb^2 + 2y.xb^2.yb + xb^2.yb^2
        // M220 + 2M210.y + M200.y.y + 2M120.x + 4M110.x.y + 2M100.x.y.y + M020.x.x + 2M010.x.x.y + M000.x.x.y.y
        mo.m220 = m.m220 + 2*m.m210*y + m.m200*y2 + 2*m.m120*x + 4*m.m110*xy
                    + 2*m.m100*xy2 + m.m020*x2 + 2*m.m010*x2y + m.m000*x2y2;
        mo.m202 = m.m202 + 2*m.m201*z + m.m200*z2 + 2*m.m102*x + 4*m.m101*xz
                    + 2*m.m100*xz2 + m.m002*x2 + 2*m.m001*x2z + m.m000*x2z2;
        mo.m022 = m.m022 + 2*m.m021*z + m.m020*z2 + 2*m.m012*y + 4*m.m011*yz
                    + 2*m.m010*yz2 + m.m002*y2 + 2*m.m001*y2z + m.m000*y2z2;


        // (x^3 + 3x^2.a + 3x.a^2 + a^3)(y+b)
        // x^3y + x^3.b + 3x^2y.a +3x^2.a.b + 3xy.a^2 + 3x.a^2.b + y.a^3 + a^3.b
        // M310 + M300.y + 3*M210.x + 3*M200.x.y + 3*M110*x*x + 3M100*x*x*y + M010*x*x*x + x*x*x*y
        mo.m310 = m.m310 + m.m300*y + 3*m.m210*x + 3*m.m200*xy + 3*m.m110*x2 + 3*m.m100*x2y + m.m010*x3 + m.m000*x3*y;
        mo.m130 = m.m130 + m.m030*x + 3*m.m120*y + 3*m.m020*xy + 3*m.m110*y2 + 3*m.m010*xy2 + m.m100*y3 + m.m000*x*y3;
        mo.m301 = m.m301 + m.m300*z + 3*m.m201*x + 3*m.m200*xz + 3*m.m101*x2 + 3*m.m100*x2z + m.m001*x3 + m.m000*x3*z;
        mo.m103 = m.m103 + m.m003*x + 3*m.m102*z + 3*m.m002*xz + 3*m.m101*z2 + 3*m.m001*xz2 + m.m100*z3 + m.m000*x*z3;
        mo.m031 = m.m031 + m.m030*z + 3*m.m021*y + 3*m.m020*yz + 3*m.m011*y2 + 3*m.m010*y2z + m.m001*y3 + m.m000*y3*z;
        mo.m013 = m.m013 + m.m003*y + 3*m.m012*z + 3*m.m002*yz + 3*m.m011*z2 + 3*m.m001*yz2 + m.m010*z3 + m.m000*y*z3;
        // (x+xbar)(y+ybar)(x+zbar)
        // (xy + x.ybar + y.xbar + xbar.ybar)(z+zbar)
        // m111 = M111 + M101.ybar + M011.xbar + M001.xbar.ybar + M110.zbar + M100.ybar.zbar + M010.xbar.zbar + xbar.ybar.zbar
        // mo.m111 = m.m111 + m.m101*y + m.m011*x + m.m110*z + m.m001*x*y + m.m100*y*z + m.m010*x*z + x*y*z;
    }

    return mo;

}

// Compute the moments of each of the processing matrices
template <typename U>
Moments3D compute_moments(int d, int order) {
    Moments3D m;
    Moments4th m1, m2, m3, m4, m5;

    m1 = drt_moments<U>(XY, order);

    // Note YZ is always iteger
    m2 = drt_moments<int>(YZ, order);

    if(order > 1)
        m3 = drt_moments<U>(XZ, order);

    if(order > 2)
        m4 = drt_moments<U>(XZY, order);

    if(order > 3)
        m5 = drt_moments<U>(X_ZY, order);

    m.m000 = m1.m00;
    m.m100 = m1.m10;
    m.m010 = m1.m01;
    m.m001 = m2.m01;

    if(order > 1) {
        m.m200 = m1.m20;
        m.m020 = m1.m02;
        m.m002 = m2.m02;

        m.m110 = m1.m11;
        m.m011 = m2.m11;
        m.m101 = m3.m11;
    }

    if(order > 2) {
        m.m300 = m1.m30;
        m.m030 = m1.m03;
        m.m003 = m2.m03;

        m.m210 = m1.m21;
        m.m201 = m3.m21;
        m.m120 = m1.m12;
        m.m021 = m2.m21;
        m.m102 = m3.m12;
        m.m012 = m2.m12;
    
        m.m111 = (m4.m21 - m1.m21 - m2.m12)/2.0;
    }

    if(order > 3) {
        m.m400 = m1.m40;
        m.m040 = m1.m04;
        m.m004 = m2.m04;

        m.m310 = m1.m31;
        m.m301 = m3.m31;
        m.m130 = m1.m13;
        m.m031 = m2.m31;
        m.m103 = m3.m13;
        m.m013 = m2.m13;
    
        m.m220 = m1.m22;
        m.m202 = m3.m22;
        m.m022 = m2.m22;

        // Adjust x position (y = 0)
        double xbar = d - 1.0;
        //u31 = M31 - 3M21xbar + 3M11xbar**2 - M01xbar**3 - M30ybar + 3M20xbarybar - 3m10xbar**2ybar + xbar**3ybar
        m5.m31 = m5.m31 - 3*m5.m21*xbar + 3*m5.m11*xbar*xbar - m5.m01*xbar*xbar*xbar;
        //u13 = M13 - 3M12ybar + 3M11ybar**2 - M10ybar**3 - M03xbar + 3M02xbarybar - 3m01ybar**2xbar + ybar**3xbar
        m5.m13 = m5.m13 - m5.m03*xbar;
        
        m.m121 = (m4.m22 - m.m220 - m.m022) / 2.0;
        m.m211 = (m4.m31 - m5.m31 - 2 * m.m013) / 6.0;
        m.m112 = (m4.m31 + m5.m31 - 2 * m.m310) / 6.0;
    }

    return m;
}

// Compute/sum along the x direction (width)
template <typename T, typename U>
inline int compute_line(const T* p, U* pXY, U* pXZ, U* pXZY, U* pX_ZY, int w, int order) {
    int sumx = 0;

    if(order == 4) {
        for(int i = 0; i < w; i++) {
            pXY[i] += p[i];
            pXZ[i] += p[i];
            pXZY[i] += p[i];
            pX_ZY[i] += p[i];
            sumx += p[i];
        }
    }
    else if(order == 3) {
        for(int i = 0; i < w; i++) {
            pXY[i] += p[i];
            pXZ[i] += p[i];
            pXZY[i] += p[i];
            sumx += p[i];
        }
    }
    else if(order == 2) {
        for(int i = 0; i < w; i++) {
            pXY[i] += p[i];
            pXZ[i] += p[i];
            sumx += p[i];
        }
    }
    else if(order == 1) {
        for(int i = 0; i < w; i++) {
            pXY[i] += p[i];
            sumx += p[i];
        }
    }

    return sumx;
}


template <typename T, typename U>
Moments3D voxel_moments(Volume<T> volume, VOI voi, int order) {
    Moments3D m;
    int slice = volume.w * volume.h;
    int step = volume.w;

    // Limit check
    if(order < 1 || order > 4) return m;

    // Origin check - and bail if necessary
    if(voi.x >= volume.w) return m;
    if(voi.y >= volume.h) return m;
    if(voi.z >= volume.d) return m;

    // Bounds check and readjust
    if(voi.x + voi.w > volume.w) voi.w = volume.w - voi.x;
    if(voi.y + voi.h > volume.h) voi.h = volume.h - voi.y;
    if(voi.z + voi.d > volume.d) voi.d = volume.d - voi.z;

    // origin voxel
    const T* pxyz = &volume.voxels[voi.z * slice + voi.y * step + voi.x];

    // Generate projections
    for(int k = 0; k < voi.d; k++) {
        U* pXZ = XZ.ptr<U>(k);
        int* pYZ = YZ.ptr<int>(k);

        const T* p = pxyz;

        for(int j = 0; j < voi.h; j++) {

            U* pXY = XY.ptr<U>(j);
            U* pXZY = &XZY.ptr<U>(j)[k];
            U* pX_ZY = &X_ZY.ptr<U>(j)[voi.d-k-1];

            pYZ[j] += compute_line<T, U>(p, pXY, pXZ, pXZY, pX_ZY, voi.w, order);

            // Move to next "y" step
            p += step;
        }

        // Move to next "z" slice
        pxyz += slice;
    }

    m = compute_moments<U>(voi.d, order);

    // Compensate for origin (if non-zero)
    if(voi.x != 0 || voi.y != 0 || voi.z != 0) {
        m = origin_compesation(m, voi.origin, order);
    }

    return m;
}

template Moments3D voxel_moments<char, ushort>(Volume<char>, VOI, int);
template Moments3D voxel_moments<uchar, ushort>(Volume<uchar>, VOI, int);
template Moments3D voxel_moments<short, ushort>(Volume<short>, VOI, int);
template Moments3D voxel_moments<ushort, ushort>(Volume<ushort>, VOI, int);
template Moments3D voxel_moments<char, int>(Volume<char>, VOI, int);
template Moments3D voxel_moments<uchar, int>(Volume<uchar>, VOI, int);
template Moments3D voxel_moments<short, int>(Volume<short>, VOI, int);
template Moments3D voxel_moments<ushort, int>(Volume<ushort>, VOI, int);


template <typename T, typename U>
Moments3D voxel_moments(Volume<T> volume, int order) {
    Moments3D m;
    VOI voi;

    if(sizeof(T) == 1 && sizeof(U) == 2) {
        int depth_runs = (int)ceil(volume.d / 256.0);
        int height_runs = (int)ceil(volume.h / 256.0);

        Cuboid cuboid(volume.w, std::min(256, volume.h), std::min(256, volume.d));
        voi.cuboid = cuboid;

        for(int k = 0; k < depth_runs; k++) {
            for(int j = 0; j < height_runs; j++) {
                voi.x = 0;
                voi.y = j * 256;
                voi.z = k * 256;

                m += voxel_moments<T, U>(volume, voi, order);

                ResetOperationalImages();
            }
        }

        return m;
    }
    
    voi.origin = Origin(0, 0, 0);
    voi.cuboid = volume.cuboid;
    return voxel_moments<T, U>(volume, voi, order);
}

template Moments3D voxel_moments<char, ushort>(Volume<char>, int order);
template Moments3D voxel_moments<uchar, ushort>(Volume<uchar>, int order);
template Moments3D voxel_moments<short, ushort>(Volume<short>, int order);
template Moments3D voxel_moments<ushort, ushort>(Volume<ushort>, int order);
template Moments3D voxel_moments<char, int>(Volume<char>, int order);
template Moments3D voxel_moments<uchar, int>(Volume<uchar>, int order);
template Moments3D voxel_moments<short, int>(Volume<short>, int order);
template Moments3D voxel_moments<ushort, int>(Volume<ushort>, int order);
