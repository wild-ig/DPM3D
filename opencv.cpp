/* 
MIT License

Copyright (c) 2020 wild-ig 
*/

#include "moments.hpp"

template <typename T>
Moments3D opencv_3dmoments(Volume<T> volume, int order)
{
    Moments3D m;

    if(order < 0 || order > 4) return m;

    for(int z = 0; z < volume.d; z++)
    { 
        double x0y0 = 0.0;
        double x0y1 = 0.0;
        double x0y2 = 0.0;
        double x0y3 = 0.0;
        double x0y4 = 0.0;
        double x1y0 = 0.0;
        double x2y0 = 0.0;
        double x3y0 = 0.0;
        double x4y0 = 0.0;
        double x1y1 = 0.0;
        double x2y2 = 0.0;
        double x2y1 = 0.0;
        double x1y2 = 0.0;
        double x3y1 = 0.0;
        double x1y3 = 0.0;

        for(int y = 0; y < volume.h; y++ )
        {
            const T* p = &volume.voxels[volume.w * y + volume.w * volume.h * z];
            double x0 = 0;   // Ep
            double x1 = 0.0; // Ex*p
            double x2 = 0.0; // Exx*p
            double x3 = 0.0; // Exxx*p
            double x4 = 0.0; // Exxx*p

            if(order == 0) {
                for(int x = 0; x < volume.w; x++ ) {
                    x0 += p[x];
                }
            }
            else if(order == 1) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x];

                    x0 += p[x];
                    x1 += xp;
                }
            }
            else if(order == 2) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = x * xp;

                    x0 += p[x];
                    x1 += xp;
                    x2 += xxp;
                }
            }
            else if(order == 3) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = x * xp;

                    x0 += p[x];
                    x1 += xp;
                    x2 += xxp;
                    x3 += xxp * x;
                }
            }
            else if(order == 4) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = x * xp, xxxp = xxp * x;

                    x0 += p[x];
                    x1 += xp;
                    x2 += xxp;
                    x3 += xxxp;
                    x4 += xxxp * x;
                }
            }

            double py = y * x0, sy = y * y, cy = sy * y;

            switch(order) {
            case 4:
                x0y4 += py * cy;
                x3y1 += x3 * y;
                x1y3 += x1 * cy;
                x4y0 += x4;
                x2y2 += x2 * sy;
            case 3:
                x0y3 += py * sy;
                x2y1 += x2 * y;
                x1y2 += x1 * sy;
                x3y0 += x3;
            case 2:
                x0y2 += x0 * sy;
                x1y1 += x1 * y;
                x2y0 += x2;
            case 1:
                x0y1 += py;
                x1y0 += x1;
            case 0:
                x0y0 += x0;
            default:
                break;
            }
        }

        double sz = z * z, cz = sz * z, pz = x0y0 * z;

        switch(order) {
        case 4:
            m.m400 += x4y0;
            m.m040 += x0y4;
            m.m004 += pz * cz;
            m.m310 += x3y1;
            m.m130 += x1y3;
            m.m301 += x3y0 * z;
            m.m103 += x1y0 * cz;
            m.m031 += x0y3 * z;
            m.m013 += x0y1 * cz;
            m.m220 += x2y2;
            m.m202 += x2y0 * sz;
            m.m022 += x0y2 * sz;
            m.m211 += x2y1 * z;
            m.m121 += x1y2 * z;
            m.m112 += x1y1 * sz;
        case 3:
            m.m300 += x3y0;
            m.m030 += x0y3;
            m.m003 += pz * sz;
            m.m210 += x2y1;
            m.m120 += x1y2;
            m.m201 += x2y0 * z;
            m.m102 += x1y0 * sz;
            m.m021 += x0y2 * z;
            m.m012 += x0y1 * sz;
            m.m111 += x1y1 * z;
        case 2:
            m.m200 += x2y0;
            m.m020 += x0y2;
            m.m002 += x0y0 * sz;
            m.m110 += x1y1;
            m.m101 += x1y0 * z;
            m.m011 += x0y1 * z;
        case 1:
            m.m100 += x1y0;
            m.m010 += x0y1;
            m.m001 += pz;
        case 0:
            m.m000 += x0y0;
        default:
            break;
        }
    }

    return m;
}

template Moments3D opencv_3dmoments<char>(Volume<char>, int);
template Moments3D opencv_3dmoments<uchar>(Volume<uchar>, int);
template Moments3D opencv_3dmoments<short>(Volume<short>, int);
template Moments3D opencv_3dmoments<ushort>(Volume<ushort>, int);

