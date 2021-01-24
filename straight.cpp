/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#include "moments.hpp"

template <typename T>
Moments3D straight_3dmoments(Volume<T> volume, int order)
{
    Moments3D m;

    if(order < 1 || order > 4) return m;

    for(int z = 0; z < volume.d; z++)
    { 
        for(int y = 0; y < volume.h; y++ )
        {   
            const T* p = &volume.voxels[volume.w * y + volume.w * volume.h * z];

            if(order == 0) {
                for(int x = 0; x < volume.w; x++ ) {
                    m.m000 += p[x];
                }
            }
            else if(order == 1) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x];
                    double yp = y * p[x];
                    double zp = z * p[x];

                    m.m000 += p[x];
                    m.m100 += xp;
                    m.m010 += yp;
                    m.m001 += zp;
                }
            }
            else if(order == 2) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = xp * x;
                    double yp = y * p[x], yyp = yp * y;
                    double zp = z * p[x], zzp = zp * z;

                    m.m000 += p[x];
                    m.m100 += xp;
                    m.m010 += yp;
                    m.m001 += zp;
                    m.m110 += xp * y;
                    m.m101 += xp * z;
                    m.m011 += yp * z;
                    m.m200 += xxp;
                    m.m020 += yyp;
                    m.m002 += zzp;
                }
            }
            else if(order == 3) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = xp * x, xxxp = xxp * x;
                    double yp = y * p[x], yy = y * y;
                    double zp = z * p[x], zz = z * z;
                    double xyp = xp * y;

                    m.m000 += p[x];
                    m.m100 += xp;
                    m.m010 += yp;
                    m.m001 += zp;
                    m.m110 += xyp;
                    m.m101 += xp * z;
                    m.m011 += yp * z;
                    m.m111 += xyp * z;
                    m.m200 += xxp;
                    m.m020 += yp * y;
                    m.m002 += zp * z;
                    m.m300 += xxxp;
                    m.m030 += yy * yp;
                    m.m003 += zz * zp;
                    m.m210 += xxp * y;
                    m.m201 += xxp * z;
                    m.m120 += xp * yy;
                    m.m021 += zp * yy;
                    m.m102 += xp * zz;
                    m.m012 += yp * zz;
                }
            }
            else if(order == 4) {
                for(int x = 0; x < volume.w; x++ ) {
                    double xp = x * p[x], xxp = xp * x, xxxp = xxp * x;
                    double yp = y * p[x], yyp = yp * y, yy = y * y, yyy = yy * y;
                    double zp = z * p[x], zzp = zp * z, zz = z * z, zzz = zz * z;
                    double xyp = xp * y, xyyp = xyp * y, xxyp = xxp * y;

                    m.m000 += p[x];
                    m.m100 += xp;
                    m.m010 += yp;
                    m.m001 += zp;
                    m.m110 += xyp;
                    m.m101 += xp * z;
                    m.m011 += yp * z;
                    m.m111 += xyp * z;
                    m.m200 += xxp;
                    m.m020 += yyp;
                    m.m002 += zzp;
                    m.m300 += xxxp;
                    m.m030 += yy * yp;
                    m.m003 += zz * zp;
                    m.m210 += xxyp;
                    m.m201 += xxp * z;
                    m.m120 += xyyp;
                    m.m021 += zp * yy;
                    m.m102 += xp * zz;
                    m.m012 += yp * zz;
                    m.m400 += xxxp * x;
                    m.m040 += yyy * yp;
                    m.m004 += zzz * zp;
                    m.m310 += xxxp * y;
                    m.m301 += xxxp * z;
                    m.m130 += xp * yyy;
                    m.m031 += zp * yyy;
                    m.m103 += xp * zzz;
                    m.m013 += yp * zzz;
                    m.m220 += xxp * yy;
                    m.m202 += xxp * zz;
                    m.m022 += yyp * zz;
                    m.m211 += xxyp * z;
                    m.m121 += xyyp * z;
                    m.m112 += xyp * zz;
                }
            }
        }
    }

    return m;
}

template Moments3D straight_3dmoments<char>(Volume<char>, int);
template Moments3D straight_3dmoments<uchar>(Volume<uchar>, int);
template Moments3D straight_3dmoments<short>(Volume<short>, int);
template Moments3D straight_3dmoments<ushort>(Volume<ushort>, int);
