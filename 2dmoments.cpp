/* 
MIT License

Copyright (c) 2020 wild-ig 
*/
#include <vector>
#include <numeric>
#include <opencv2\opencv.hpp>
#include "moments.hpp"

using namespace cv;
using namespace std;

//power arrays - diagonal
double *d1, *d2, *d3, *d4;
//power arrays - anti-diagonal
double *a3, *a4;
//power array - slope
double *s4;
int crossing = 0;

// Compute the product sum of a vector and a power array
double product(const vector<long> &mat, const double power[], int many)
{
    double sum = 0.0;
    for(int i = 0; i < many; i++)
        sum += static_cast<double>(mat[i]) * power[i];

    return sum;
}

// Pre-compute the power arrays
void pre_compute_power_arrays(const Size s) {
    const int w = s.width;
    const int h = s.height;

    //power arrays
    d1 = new double [w + h];
    d2 = new double [w + h];
    d3 = new double [w + h];
    d4 = new double [w + h];
    a3 = new double [w + h];
    a4 = new double [w + h];
    s4 = new double [w + h * 2];

    for (int k=0; k< w + h; ++k)
    {
        d1[k] = k;
        double k2 = static_cast<double>(k) * static_cast<double>(k);
        d2[k] = k2;
        d3[k] = k2 * static_cast<double>(k);
        d4[k] = k2 * k2;
        // Anti-diagonal
        a3[k] = pow(static_cast<double>(k - h + 1), 3);
        a4[k] = pow(static_cast<double>(k - h + 1), 4);
    }

    // remember the crossing point
    crossing = h;

    // Slope based power array
    for (int k = 0; k < w + h * 2; ++k)
    {
        s4[k] = pow(static_cast<double>(k), 4);
    }
}

template <typename U>
Moments4th drt_moments(const Mat& image, int order)
{
    Size s = image.size();
    const int width = s.width;
    const int height = s.height;

    Moments4th m;

    if(order < 1 || order > 4) return m;

    // projection arrays
    vector<long> vert(width, 0);
    vector<long> horz(height, 0);
    vector<long> diag(width+height, 0);
    vector<long> anti(width+height, 0);
    vector<long> x_2y(width+height*2, 0);

    long* hptr = &horz[0];
    long* vptr = &vert[0];
    long* dptr = &diag[0];
    long* aptr = &anti[height - 1];
    long* x2yptr = &x_2y[0];

    for (int i = 0; i < height; i++)
    {
        const U* p = image.ptr<U>(i);

        if(order == 1) {
            for(int j = 0; j < width; j++)
            {
                vptr[j] += p[j];
                hptr[i] += p[j];
            }
        }
        else if(order == 2) {
            for(int j = 0; j < width; j++)
            {
                vptr[j] += p[j];
                hptr[i] += p[j];
                dptr[j] += p[j];
            }

            dptr++;
        }
        else if(order == 3) {
            for(int j = 0; j < width; j++)
            {
                vptr[j] += p[j];
                hptr[i] += p[j];
                dptr[j] += p[j];
                aptr[j] += p[j];
            }

            dptr++;
            aptr--;
        }
        else if(order == 4) {
            for(int j = 0; j < width; j++)
            {
                vptr[j] += p[j];
                hptr[i] += p[j];
                dptr[j] += p[j];
                aptr[j] += p[j];
                x2yptr[j] += p[j];
            }

            x2yptr+=2;
            dptr++;
            aptr--;
        }
    }

    m.m00 = accumulate(begin(vert), end(vert), 0.0);
    m.m10 = product(vert, d1, width);
    m.m01 = product(horz, d1, height);

    if(order > 1) {
        m.m20 = product(vert, d2, width);
        m.m02 = product(horz, d2, height);
        m.m11 = (product(diag, d2, width+height) - m.m02 - m.m20) / 2.0;
    }
    if(order > 2) {
        m.m30 = product(vert, d3, width);
        m.m03 = product(horz, d3, height);
        double temp_1 = product(diag, d3, width+height) / 6.0;
        double temp_2 = product(anti, &a3[crossing - height], width+height) / 6.0;
        m.m12 = temp_1 + temp_2 - m.m30/3.0;
        m.m21 = temp_1 - temp_2 - m.m03/3.0;
    }
    if(order > 3) {
        m.m40 = product(vert, d4, width);
        m.m04 = product(horz, d4, height);

        // 4th order diagonal and anti-diagonal projection moments
        double md_4 = product(diag, d4, width+height);
        double ma_4 = product(anti, &a4[crossing - height], width+height);

        m.m22 = (md_4 + ma_4 - 2*m.m40 - 2*m.m04) / 12.0;

        // 4th order moment along x+2y projection
        double ms_4 = product(x_2y, s4, width+height*2);

        m.m13 = (ms_4 - 2*md_4 + m.m40 - 14.0 * m.m04 - 12.0  * m.m22) / 24.0;
        m.m31 = (ms_4 - m.m40 - 24 * m.m22 - 32*m.m13 - 16*m.m04) / 8.0;
    }

    return m;
}

template Moments4th drt_moments<char>(const Mat& , int );
template Moments4th drt_moments<uchar>(const Mat& , int );
template Moments4th drt_moments<short>(const Mat& , int );
template Moments4th drt_moments<ushort>(const Mat& , int );
template Moments4th drt_moments<int>(const Mat& , int );
template Moments4th drt_moments<uint>(const Mat& , int );
