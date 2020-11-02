// #pragma once


#ifndef __EIGEN__
#define __EIGEN__

class Eigen
{
public:
    static int PrincipalVector2X2(float sxx, float sxy, float syy, float* ux, float* uy);
    static int PrincipalValue2X2(float sxx, float sxy, float syy, float* lambda);
    static int EigenVector(float sxx, float syy, float sxy, float lambda, float* ux, float* uy);
};

#endif