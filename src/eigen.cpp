
#include <math.h>
#include <eigen.h>





int Eigen::PrincipalVector2X2(float sxx, float sxy, float syy, float *ux, float *uy)
{
	float norm, lambda;
	norm = (float)sqrtf(sxx * sxx + syy * syy + 2.0f * (sxy * sxy));
	if (norm > 0.0f)
	{
		sxx /= norm;
		sxy /= norm;
		syy /= norm;
	}
	
    PrincipalValue2X2(sxx, sxy, syy, &lambda);
    EigenVector(sxx, syy, sxy, lambda, ux, uy);
	return 0;
}

int  Eigen::PrincipalValue2X2(float sxx, float sxy, float syy, float *lambda)
{
  float delta,  epsilon = 1e-10f, sq, r;
  delta = (sxx - syy)*(sxx - syy) + 4.0f*sxy*sxy;
  if (delta < -epsilon) *lambda = 0.0f;
  else {
    sq = (float)sqrt(fabsf(delta)); 
    r = sxx + syy; 
    *lambda = (r + sq) / 2.0f; 
    }
return 0;
}

int  Eigen::EigenVector(float sxx, float syy, float sxy, float lambda, float *ux, float *uy)
{
        float det; 
        det = sqrtf(sxy * sxy + (sxx - lambda) * (sxx - lambda));
        if (det > 0.0f) \
        { 
            *ux = -sxy / det; 
            *uy = (sxx - lambda) / det; 
        } 
        else 
        { 
            det = sqrtf(sxy * sxy + (syy - lambda) * (syy - lambda)); 
            if (det > 0.0f) 
            { 
                *ux = (syy - lambda) / det; 
                *uy = -sxy / det; 
            } 
        }
        return 0;
}


