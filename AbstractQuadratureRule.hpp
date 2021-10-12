#ifndef ABSTRACTQUADRATURERULEHEADERDEF
#define ABSTRACTQUADRATURERULEHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class AbstractQuadratureRule
{
    public:


        //Pure virtual methods for calculating various integrals.
         virtual double IntegrateFunction() = 0;

         virtual double IntegrateRHSProduct
            (const int i, int npoints, Vector* pPoints) = 0;

         virtual double IntegrateMatrixProduct
            (const int i, const int j, int npoints, Vector* pPoints) = 0;

        void SetInterval(const double xmin, const double xmax);


        //I made this function public to be used by my
        //BestL2Fit class
        double EvaluateLagrangeBasis(const double xval, const int j,
        const int npoints, const Vector* pPoints) const;


    protected:



        double mXmin;
        double mXmax;

        double (*mpFunction) (double);

        int mNoValues;



};
#endif
