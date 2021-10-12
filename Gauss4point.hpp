#ifndef GAUSS4POINTHEADERDEF
#define GAUSS4POINTHEADERDEF

#include "AbstractQuadratureRule.hpp"

class Gauss4point:
    public AbstractQuadratureRule
{

    public:

        //Specialised constructor
        Gauss4point(double (*pFunction) (double), const double xmin, const double xmax);

        double IntegrateFunction();

        double IntegrateRHSProduct
            (const int i, int npoints, Vector* pPoints);

        double IntegrateMatrixProduct
            (const int i, const int j, int npoints, Vector* pPoints);

    private:

        //Hidden default constructor
        Gauss4point();


};

#endif
