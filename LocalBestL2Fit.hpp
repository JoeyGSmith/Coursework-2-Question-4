#ifndef BESTL2FITHEADERDEF
#define BESTL2FITHEADERDEF

#include "AbstractApproximator.hpp"
#include "AbstractQuadratureRule.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"


class LocalBestL2Fit:
    public AbstractApproximator
{

public:

    //specialised constructor
    LocalBestL2Fit(double (*pFunction) (double), AbstractQuadratureRule* pIntegrator,
              const double xmin, const double xmax, const int npoints,
               const int nintervals,  const std::string outputFileName);

    //approximation method
    void Approximate(const int nxvalues);

    void GaussianElimination(Matrix* p_A, Vector* p_P, Vector* p_f);

private:

    //hidden default constructor
    LocalBestL2Fit();

    AbstractQuadratureRule* mpIntegrator;

    //number of intervals
    int mnIntervals;


};

#endif
