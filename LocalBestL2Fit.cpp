#include <iostream>
#include <cassert>

#include "LocalBestL2Fit.hpp"



//Specialised constructor
LocalBestL2Fit::LocalBestL2Fit(double (*pFunction) (double), AbstractQuadratureRule* pIntegrator,
const double xmin, const double xmax, const int npoints, const int nintervals,
const std::string outputFileName)
{
        //Initialise
    SetInterval(xmin,xmax);
    assert(npoints > 1);
    SetNoPoints(npoints);
    Setfunction(pFunction);
    mOutputFileName = outputFileName + ".dat";
    mpIntegrator = pIntegrator;
    mnIntervals = nintervals;
}

void LocalBestL2Fit::Approximate(const int nxvalues)
{

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Using LocalBestL2Fit method with the following values: " << std::endl;
    std::cout << "Lower interval bound = " << mXmin << std::endl;
    std::cout << "Upper interval bound = " << mXmax << std::endl;
    std::cout << "Number of intervals = " << mnIntervals << std::endl;
    std::cout << "Number of points in each interval = " << nxvalues << std::endl;
    std::cout << "Degree of polynomials = " << mNpoints-1 << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;


    Vector* pPoints
    = new Vector(mNpoints);

    //Matrix of integrated Lagrange polynomials
    Matrix* A
    = new Matrix(mNpoints,mNpoints);

    //Vector for coefficients of approximated polynomial
    Vector* p
    = new Vector(mNpoints);

    //integrated RHS of equation (5)
    Vector* f
    = new Vector(mNpoints);

        //Total number of points approximated. Used at end to fill vectors
    int totalpoints= mnIntervals*(nxvalues-1)+1;

    //Boundaries of each interval
    Vector* xintervalpoints
    = new Vector(mnIntervals+1);

    //Vector of points to be approximated
    Vector* xpoints
    = new Vector(totalpoints);

    //vector of approximated values
    Vector* approx
    = new Vector(totalpoints);

    Vector* fexact
    = new Vector(totalpoints);





    double stepsize = (mXmax-mXmin)/((double)(mNpoints)-1.0);

    double stepsize2 = (mXmax-mXmin)/((double)(mnIntervals));

    double stepsize3 = stepsize2/(double)(nxvalues-1);

    //variable for the value of the function at a given point,
    //as well as error calculators
    double functionvalue, MaxError, currentError;

    for (int i=0; i< mnIntervals+1;i++)
    {
        (*xintervalpoints)[i] = mXmin + i*stepsize2;
    }


    //constructs vector pPoints for Lagrange polynomials
    for (int i=0; i< mNpoints; i++)
    {
        (*pPoints)[i]= mXmin + i*stepsize;
    }


    for(int k=0; k < mnIntervals; k++)
    {
        mpIntegrator->SetInterval(xintervalpoints->Read(k)
                                   ,xintervalpoints->Read(k+1));

        //Fills the vector f and matrix A with appropriate integrals
        for(int i =1;i < mNpoints+1;i++)
        {
            (*f)(i) = mpIntegrator->IntegrateRHSProduct(i-1,mNpoints,pPoints);

            for(int j=1;j < mNpoints+1;j++)
            {
                (*A)(i,j) = mpIntegrator->IntegrateMatrixProduct(i-1,j-1,mNpoints,pPoints);
            }
        }

        //This will fill the p vector with coefficients
        GaussianElimination(A,p,f);


        for (int i=0;i<nxvalues;i++)
        {
            (*xpoints)[k*(nxvalues-1)+i] = xintervalpoints->Read(k) + i*stepsize3;
            (*fexact)[k*(nxvalues-1)+i] = mpFunction(xpoints->Read(k*(nxvalues-1)+i));

            functionvalue = 0.0;

            for (int j=0;j<mNpoints;j++)
            {
                functionvalue += p->Read(j)*
                mpIntegrator->EvaluateLagrangeBasis(xpoints->Read(k*(nxvalues-1)+i),j,mNpoints,pPoints);
            }


            (*approx)[k*(nxvalues-1)+i] = functionvalue;
            currentError = approx->Read(k*(nxvalues-1)+i)-fexact->Read(k*(nxvalues-1)+i);


            if ( MaxError < std::abs(currentError))
            {
                MaxError = std::abs(currentError);
            }

    }




    }



    std::cout << "Function approximated with maximum error " << MaxError
    << std::endl;

    SaveFile(mOutputFileName,xpoints,fexact,approx);




    delete pPoints;
    delete A;
    delete p;
    delete f;
    delete pPoints;
    delete approx;
    delete fexact;
}

//This is the GaussianElimination function from question 2. I did not write this.
void LocalBestL2Fit::GaussianElimination(Matrix* p_A, Vector* p_p, Vector* p_f)
{
  int systemSize = p_p->GetSize();

  // Forward elimination step
  for (int i=1; i<systemSize; i++)
  {
    // Elimination
    double diagonal = 1.0 / (*p_A)(i,i);
    for (int j=i+1; j<=systemSize; j++)
    {
      double multiplier = diagonal * (*p_A)(j,i);
      for (int k=i; k<=systemSize; k++)
      {
        (*p_A)(j,k) = (*p_A)(j,k) - multiplier * (*p_A)(i,k);
      }
      (*p_f)(j) = (*p_f)(j) - multiplier * (*p_f)(i);
    }
  }

  // Backward substitution step
  (*p_p)(systemSize) = (*p_f)(systemSize) / (*p_A)(systemSize,systemSize);
  for (int i=systemSize-1; i>=1; i--)
  {
    (*p_p)(i) = (*p_f)(i);
    for (int k=i+1; k<=systemSize; k++)
    {
      (*p_p)(i) = (*p_p)(i) - (*p_A)(i,k) * (*p_p)(k);
    }
    (*p_p)(i) = (*p_p)(i) / (*p_A)(i,i);
  }
}
