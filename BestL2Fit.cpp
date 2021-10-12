#include <iostream>

#include "BestL2Fit.hpp"



//Specialised constructor
BestL2Fit::BestL2Fit(double (*pFunction) (double), AbstractQuadratureRule* pIntegrator,
const double xmin, const double xmax, const int npoints, const std::string outputFileName)
{
        //Initialise
    SetInterval(xmin,xmax);
    SetNoPoints(npoints);
    Setfunction(pFunction);
    mOutputFileName = outputFileName + ".dat";
    mpIntegrator = pIntegrator;
}

void BestL2Fit::Approximate(const int nxvalues)
{

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Using BestL2Fit method with the following values: " << std::endl;
    std::cout << "Lower interval bound = " << mXmin << std::endl;
    std::cout << "Upper interval bound = " << mXmax << std::endl;
    std::cout << "Number of points = " << nxvalues << std::endl;
    std::cout << "Degree of polynomials = " << mNpoints-1 << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;


    Vector* pPoints
    = new Vector(mNpoints);



    double stepsize = (mXmax-mXmin)/((double)(mNpoints)-1.0);

    double stepsize2 = (mXmax-mXmin)/((double)(nxvalues)-1.0);

    //variable for the value of the function at a given point
    double functionvalue, MaxError, currentError;


    //constructs vector pPoints for Lagrange polynomials
    for (int i=0; i< mNpoints; i++)
    {
        (*pPoints)[i]= mXmin + i*stepsize;
    }


    //Matrix of integrated Lagrange polynomials
    Matrix* A
    = new Matrix(mNpoints,mNpoints);


    //Vector for coefficients of approximated polynomial
    Vector* p
    = new Vector(mNpoints);



    //integrated RHS of equation (5)
    Vector* f
    = new Vector(mNpoints);

    Vector* xpoints
    = new Vector(nxvalues);

    //vector of approximated values
    Vector* approx
    = new Vector(nxvalues);

    Vector* fexact
    = new Vector(nxvalues);



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
        (*xpoints)[i] = mXmin + i*stepsize2;
        (*fexact)[i] = mpFunction(xpoints->Read(i));

        functionvalue = 0.0;

        for (int j=0;j<mNpoints;j++)
        {
            functionvalue += p->Read(j)*
            mpIntegrator->EvaluateLagrangeBasis(xpoints->Read(i),j,mNpoints,pPoints);
        }


        (*approx)[i] = functionvalue;
        currentError = approx->Read(i)-fexact->Read(i);


        if ( MaxError < std::abs(currentError))
        {
            MaxError = std::abs(currentError);
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
void BestL2Fit::GaussianElimination(Matrix* p_A, Vector* p_p, Vector* p_f)
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
