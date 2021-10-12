#include <cmath>


#include "Gauss4point.hpp"


//Specialised constructor
Gauss4point::Gauss4point(double (*pFunction) (double), const double xmin, const double xmax)
{
    SetInterval(xmin, xmax);
    mpFunction = pFunction;
}


double Gauss4point::IntegrateFunction()
{

    //Output variable
    double IntValue = 0.0;

    //Vector containing the weights, adjusted for new interval
    Vector* weights
    = new Vector(4);

    //Vector of 4-point Gaussian nodes adjusted for interval
    Vector* xpoints
    = new Vector(4);



    //Fills vector with appropriate weights. I did the multiplying
    // outside the sum here, instead of at the end.
    (*weights)(1) = (mXmax-mXmin)*(18-sqrt(30))/72;
    (*weights)(2) = (mXmax-mXmin)*(18+sqrt(30))/72;
    (*weights)(3) = weights->Read(1);
    (*weights)(4) = weights->Read(0);

    //Fills vector with appropriate nodes
    (*xpoints)(1) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));
    (*xpoints)(2) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(3) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(4) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));


    //Applies standard quadrature rule
    for (int k=0; k<4;k++)
    {
        IntValue += weights->Read(k)*mpFunction(xpoints->Read(k));
    }

    //clear
    delete weights;
    delete xpoints;

    return IntValue;
}


//Very similar to function integration
double Gauss4point::IntegrateRHSProduct
(const int i, int npoints, Vector* pPoints)
{

    //Output variable
    double IntValue = 0.0;

    Vector* weights
    = new Vector(4);

    Vector* xpoints
    = new Vector(4);



    //Fills vector with weights
    (*weights)(1) = (mXmax-mXmin)*(18-sqrt(30))/72;
    (*weights)(2) = (mXmax-mXmin)*(18+sqrt(30))/72;
    (*weights)(3) = weights->Read(1);
    (*weights)(4) = weights->Read(0);

    //Essentially the same as before
    (*xpoints)(1) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));
    (*xpoints)(2) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(3) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(4) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));

    for (int k =0;k<4;k++)
    {
        //Now also multiplies by Lagrange Polynomial
        IntValue += weights->Read(k)*mpFunction(xpoints->Read(k))
        *EvaluateLagrangeBasis(xpoints->Read(k),i,npoints,pPoints);
    }

    //Clear
    delete weights;
    delete xpoints;


    return IntValue;
}

//Also similar to before, but now just uses Lagrange functions
double Gauss4point::IntegrateMatrixProduct
(const int i, const int j, int npoints, Vector* pPoints)
{

    //Output variable
    double IntValue = 0.0;

    Vector* weights
    = new Vector(4);

    Vector* xpoints
    = new Vector(4);



    (*weights)(1) = (mXmax-mXmin)*(18-sqrt(30))/72;
    (*weights)(2) = (mXmax-mXmin)*(18+sqrt(30))/72;
    (*weights)(3) = weights->Read(1);
    (*weights)(4) = weights->Read(0);


    (*xpoints)(1) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));
    (*xpoints)(2) = ((mXmax+mXmin)/2)- ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(3) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3-2*sqrt(6/5))/7));
    (*xpoints)(4) = ((mXmax+mXmin)/2)+ ((mXmax-mXmin)/2)*
    (sqrt((3+2*sqrt(6/5))/7));

        for (int k =0;k<4;k++)
    {
        IntValue += (weights->Read(k))
        *EvaluateLagrangeBasis(xpoints->Read(k),j,npoints,pPoints)
        *EvaluateLagrangeBasis(xpoints->Read(k),i,npoints,pPoints);
    }

    //Clear
    delete weights;
    delete xpoints;


    return IntValue;
}
