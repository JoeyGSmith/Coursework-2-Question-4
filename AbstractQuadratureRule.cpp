#include <cassert>

#include "AbstractQuadratureRule.hpp"


double AbstractQuadratureRule::EvaluateLagrangeBasis(const double xval,
const int j, const int npoints, const Vector* pPoints) const
{

    //Makes sure the number of points is the same as the size of the Vector.
    assert (npoints == pPoints->GetSize());

    //Initialises output variable
    double Lvalue = 1.0;


    for (int i=0; i<npoints;i++)
    {
        if (i == j)
        {
            continue; //An undefined number is hard to integrate.
        }
        else
        {
            Lvalue *= (xval- pPoints->Read(i))/(pPoints->Read(j)-pPoints->Read(i));
        }

    }

    return Lvalue;
}

//Sets interval values
void AbstractQuadratureRule::SetInterval(const double xmin, const double xmax)
{
  assert(xmax>xmin);
  mXmin = xmin;
  mXmax = xmax;
}
