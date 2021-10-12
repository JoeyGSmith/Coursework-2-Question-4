#include <cassert>
#include <iomanip>
#include <fstream>

#include "AbstractApproximator.hpp"


//set number of points
void AbstractApproximator::SetNoPoints(const int Npoints)
{
    mNpoints = Npoints;
}


//set interval
void AbstractApproximator::SetInterval(const double xmin, const double xmax)
{
  assert(xmax>xmin);
  mXmin = xmin;
  mXmax = xmax;
}

//set function
void AbstractApproximator::Setfunction(double (*pFunction)(double))
{
    mpFunction = pFunction;
}

void AbstractApproximator::SaveFile(std::string fileName,Vector* xPoints ,
Vector* Exact, Vector* Approx)
{
    int nxvalue = xPoints->GetSize();
    std::ofstream writeFile;
    writeFile.open(fileName);
    assert(writeFile.is_open());

      // Write data
    for (int i=0; i<nxvalue; i++)
    {
        writeFile << std::setw(15) << xPoints->Read(i) << std::endl;
    }

        for (int i=0; i<nxvalue; i++)
    {
        writeFile << std::setw(15) << Exact->Read(i) << std::endl;
    }

        for (int i=0; i<nxvalue; i++)
    {
        writeFile << std::setw(15) << Approx->Read(i) << std::endl;
    }
    // Close file
    writeFile.close();



}
