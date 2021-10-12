#ifndef ABSTRACTDIFFERENTIATORHEADERDEF
#define ABSTRACTDIFFERENTIATORHEADERDEF

#include <iostream>
#include "Vector.hpp"


class AbstractApproximator
{
    public:

        //Pure virtual method
        virtual void Approximate (const int nxvalues) = 0;

        //Sets the number of points n
        void SetNoPoints(const int Npoints);

        //Sets interval
        void SetInterval(const double Xmin, const double Xmax);

        //Sets the function to be approximated
        void Setfunction(double (*pFunction) (double));

        //Function for saving a vector to a file

        void SaveFile(std::string fileName,Vector* xPoints ,Vector* Exact, Vector* Approx);

    protected:

        double mXmin, mXmax;
        int mNpoints;

        double (*mpFunction) (double);

        //Output file name
        std::string mOutputFileName;


};

#endif
