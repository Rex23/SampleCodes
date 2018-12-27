#include "Cox-de_Boor_Recursion.h"

int main(int argc, char* argv[])
{
	int k = 3;
	std::vector <double> Knot_Vector{ 0,0,0,1,1,3,3,3 };
	std::vector <std::vector <std::vector <double> > > Coefficients;
	std::vector <int> Initial_Intervals;

	Cox_de_Boor_Recursion(k, Knot_Vector, Initial_Intervals, Coefficients);

	PrintSplines(2, Knot_Vector, Initial_Intervals, Coefficients);
	
	std::string FileName = "BSplines.vtk";
	Paraview_Visualization(FileName, 1000, Knot_Vector, Initial_Intervals, Coefficients);

	system("pause");
	return 1;
}