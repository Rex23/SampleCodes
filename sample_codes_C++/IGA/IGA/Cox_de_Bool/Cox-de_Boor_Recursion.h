#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

int Cox_de_Boor_Recursion(const int& k, const std::vector <double>& Knot_Vector, std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients);

int Cox_de_Boor_Recursion_Single(const int& k, const std::vector <std::vector <std::vector <double> > >& Coefficients,
	const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients_New);

int PrintSplines(const int& Mode, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients);

int Paraview_Visualization(const std::string& FileName, const int& Intervals, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients);