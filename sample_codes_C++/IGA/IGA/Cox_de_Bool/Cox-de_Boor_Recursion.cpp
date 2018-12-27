#include "Cox-de_Boor_Recursion.h"

//Implemented for the Cox-de Boor recursion formula (See Pg. 44 of the book "An Introduction to NURBS with Historical Perspective")
//Knot vectors: 1, ... , p
//Coefficients: (n,p - 1,k)
//Goal: Obtain Coefficients for the (k+1) order

int Cox_de_Boor_Recursion(const int& k, const std::vector <double>& Knot_Vector, std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	//std::vector <std::vector <std::vector <double> > > Coefficients;
	std::vector <std::vector <std::vector <double> > > Coefficients_New;

	//The initial one:

	//int Num_of_Splines = 0;
	//std::vector <int> Initial_Intervals;
	for (size_t m = 0; m < Knot_Vector.size() - 1; m++) {
		if (Knot_Vector[m] < Knot_Vector[m + 1]) {
			//Num_of_Splines = m + 1;
			Initial_Intervals.push_back(m + 1);
		}
	}

	Coefficients.resize(Knot_Vector.size() - 1); //The number of splines for k = 1
	for (size_t m = 0; m < Coefficients.size(); m++) {
		Coefficients[m].resize(Knot_Vector.size() - 1); //The number of intervals
		for (size_t m1 = 0; m1 < Coefficients[m].size(); m1++) {
			Coefficients[m][m1].resize(1);
			//Coefficients[m][m1][0] = 1.0;
		}
	}

	for (size_t m = 0; m < Initial_Intervals.size(); m++) {
		Coefficients[Initial_Intervals[m] - 1][Initial_Intervals[m] - 1][0] = 1.0;
	}

	for (int m = 0; m < k - 1; m++) {
		Cox_de_Boor_Recursion_Single(m+1, Coefficients, Knot_Vector, Initial_Intervals, Coefficients_New);
		Coefficients = Coefficients_New;
	}

	return 1;
}

int Cox_de_Boor_Recursion_Single(const int& k, const std::vector <std::vector <std::vector <double> > >& Coefficients, 
	                      const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients_New)
{
	Coefficients_New.resize(Coefficients.size() - 1);
	for (size_t m = 0; m < Coefficients_New.size(); m++) {
		Coefficients_New[m].resize(Knot_Vector.size() - 1);
		for (size_t m1 = 0; m1 < Coefficients_New[m].size(); m1++) {
			Coefficients_New[m][m1].resize(Coefficients[m][m1].size() + 1);
		}
	}

	for (size_t m = 0; m < Coefficients.size() - 1; m++) { // (xi <= t < xi+1 )

		for (size_t mm = 0; mm < Knot_Vector.size() - 1; mm++) {
			
			double x_i = Knot_Vector[m];
			double x_i_1;
			double x_i_k;
			double x_i_k_1;

			x_i_1 = Knot_Vector[m + 1];
			x_i_k = Knot_Vector[m + k];
			x_i_k_1 = Knot_Vector[m + k + 1];

			for (size_t m1 = 0; m1 < Coefficients[m][mm].size(); m1++) { // Orders 0 to k - 1

				double t1, t2;
				if (x_i_k - x_i == 0.0)
					t1 = 0.0;
				else
					t1 = (-x_i) / (x_i_k - x_i);

				if (x_i_k_1 - x_i_1 == 0.0)
					t2 = 0.0;
				else
					t2 = x_i_k_1 / (x_i_k_1 - x_i_1);

				Coefficients_New[m][mm][m1] = Coefficients[m][mm][m1] * t1 + Coefficients[m + 1][mm][m1] * t2;
			}
		}

		for (size_t mm = 0; mm < Knot_Vector.size() - 1; mm++) {

			double x_i = Knot_Vector[m];
			double x_i_1;
			double x_i_k;
			double x_i_k_1;

			x_i_1 = Knot_Vector[m + 1];
			x_i_k = Knot_Vector[m + k];
			x_i_k_1 = Knot_Vector[m + k + 1];

			double t1, t2;
			if (x_i_k - x_i == 0.0)
				t1 = 0.0;
			else
				t1 = 1.0 / (x_i_k - x_i);

			if (x_i_k_1 - x_i_1 == 0.0)
				t2 = 0.0;
			else
				t2 = (-1.0) / (x_i_k_1 - x_i_1);

			for (size_t m1 = 0; m1 < Coefficients[m][mm].size(); m1++) { // Orders 0 to k - 1
				Coefficients_New[m][mm][m1 + 1] += Coefficients[m][mm][m1] * t1 + Coefficients[m + 1][mm][m1] * t2;
			}
		}
	}

	return 1;
}

int PrintSplines(const int& Mode, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	if (Mode == 1){
		for (size_t m = 0; m < Initial_Intervals.size(); m++) {
			if (m != 0) {
				std::cout << "**********************************************\n";
				std::cout << "\n" << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m] - 1 + 1] << ": " << std::endl;
			}
			else {
				std::cout << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m] - 1 + 1] << ": " << std::endl;
			}
			for (size_t m1 = 0; m1 < Coefficients.size(); m1++) {
				int Temp = 0;
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
						Temp++;
					}
				}
				if (Temp > 0) {
					std::cout << "\nSpline: " << m1 + 1 << std::endl;
				}
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (m2 != Coefficients[m1][Initial_Intervals[m] - 1].size() - 1) {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							std::cout << " + " << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
					else {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							std::cout << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
				}
			}
			std::cout << std::endl;
		}
	}
	else if (Mode == 2){
		std::cout << Coefficients.size() << " Splines (Order " << Knot_Vector.size() - Coefficients.size() << ") ";
		std::cout << "for knot vector [";
		for (size_t m = 0; m < Knot_Vector.size(); m++) {
			if (m < Knot_Vector.size() - 1)
				std::cout << Knot_Vector[m] << ",";
			else
				std::cout << Knot_Vector[m];
		}
		std::cout << "]: ";
		std::cout<<std::endl;
		for (size_t m1 = 0; m1 < Coefficients.size(); m1++) { //Number of Splines
			std::cout << "**********************************************\n";
			std::cout << "\nSpline: " << m1 + 1 << std::endl;
			for (size_t m = 0; m < Initial_Intervals.size(); m++) {
				std::cout << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m]] << ": "<<std::endl;
				bool Temp = false;
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (m2 != Coefficients[m1][Initial_Intervals[m] - 1].size() - 1) {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							Temp = true;
							std::cout << " + " << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
					else {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							Temp = true;
							std::cout << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
				}
				if (Temp == false) {
					std::cout.precision(2);
					std::cout << 0.0;
				}
				std::cout << std::endl;
			}
		}
	}

	return 1;
}

int Paraview_Visualization(const std::string& FileName, const int& Intervals, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	system("del *.vtk");
	size_t Num_Splines = Coefficients.size();
	for (size_t m = 0; m < Num_Splines; m++) {
		int NumberofPoints = Intervals + 2;
		std::ofstream FileOutput;

		size_t Location = FileName.find(".");
		std::string NewFileName = FileName.substr(0, Location) + "_" + std::to_string(m + 1) + ".vtk";
		FileOutput.open(NewFileName);

		FileOutput << "# vtk DataFile Version 4.1.0" << "\r\n";
		FileOutput << "XFEM Outputs" << "\r\n";
		FileOutput << "ASCII" << std::endl;
		FileOutput << "DATASET UNSTRUCTURED_GRID" << std::endl;
		FileOutput << std::endl;
		FileOutput << "POINTS" << " " << NumberofPoints << " " << "float" << std::endl;
		FileOutput << std::fixed;
		FileOutput << std::setprecision(8);

		std::vector <double> Point_Temp(3);
		for (int i2 = 0; i2 < NumberofPoints; i2++) {

			std::vector <double> Point(3);

			if (i2 <= NumberofPoints - 2) {
				double t = Knot_Vector.back() / (NumberofPoints - 1) * i2;
				Point[0] = t;

				int Index = 0;
				for (size_t mm = 0; mm < Initial_Intervals.size(); mm++) {
					double Inter_1 = Knot_Vector[Initial_Intervals[mm] - 1];
					double Inter_2 = Knot_Vector[Initial_Intervals[mm] - 1 + 1];

					if (t >= Inter_1 && t < Inter_2) {
						Index = Initial_Intervals[mm];
						break;
					}
				}

				double Spline_Value = 0.0;
				for (size_t mm = 0; mm < Coefficients[m][Index - 1].size(); mm++) {
					Spline_Value += Coefficients[m][Index - 1][mm] * pow(t, mm);
				}

				Point[2] = Spline_Value;

				if (i2 == NumberofPoints - 2) {
					Point_Temp = Point;
				}
			}
			else{
				Point = Point_Temp;
			}

			FileOutput << Point[0] << " " << Point[1] << " " << Point[2] << std::endl;
		}

		FileOutput << "CELLS" << " " << NumberofPoints << " " << 2 * NumberofPoints << std::endl;

		for (int i2 = 0; i2 < NumberofPoints; i2++) {
			FileOutput << 1 << " " << i2 << std::endl;
		}

		FileOutput << "CELL_TYPES" << " " << NumberofPoints << std::endl;
		for (int i2 = 0; i2 < NumberofPoints; i2++) {
			FileOutput << 1 << std::endl;
		}
		FileOutput.close();
	}

	return 1;
}