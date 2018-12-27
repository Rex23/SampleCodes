#pragma once
#include "stdio.h"
#include "Utilities.h"
#include <vector>

class Class_C3D8 : virtual New_Class_Utilities::Class_Utilities
{
public:
	Class_C3D8(int Nnode_In = 8, int NumofNaturalCoords_In = 3, int NumberofIP_Normal_In = 8);
	virtual ~Class_C3D8() {};

	int Obtain_C3D8_Local_Normal_Assembly(const std::vector <std::vector <double> >& Coords);

public:
	std::vector <std::vector <double> > JJ;
	std::vector <std::vector <double> > InverseJJ;
	std::vector <double> NN;
	std::vector <std::vector <double> > DNDX;
	std::vector <std::vector <double> > BB;
	std::vector <std::vector <double> > Natural_Coords_IP;
	std::vector <double> W0;
	std::vector <std::vector <std::vector <double> > > DNN_DXi;

	int Nnode;
	int NumofNaturalCoords;
	int NumberofIP_Normal;

	double Volume;
};
