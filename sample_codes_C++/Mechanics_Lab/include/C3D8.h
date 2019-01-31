#pragma once
#include <vector>
#include "stdio.h"
#include "Mesh_3D.h"
#include "Material.h"
#include "Section.h"
#include "Surface.h"
#include "Load.h"
#include "Boundary.h"
#include "Step.h"
#include "Utilities.h"

class Class_C3D8 : virtual public Class_Mesh, virtual public Class_Material, virtual public struct_Section, 
				   virtual public Class_Surface, virtual public Class_Load, virtual public Class_Boundary, virtual public Class_Step
{
public:
	Class_C3D8(int Nnode_In = 8, int NumofNaturalCoords_In = 3, int NumberofIP_Normal_In = 8);
	virtual ~Class_C3D8() {}

	int Obtain_C3D8_Local_Normal_Assembly(const std::vector <std::vector <double> >& Coords, const std::vector <std::vector <double> >& CC,
		const std::vector <std::vector <double> >& U, std::vector <double>& Rhs, std::vector <double> Sig, std::vector <std::vector <double> >& Amatrx);

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
