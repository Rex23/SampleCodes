#pragma once
#include "Utilities.h"
#include "C3D8.h"
#include "C3D10.h"

class Class_Model : virtual public Class_C3D8, virtual public Class_C3D10
{
public:

	Class_Model();

	virtual ~Class_Model() {};

	int ObtainElementsMaterialName();
	int ObtainElementsSurfacePressure();
	int ObtainNodesForce();

public:

	std::vector < std::string > Connectivities_Material_Name;
	std::vector < std::vector <double> > Connectivities_Sig;
	std::vector < std::vector <double> > Nodes_Disp;
	
	/*! 
	Each element has direction and pressure magnitude 
	*/
	std::map <int, std::map <int, double> > Connectivities_Pressure;
	
	/*!
	Each node has direction and pressure magnitude
	*/
	std::map <int, std::map <int, double> > Nodes_Force;

	std::vector <std::vector <std::vector <double> > > Connectivities_Amatrx;

	std::vector < std::map <int, double> > Global_Stiffness_Row;
	std::vector < std::map <int, double> > Global_Stiffness_Column;
};
