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

	/*!
	For obtaning Connectivities_Pressure
	*/
	int ObtainElementsSurfacePressure();
	
	/*!
	For obtaining Nodes_Force
	*/
	int ObtainNodesForce();

	/*!
	For obtaining Nodes_BC
	*/
	int ObtainNodesAppliedBCs();

public:

	std::vector < std::string > Connectivities_Material_Name;
	std::vector < std::vector <double> > Connectivities_Sig;
	std::vector < std::vector <double> > Nodes_Disp;
	
	/*! 
	First: StepName, Each element has direction and pressure magnitude 
	*/
	std::map<std::string, std::map <int, std::map <int, double> > > Connectivities_Pressure;
	
	/*!
	First: StepName, Each node has direction and pressure magnitude
	*/
	std::map<std::string, std::map <int, std::map <int, double> > > Nodes_Force;

	/*!
	Stiffness matrix for all of the elements
	*/
	std::vector <std::vector <std::vector <double> > > Connectivities_Amatrx;

	/*!
	Rhs for all of the elements
	*/
	std::vector < std::vector <double> > Connectivities_Rhs;

	/*!
	First: StepName, Each node has direction and displacement magnitude
	*/
	std::map <std::string, std::map <int, std::map<int, double> > > Nodes_BC;

	std::vector < std::map <int, double> > Global_Stiffness_Row;
	std::vector < std::map <int, double> > Global_Stiffness_Column;
	std::vector < double > Global_Rhs;
};
