#pragma once
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include "Utilities.h"

class Class_Mesh : virtual public New_Class_Utilities::Class_Utilities
{
	public:
		Class_Mesh() {};
		virtual ~Class_Mesh() {};

		//int Push_Back_Nodes(const std::vector <double>& Coords);
		//int Push_Back_Connectivities(const std::vector <int>& Connec_Single_Line);
		int ObtainConnectivitiesSurfaceUnitNormals();

	public:
		std::vector < std::vector <double> > Nodes;
		std::vector < std::vector <int> > Connectivities;
		std::map < std::string, std::vector <int> > Node_Set;
		std::map < std::string, std::vector <int> > Element_Set;
		std::vector < std::string > Connectivities_Type;

		/*!
		First: Connectivities Index; Second: 1-6 Surfaces; Third: Three unit normals; Fourth: unit normal vector

		C3D8: 
		Surface 1: 1, 4, 3, 2
		Surface 2: 5, 6, 7, 8
		Surface 3: 1, 2, 6, 5
		Surface 4: 2, 3, 7, 6
		Surface 5: 8, 7, 3, 4
		Surface 6: 5, 8, 4, 1
		*/
		std::vector < std::vector <std::vector <std::vector <double> > > > Connectivities_Surface_Unit_Normals;

		bool C3D8_Flag = false;
		bool C3D10_Flag = false;
};