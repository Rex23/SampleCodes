#pragma once
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include "Utilities.h"
#pragma once

class Class_Mesh : virtual New_Class_Utilities::Class_Utilities
{
	public:
		Class_Mesh() {};
		virtual ~Class_Mesh() {};

		//int Push_Back_Nodes(const std::vector <double>& Coords);
		//int Push_Back_Connectivities(const std::vector <int>& Connec_Single_Line);

	public:
		std::vector < std::vector <double> > Nodes;
		std::vector < std::vector <int> > Connectivities;
		std::map < std::string, std::vector <int> > Node_Set;
		std::map < std::string, std::vector <int> > Element_Set;
		std::vector < std::string > Connectivities_Type;

		bool C3D8_Flag = false;
		bool C3D10_Flag = false;
};