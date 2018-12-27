#pragma once
#include "Utilities.h"
#include <string>
#include <vector>
#include <map>

class Class_Material : virtual New_Class_Utilities::Class_Utilities
{
public:
	virtual ~Class_Material() {};

public:

	int Obtain_MaterialName_Properties(std::string Material_Name, std::vector <std::vector <double> >& CC);

	std::map <std::string, std::string> Material_Name_Type;
	std::map <std::string, std::vector <double> > Material_Name_Properties;
};
