#pragma once
#include <map>
#include <vector>
#include <string>

struct struct_Section
{
public:
	virtual ~struct_Section() {};

	std::map <std::string, std::vector <std::string> > SectionSet_Material_Type;
};
