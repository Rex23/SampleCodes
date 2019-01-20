#pragma once
#include <map>
#include <vector>
#include <string>

class Class_Boundary
{
public:
	virtual ~Class_Boundary() {};

	std::map < std::string, std::vector <int> > BoundarySet_Directions;
	std::map < std::string, double > BoundarySet_Magnitude;
};
