#pragma once
#include <map>
#include <vector>
#include <string>

class Class_Load
{
public:
	virtual ~Class_Load() {};

	std::map < std::string, std::string > DsloadSurface_LoadType;
	std::map < std::string, double > DsloadSurface_Magnitude;
};
