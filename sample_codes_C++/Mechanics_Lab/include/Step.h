#pragma once
#include <map>
#include <string>
#include <vector>

class Class_Step
{
public:
	virtual ~Class_Step() {};
	std::map < std::string, std::vector <std::string> > StepName_Boundaries;

	/*!
	Step name and Surface name
	*/
	std::map < std::string, std::vector <std::string> > StepName_Dsloads; 
	
	std::map < std::string, std::vector <double> > StepName_Static_Parameters;
	int Number_of_Steps = 0;
};
