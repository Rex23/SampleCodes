#include "Mesh.h"
#include "Knots.h"
#include "Utilities.h"
#include <fstream>
#pragma once

class IGA_Object : virtual New_Class_Utilities::Class_Utilities
{
public:
	~IGA_Object();
	int ReadMesh(const std::string& FileName);
	int ProcessMesh();
	int ObtainKnotInformation();

private:
	Class_Mesh Original_Mesh;
	Class_Knots* Knots;
};