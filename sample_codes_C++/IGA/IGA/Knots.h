#include "Mesh.h"
#pragma once

class Class_Knots
{
public:
	Class_Knots(Class_Mesh& Mesh_Input) : Mesh(Mesh_Input){};
	int ObtainKnotIntervals_on_edges();

protected:
	Class_Mesh& Mesh;
};