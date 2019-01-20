#pragma once
#include <vector>
#include "stdio.h"
#include "Mesh_3D.h"
#include "Material.h"
#include "Section.h"
#include "Surface.h"
#include "Load.h"
#include "Boundary.h"
#include "Step.h"
#include "Utilities.h"

class Class_C3D10 : virtual public Class_Mesh, virtual public Class_Material, virtual public struct_Section, 
					virtual public Class_Surface, virtual public Class_Load, virtual public Class_Boundary, virtual public Class_Step
{
public:
	Class_C3D10(int Nnode_In = 10, int NumofNaturalCoords_In = 4, int NumberofIP_Normal_In = 5);
	virtual ~Class_C3D10() {};
};
