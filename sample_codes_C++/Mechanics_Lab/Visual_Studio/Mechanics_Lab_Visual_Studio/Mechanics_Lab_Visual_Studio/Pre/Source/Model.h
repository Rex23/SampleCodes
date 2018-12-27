#pragma once
#include "Utilities.h"
#include "Mesh_3D.h"
#include "Material.h"
#include "Section.h"
#include "Surface.h"
#include "Load.h"
#include "Boundary.h"
#include "Step.h"

class Class_Model : virtual New_Class_Utilities::Class_Utilities
{
public:
	Class_Model(Class_Mesh& Mesh, Class_Material& Material, struct_Section& Section, Class_Surface& Surface,
		Class_Load& Load, Class_Boundary& Boundary, Class_Step& Step);

	virtual ~Class_Model() {};

	int ObtainElementsMaterialName(Class_Mesh& Mesh, Class_Material& Material, struct_Section& Section);

public:
	Class_Mesh& Mesh;
	Class_Material& Material; 
	struct_Section& Section; 
	Class_Surface& Surface;
	Class_Load& Load;
	Class_Boundary& Boundary;
	Class_Step& Step;

	std::vector < std::string > Connectivities_Material_Name;
};
