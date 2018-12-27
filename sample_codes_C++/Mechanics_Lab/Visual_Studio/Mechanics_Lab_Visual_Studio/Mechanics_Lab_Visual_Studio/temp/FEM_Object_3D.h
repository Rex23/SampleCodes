#pragma once
#include "Utilities.h"
#include "Mesh_3D.h"
#include "Material.h"
#include "Section.h"
#include "Surface.h"
#include "Load.h"
#include "Boundary.h"
#include "Step.h"
#include "Model.h"
#include "Local_Parameters.h"
#include <fstream>

class FEM_Object_3D : virtual New_Class_Utilities::Class_Utilities
{
public:
	FEM_Object_3D() { Model = new Class_Model(); };
	virtual ~FEM_Object_3D() { delete Model; };
	int Read_Input_File(const std::string& FileName);
	
	int Solver();
	int Obtain_Local_Assemblies(const int& Inc);
	int Obtain_Global_Assembly(const int& Inc);
	int Obtain_Local_Assembly(const int& Inc, const int& Ele_Index, std::vector <int> Connectivity_Local);
	int Obtain_Local_Normal_Assembly(struct_Local_Para& Local_Para);
	
	int Solve_Equation();
	int Post_Processor_Inc();

private:
	//Class_Mesh Mesh;
	//Class_Material Material;
	//struct_Section Section;
	//Class_Surface Surface;
	//Class_Load Load;
	//Class_Boundary Boundary;
	//Class_Step Step;
	
	Class_Model* Model;
	/*!
	Model is derived from C3D8
	*/
	Class_C3D8* C3D8;
};