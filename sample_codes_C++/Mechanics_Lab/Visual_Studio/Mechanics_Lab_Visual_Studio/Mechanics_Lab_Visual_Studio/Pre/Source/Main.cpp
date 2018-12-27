#include <iostream>
#include "Utilities.h"
#include "FEM_Object_3D.h"

int main(int argc, char* argv[])
{
	FEM_Object_3D FEM_Instance_3D;
	FEM_Instance_3D.Read_Input_File(argv[1]);
	FEM_Instance_3D.Solver();

	system("pause");
	return 1;
}