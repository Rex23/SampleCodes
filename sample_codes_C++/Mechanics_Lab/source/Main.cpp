#include <iostream>
#include "Utilities.h"
#include "FEM_Object_3D.h"

using namespace std;

int main(int argc, char* argv[])
{
	std::cout << "Enter main...\n" << std::endl;

	FEM_Object_3D FEM_Instance_3D;
	FEM_Instance_3D.Read_Input_File(argv[1]);
	FEM_Instance_3D.Solver();

	system("pause");
	return 1;
}
