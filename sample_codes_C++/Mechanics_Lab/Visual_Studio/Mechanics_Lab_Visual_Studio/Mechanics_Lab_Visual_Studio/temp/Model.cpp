#include "Model.h"

Class_Model::Class_Model()
{
	//Global_Stiffness_Row.resize(3 * Mesh.Nodes.size());
	//Global_Stiffness_Column.resize(3 * Mesh.Nodes.size());
}

int Class_Model::ObtainNodesForce()
{
	for (auto iter1 = Connectivities_Pressure.begin(); iter1 != Connectivities_Pressure.end(); iter1++) {

		int Ele_Index = iter1->first;

		std::string Ele_Type = Connectivities_Type[Ele_Index - 1];

		std::map <int, double> Dir_Mag = iter1->second;

		for (auto iter2 = Dir_Mag.begin(); iter2 != Dir_Mag.end(); iter2++) {
			int Dir = iter2->first;
			double Mag = iter2->second;

			std::vector <int> Nodes_Index;
			if (Ele_Type == "C3D8") {
				if (Dir == 1) {
					Nodes_Index = { 1,4,3,2 };
				}
				else if (Dir == 2) {
					Nodes_Index = { 5,6,7,8 };
				}
				else if (Dir == 3) {
					Nodes_Index = { 1,2,6,5 };
				}
				else if (Dir == 4) {
					Nodes_Index = { 2,3,7,6 };
				}
				else if (Dir == 5) {
					Nodes_Index = { 8,7,3,4 };
				}
				else if (Dir == 6) {
					Nodes_Index = { 5,8,4,1 };
				}
			}

			std::vector <std::vector <double> > Local_Directions, Rotation_Matrix;
			Local_Directions = Connectivities_Surface_Unit_Normals[Ele_Index - 1][Dir - 1];

			ObtainLocalToGlobalRotationMatrix(Local_Directions, Rotation_Matrix);

			std::vector <double> Local_Force(3), Global_Force(3);

			std::vector <std::vector <double> > Points;

			for (auto mm = 0; mm < Nodes_Index.size(); mm++) {
				Points.push_back(Nodes[Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1] - 1]);
			}

			double Surf_Area = ObtainPolygonArea(Points);

			Local_Force[2] = Surf_Area * Mag;

			Global_Force = RotationOfVector(Rotation_Matrix, Local_Force);

			//std::cout << Mag << ", " << Surf_Area << ", " << Global_Force[0] << ", " << Global_Force[1] << ", " << Global_Force[2] << std::endl;

			for (auto mm = 0; mm < Nodes_Index.size(); mm++) {
				for (auto mm1 = 0; mm1 < 3; mm1++) {
					if (Global_Force[mm1] != 0.0E0) {

						//////////if (Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1]) != Nodes_Force.end()) {
						//////////	
						//////////	//std::map<int, double> Dir_Force = Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1])->second;
						//////////	
						//////////	if (Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1])->second.find(mm1 + 1) != 
						//////////		Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1])->second.end()) {
						//////////		Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1])->second.find(mm1+1) -> second += Global_Force[mm1] / Nodes_Index.size();
						//////////	}
						//////////	else {
						//////////		Nodes_Force.find(Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1])->second[mm1 + 1] = Global_Force[mm1] / Nodes_Index.size();
						//////////	}
						//////////}
						//////////else {
						//////////	Nodes_Force[Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1]][mm1 + 1] = Global_Force[mm1] / Nodes_Index.size();
						//////////}

						Nodes_Force[Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1]][mm1 + 1] += Global_Force[mm1] / Nodes_Index.size();

						//std::cout << "Node Index: " << Mesh.Connectivities[Ele_Index - 1][Nodes_Index[mm] - 1] << std::endl;
						//std::cout << "Direction: " << mm1 + 1 << std::endl;
						//std::cout << "Global Force: " << Global_Force[mm1] / Nodes_Index.size() << std::endl;
					}
				}
			}
		}

	}

	//for (auto m = Nodes_Force.begin(); m != Nodes_Force.end(); m++)
	//{
	//	for (auto m1 = (m->second).begin(); m1 != (m->second).end(); m1++) {

	//		int Node_Index = m->first;
	//		int Dir = m1->first;
	//		double Mag = m1->second;

	//		std::cout << "N " << Node_Index << ", " << Dir << ", " << Mag << std::endl;

	//	}
	//}

	return 1;
}

int Class_Model::ObtainElementsSurfacePressure()
{
	for (auto iter = DsloadSurface_Magnitude.begin(); iter != DsloadSurface_Magnitude.end(); iter++) {

		std::string Surface_Name = iter->first;
		double Pressure_Magnitude = iter->second;

		//std::map < std::string, std::string > SurfaceName_Type;
		//std::map < std::string, std::map < std::string, std::string > > SurfaceName_Sets;

		std::string Surface_Type = SurfaceName_Type.find(Surface_Name)->second;
		std::map < std::string, std::string > Sets_Direction = SurfaceName_Sets.find(Surface_Name)->second;

		for (auto iter1 = Sets_Direction.begin(); iter1 != Sets_Direction.end(); iter1++) {
			std::string Set_Name = iter1->first;
			std::string Direction = iter1->second;

			if (Surface_Type == "ELEMENT") {

				std::vector <int> Elements = Element_Set.find(Set_Name)->second;

				for (auto m = 0; m < Elements.size(); m++) {
					if (Direction == "S1") {
						Connectivities_Pressure[Elements[m]][1] = Pressure_Magnitude;
					}
					else if (Direction == "S2") {
						Connectivities_Pressure[Elements[m]][2] = Pressure_Magnitude;
					}
					else if (Direction == "S3") {
						Connectivities_Pressure[Elements[m]][3] = Pressure_Magnitude;
					}
					else if (Direction == "S4") {
						Connectivities_Pressure[Elements[m]][4] = Pressure_Magnitude;
					}
					else if (Direction == "S5") {
						Connectivities_Pressure[Elements[m]][5] = Pressure_Magnitude;
					}
					else if (Direction == "S6") {
						Connectivities_Pressure[Elements[m]][6] = Pressure_Magnitude;
					}
				}
				
				//std::cout << "Test Set_Name: " << Set_Name << std::endl;
				//for (auto m = 0; m < Elements.size(); m++) {
				//	std::cout << "Ele_Indexes: " << Elements[m] << std::endl;
				//}
				//int Pause; std::cin >> Pause;
			}
			else { //Has not been implemented yet

			}
		}
	}

	//for (auto iter1 = Connectivities_Pressure.begin(); iter1 != Connectivities_Pressure.end(); iter1++) {
	//	std::map <int, double> temp = iter1 ->second;
	//	for (auto iter2 = temp.begin(); iter2 != temp.end(); iter2++) {
	//		std::cout << "Ele: " << iter1 -> first << std::endl;
	//		std::cout << "Direction: " << iter2->first << ", " << "Magnitude: " << iter2->second << std::endl;
	//	}
	//}
	//system("pause");

	return 1;
}

int Class_Model::ObtainElementsMaterialName()
{
	//std::map <std::string, std::vector <std::string> > SectionSet_Material_Type;

	Connectivities_Material_Name.resize(Connectivities.size());

	for (std::map <std::string, std::vector <std::string> >::const_iterator iter1 = SectionSet_Material_Type.begin();
		iter1 != SectionSet_Material_Type.end(); iter1++) {
		std::string Set_Name = iter1->first;
		std::vector <std::string> Vector_Material_Type = iter1->second;
		std::string Material_Name = Vector_Material_Type[0];
		std::string  Material_Type = Material_Name_Type.find(Material_Name)->second;
		std::vector <double> Material_Properties = Material_Name_Properties.find(Material_Name) -> second;
		std::vector <int> Ele_Set = Element_Set.find(Set_Name)->second;
		for (int m1 = 0; m1 < Ele_Set.size(); m1++) {
			int Ele_Index = Ele_Set[m1];
			Connectivities_Material_Name[Ele_Index - 1] = Material_Name;
		}
	}

	return 1;
}