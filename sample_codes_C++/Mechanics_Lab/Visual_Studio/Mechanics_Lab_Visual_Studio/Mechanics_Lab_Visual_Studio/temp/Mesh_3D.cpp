#include "Mesh_3D.h"

int Class_Mesh::ObtainConnectivitiesSurfaceUnitNormals()
{
	//Connectivities_Surface_Unit_Normals.resize(Connectivities.size(), std::vector <std::vector <std::vector <double> > >(6, std::vector <std::vector <double> >(3, std::vector <double>(3))));

	Connectivities_Surface_Unit_Normals.resize(Connectivities.size());

	//std::vector < std::string > Connectivities_Type;

	std::vector <double> N1(3), N2(3), N3(3), N4(3);
	std::vector <std::vector <double> > UnitDirections;

	for (auto m = 0; m < Connectivities_Surface_Unit_Normals.size(); m++) {

		std::string Ele_Type = Connectivities_Type[m];

		if (Ele_Type == "C3D8") {
			Connectivities_Surface_Unit_Normals[m].resize(6, std::vector <std::vector <double> >(3, std::vector <double>(3)));
		
			for (auto m1 = 0; m1 < Connectivities_Surface_Unit_Normals[m].size(); m1++) { //Six Surfaces

				if (m1 == 0) {
					N1 = Nodes[Connectivities[m][1 - 1] - 1];
					N2 = Nodes[Connectivities[m][4 - 1] - 1];
					N3 = Nodes[Connectivities[m][3 - 1] - 1];
					N4 = Nodes[Connectivities[m][2 - 1] - 1];
				}
				else if (m1 == 1) {
					N1 = Nodes[Connectivities[m][5 - 1] - 1];
					N2 = Nodes[Connectivities[m][6 - 1] - 1];
					N3 = Nodes[Connectivities[m][7 - 1] - 1];
					N4 = Nodes[Connectivities[m][8 - 1] - 1];
				}
				else if (m1 == 2) {
					N1 = Nodes[Connectivities[m][1 - 1] - 1];
					N2 = Nodes[Connectivities[m][2 - 1] - 1];
					N3 = Nodes[Connectivities[m][6 - 1] - 1];
					N4 = Nodes[Connectivities[m][5 - 1] - 1];
				}
				else if (m1 == 3) {
					N1 = Nodes[Connectivities[m][2 - 1] - 1];
					N2 = Nodes[Connectivities[m][3 - 1] - 1];
					N3 = Nodes[Connectivities[m][7 - 1] - 1];
					N4 = Nodes[Connectivities[m][6 - 1] - 1];
				}
				else if (m1 == 4) {
					N1 = Nodes[Connectivities[m][8 - 1] - 1];
					N2 = Nodes[Connectivities[m][7 - 1] - 1];
					N3 = Nodes[Connectivities[m][3 - 1] - 1];
					N4 = Nodes[Connectivities[m][4 - 1] - 1];
				}
				else if (m1 == 5) {
					N1 = Nodes[Connectivities[m][5 - 1] - 1];
					N2 = Nodes[Connectivities[m][8 - 1] - 1];
					N3 = Nodes[Connectivities[m][4 - 1] - 1];
					N4 = Nodes[Connectivities[m][1 - 1] - 1];
				}

				std::vector <double> Unit_Normal;

				PlaneUnitNormalVector2(N1, N2, N3, Unit_Normal); //Local Direction 3
				
				std::vector <double> N_Mid_1(3), N_Mid_2(3), N_Mid_Dir(3);
				for (auto nn = 0; nn < 3; nn++) {
					N_Mid_1[nn] = (N1[nn] + N2[nn]) / 2.0;
					N_Mid_2[nn] = (N3[nn] + N4[nn]) / 2.0;
					N_Mid_Dir[nn] = (N_Mid_2[nn] - N_Mid_1[nn]);
				}
				VectorNormalization2(N_Mid_Dir); //Local Direction 1

				std::vector <double> Local_E2(3);

				CrossProduct(Unit_Normal, N_Mid_Dir, Local_E2); //Local Direction 2

				UnitDirections.clear();

				UnitDirections.push_back(N_Mid_Dir);
				UnitDirections.push_back(Local_E2);
				UnitDirections.push_back(Unit_Normal);
				Connectivities_Surface_Unit_Normals[m][m1] = UnitDirections;
			}
		}
	}

	return 1;
}