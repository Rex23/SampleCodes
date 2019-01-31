#include "FEM_Object_3D.h"

int FEM_Object_3D::Solve_Equation()
{
	return 1;
}

int FEM_Object_3D::Post_Processor_Inc()
{
	return 1;
}

int FEM_Object_3D::Solver()
{
	//int Inc = 10;

	//for (int m = 0; m < Inc; m++){
	
	//double initial_del_t = 

	//for (int m = 0; m < Step.Number_of_Steps; m++) {

	//std::map < std::string, std::vector <double> > StepName_Static_Parameters;

	for (auto iter1 = Model->StepName_Static_Parameters.begin();
		 iter1 != Model->StepName_Static_Parameters.end(); iter1++)
	{
		std::string Step_Name = iter1->first;
		std::vector <double> Step_Parameters = iter1->second;

		double initial_del_t = Step_Parameters[0];
		double total_t = Step_Parameters[1];
		double min_del_t = Step_Parameters[2];
		double max_del_t = Step_Parameters[3];

		double t = 0.0E0;

		bool break_bool = false;

		int Inc = 0;

		while (true) {
			t += initial_del_t;
			if (t >= total_t) {
				t = total_t;
				break_bool = true;
			}
			else {
				break_bool = false;
			}

			Inc++;
			Obtain_Local_Assemblies(Step_Name, Inc);
			Obtain_Global_Assembly_RHS(Step_Name, Inc);
			ApplyBCs(Step_Name, Inc);
			Solve_Equation();
			Post_Processor_Inc();
			std::cout << "Time: " << t << std::endl;
			if (break_bool == true) {
				break;
			}
		}

	}
	
	//}

	return 1;
}

int FEM_Object_3D::ApplyBCs(const std::string& Step_Name, const int& Inc)
{
	//std::map <std::string, std::map <int, std::map<int, double> > > Nodes_BC;
	//std::vector < std::map <int, double> > Global_Stiffness_Row;
	//std::vector < std::map <int, double> > Global_Stiffness_Column;
	//std::vector < double > Global_Rhs;

	std::map <int, std::map<int, double> > Nodes_BC_Step = Model->Nodes_BC.find(Step_Name)->second;

	for (auto m = Nodes_BC_Step.begin(); m != Nodes_BC_Step.end(); m++) {
		int Node_Index = m->first;
		std::map<int, double> Dir_Mag = m->second;

		for (auto m1 = Dir_Mag.begin(); m1 != Dir_Mag.end(); m1++) {
			int Dir = m1->first;
			double Mag = m1->second;

			//For Rhs
			if (Mag != 0.0E0) {
				for (auto mm = 0; mm < Model->Global_Rhs.size(); mm++) {
					int Row_Index = mm + 1;
					int Column_Index = 3 * (Node_Index - 1) + Dir;

					if (Column_Index >= Row_Index) {
						if (Model->Global_Stiffness_Row[Row_Index].find(Column_Index) != Model->Global_Stiffness_Row[Row_Index].end()) {
							double Stiffness_Term = Model->Global_Stiffness_Row[Row_Index].find(Column_Index)->second;
							Model->Global_Rhs[mm] -= Stiffness_Term * Mag;
						}
					}
					else {
						if (Model->Global_Stiffness_Row[Column_Index].find(Row_Index) != Model->Global_Stiffness_Row[Column_Index].end()) {
							double Stiffness_Term = Model->Global_Stiffness_Row[Column_Index].find(Row_Index)->second;
							Model->Global_Rhs[mm] -= Stiffness_Term * Mag;
						}
					}
				}
			}

			//For Stiffness
			//Changing row:
			int Row_Index = 3 * (Node_Index - 1) + Dir;
			for (auto mm = Model->Global_Stiffness_Row[Row_Index-1].begin(); mm != Model->Global_Stiffness_Row[Row_Index-1].end(); mm++) {
				int Column_Index = mm -> first;
				if (Column_Index != Row_Index) {
					mm->second = 0.0E0;
				}
				else {
					mm->second = 1.0E0;
				}
			}
			//Chaning Column:
			int Column_Index = 3 * (Node_Index - 1) + Dir;
			for (auto mm = Model->Global_Stiffness_Column[Column_Index-1].begin(); mm != Model->Global_Stiffness_Column[Column_Index-1].end(); mm++) {
				int Row_Index = mm->first;
				if (Row_Index != Column_Index) {
					mm->second = 0.0E0;
				}
				else {
					mm->second = 1.0E0;
				}
			}
		}

	}

	std::ofstream stream_stiffness;

	stream_stiffness.open("stiffness_test.out");

	for (auto m = Model->Global_Stiffness_Row.begin(); m != Model->Global_Stiffness_Row.end(); m++) {

		stream_stiffness << "Row: " << m - Model->Global_Stiffness_Row.begin() + 1;

		for (auto m1 = (*m).begin(); m1 != (*m).end(); m1++) {
			stream_stiffness << " Index: " << m1->first << ", Value: " << m1->second << ", ";
		}

		stream_stiffness << "\n" << std::endl;

	}

	stream_stiffness.close();

	return 1;
}

int FEM_Object_3D::Obtain_Local_Assemblies(const std::string& Step_Name, const int& Inc)
{
	for (auto m = 0; m < Model->Connectivities.size(); m++){
		Obtain_Local_Assembly(Step_Name, Inc, m+1, Model->Connectivities[m]);
	}

	return 1;
}

int FEM_Object_3D::Obtain_Global_Assembly_RHS(const std::string& Step_Name, const int& Inc)
{
	Model->Global_Stiffness_Row.clear();
	Model->Global_Stiffness_Column.clear();
	Model->Global_Stiffness_Row.resize(3 * Model->Nodes.size());
	Model->Global_Stiffness_Column.resize(3 * Model->Nodes.size());

	for (auto m = 0; m < Model -> Connectivities_Amatrx.size(); m++) {
		std::vector <std::vector <double> > Amatrx = Model -> Connectivities_Amatrx[m];
		for (auto n1 = 0; n1 < Amatrx.size(); n1++) {
			for (auto n2 = 0; n2 < Amatrx[n1].size(); n2++) {
				int Local_Index1 = n1 / 3 + 1;
				int Local_Index2 = n2 / 3 + 1;
				int Global_Index1 = Model->Connectivities[m][Local_Index1 - 1];
				int Global_Index2 = Model->Connectivities[m][Local_Index2 - 1];

				if (Global_Index2 >= Global_Index1) {
					double Value = Amatrx[n1][n2];
					if (Value != 0.0E0) {
						int Dir1 = n1 % 3 + 1;
						int Dir2 = n2 % 3 + 1;
						Model->Global_Stiffness_Row[3 * (Global_Index1 - 1) + (Dir1 - 1)][3 * (Global_Index2 - 1) + (Dir1)] += Value;
						Model->Global_Stiffness_Column[3 * (Global_Index2 - 1) + (Dir1 - 1)][3 * (Global_Index1 - 1) + (Dir1)] += Value;
					}
				}
			}
		}
	}

	Model->Global_Rhs.clear();
	Model->Global_Rhs.resize(3 * Model->Nodes.size());

	for (auto m = 0; m < Model->Connectivities_Rhs.size(); m++) {
		std::vector <double> Rhs = Model->Connectivities_Rhs[m];
		for (auto n1 = 0; n1 < Rhs.size(); n1++) {
			int Local_Index1 = n1 / 3 + 1;
			int Global_Index1 = Model->Connectivities[m][Local_Index1 - 1];
			double Value = Rhs[n1];
			if (Value != 0.0E0) {
				int Dir1 = n1 % 3 + 1;
				Model->Global_Rhs[3 * (Global_Index1 - 1) + (Dir1 - 1)] += Value;
			}
		}
	}

	return 1;
}

int FEM_Object_3D::Obtain_Local_Assembly(const std::string& Step_Name, const int& Inc, const int& Ele_Index, std::vector <int> Connectivity_Local)
{
	struct_Local_Para Local_Para;

	C3D8 = Model;

	Local_Para.Ele_Type = Model->Connectivities_Type[Ele_Index - 1];

	Local_Para.Step_Name = Step_Name;
	Local_Para.Kinc = Inc;
	Local_Para.Jelem = Ele_Index;
	Local_Para.Nnode = Connectivity_Local.size();
	Local_Para.Sig = Model->Connectivities_Sig[Ele_Index - 1];

	Local_Para.Coords.resize(Local_Para.Nnode, std::vector <double>(3));
	for (auto m1 = 0; m1 < Local_Para.Nnode; m1++) {
		for (auto m2 = 0; m2 < 3; m2++) {
			Local_Para.Coords[m1][m2] = Model->Nodes[Connectivity_Local[m1] - 1][m2];
		}
	}

	Local_Para.U.resize(Local_Para.Nnode);
	for (auto m1 = 0; m1 < Local_Para.Nnode; m1++) {
			Local_Para.U[m1] = Model -> Nodes_Disp[Connectivity_Local[m1] - 1];
	}

	Local_Para.Material_Name = Model->Connectivities_Material_Name[Ele_Index - 1];

	Model->Obtain_MaterialName_Properties(Local_Para.Material_Name, Local_Para.CC);

	//PrintScreen_Vector2("CC", Local_Para.CC, 0);

	Obtain_Local_Normal_Assembly(Local_Para);

	Model->Connectivities_Amatrx[Ele_Index - 1] = Local_Para.Amatrx;
	Model->Connectivities_Rhs[Ele_Index - 1] = Local_Para.Rhs;

	return 1;
}

int FEM_Object_3D::Obtain_Local_Normal_Assembly(struct_Local_Para& Local_Para)
{
	Local_Para.Amatrx.resize(3 * Local_Para.Nnode, std::vector <double>(3 * Local_Para.Nnode, 0.0E0));

	if (Local_Para.Ele_Type == "C3D8") {
		C3D8 -> Obtain_C3D8_Local_Normal_Assembly(Local_Para.Coords, Local_Para.CC, Local_Para.U, Local_Para.Rhs, Local_Para.Sig, Local_Para.Amatrx);
	}

	//const std::vector <std::vector <double> >& Coords, const std::vector <std::vector <double> >& CC,
	//const std::vector <double> U, std::vector <double>& Rhs, std::vector <double> Sig, std::vector <std::vector <double> >& Amatrx
	
	return 1;
}

int FEM_Object_3D::Read_Input_File(const std::string& FileName)
{
	std::string tmp2;
	std::vector< int > row;
	std::vector< double > row2;
	std::vector< std::string > row_string;
	std::vector< std::string > Terms;

	std::ifstream infile(FileName, std::ios_base::binary);

	if (infile.fail())
	{
		std::cout << "Error in opening the input file.\n" << std::endl;
		system("pause");
		exit(-3);
	}

	std::string Step_Name = std::string();
	std::string tmp;

	while (getline(infile, tmp)) {

		if (tmp[0] == '*' && tmp[1] != '*') {
			Terms.clear();
			DivideTerms(tmp, &Terms);

			//**************************Mesh**************************
			if (MatchKeyword(Terms[0].c_str(), "*NODE")) {
				while (true) {
					row2.clear();
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						ReadLine(tmp2, row2);
						std::vector<double> row_Temp(row2.begin() + 1, row2.end());
						Model -> Nodes.push_back(row_Temp);
					}
				}
			}
			else if ( MatchKeyword(Terms[0].c_str(), "*ELEMENT") ) {

				std::string Ele_Type;

				if (Terms.size() >= 2) {
					Ele_Type = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("TYPE=") + 5);
					
					if ( Ele_Type == "C3D8" ) {
						Model->C3D8_Flag = true;
					}
					else if (Ele_Type == "C3D10") {
						Model->C3D10_Flag = true;
					}
				}
				else {
					std::cerr << "Error: element type is undefined for the keyword '*ELEMENT'.\n";
					system("pause");
					exit(-1);
				}

				while (true) {
					row.clear();
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						ReadLine(tmp2, row);
						std::vector<int> row_Temp(row.begin() + 1, row.end());
						Model->Connectivities.push_back(row_Temp);
						Model->Connectivities_Type.push_back(Ele_Type);
					}
				}
			}
			//********************************************************

			//**************************Material**************************
			else if ( MatchKeyword(Terms[0].c_str(), "*MATERIAL") ){

				std::string Material_Name = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("=") + 1);
				
				std::string Material_Type;
				while (true) {
					getline(infile, Material_Type);
					if ((Material_Type[0] == '*' && Material_Type[1] == '*')) {
						continue;
					}
					else if ((Material_Type[0] == '*' && Material_Type[1] != '*')) {
						Material_Type = UpperString(Material_Type.substr(1));
						if (!MatchKeyword(Material_Type.c_str(), "ELASTIC")) {
							std::cerr << "Error: material type is not defined.\n";
							system("pause");
							exit(-1);
						}
						break;
					}
				}

				std::string string_Material_Properties;

				while (true) {
					getline(infile, string_Material_Properties);
					if ((string_Material_Properties[0] == '*' && string_Material_Properties[1] == '*')) {
						continue;
					}
					else if ((string_Material_Properties[0] == '*' && string_Material_Properties[1] != '*')) {
						std::cerr << "Error: material properties are not defined.\n";
						system("pause");
						exit(-1);
					}
					else {
						break;
					}
				}
				std::vector <double> Material_Properties;
				ReadLine(string_Material_Properties, Material_Properties);
				//PrintScreen_Vector("Test Material_Properties", Material_Properties, 0);

				Model->Material_Name_Type[Material_Name] = Material_Type;

				//std::cout << "sfsdfsfdMaterial Name Type: " << Material_Name << ", " << Material_Type << std::endl;

				Model->Material_Name_Properties[Material_Name] = Material_Properties;
			}
			//************************************************************

			//**************************Set**************************
			else if ( MatchKeyword(Terms[0].c_str(), "*ELSET") ) {
				//std::cout << Terms[0] << std::endl;
				std::string elset_name = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("=") + 1);
				std::string gen_keyword;
				bool gen_keyword_bool = false;
				if (Terms.size() >= 3) {
					gen_keyword = Terms[2];
					if ((MatchKeyword(gen_keyword.c_str(), "GENERATE") || MatchKeyword(gen_keyword.c_str(), "GEN"))) {
						gen_keyword_bool = true;
						//std::cout << elset_name << std::endl;
					}
				}

				if (!gen_keyword_bool) {

					std::vector <int> ele_set_temp;

					while (true) {
						std::string tmp2;
						std::streamoff previous = infile.tellg();
						getline(infile, tmp2);

						//std::cout << tmp << std::endl;

						if ((tmp2[0] == '*' && tmp2[1] != '*')) {
							infile.seekg(previous, std::ios::beg);
							//std::cout << "Test tmp2" << tmp2 << std::endl;
							//int Pause; std::cin >> Pause;
							break;
						}
						else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
							continue;
						}
						else {
							ReadLine(tmp2, row);

							//for (int m = 0; m < row.size(); m++) {
							//	std::cout << "Test row: " << row[m] << std::endl;
							//}

							ele_set_temp.insert(ele_set_temp.end(), row.begin(), row.end());
						}
					}
					Model->Element_Set[elset_name] = ele_set_temp;

					//std::cout << "Test elset name: " << elset_name << std::endl;

				}
				else {
					std::vector <int> ele_set_temp;
					while (true) {
						std::string tmp2;
						std::streamoff previous = infile.tellg();
						getline(infile, tmp2);

						if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
							infile.seekg(previous, std::ios::beg);
							break;
						}
						else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
							continue;
						}
						else {
							ReadLine(tmp2, row);
							//ele_set_temp.insert(ele_set_temp.end(), row.begin(), row.end());
							for (auto m = row[0]; m <= row[1]; m += row[2]) {
								ele_set_temp.push_back(m);
							}
						}
					}
					Model->Element_Set[elset_name] = ele_set_temp;
				}

				//PrintScreen_Vector(elset_name, Mesh.Element_Set[elset_name], 1);

			}
			else if ( MatchKeyword(Terms[0].c_str(), "*NSET") ) {

				std::string nset_name = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("=") + 1);
				std::string gen_keyword;
				bool gen_keyword_bool = false;
				if (Terms.size() >= 3) {
					gen_keyword = Terms[2];
					if ( (MatchKeyword(gen_keyword.c_str(), "GENERATE") || MatchKeyword(gen_keyword.c_str(), "GEN")) ) {
						gen_keyword_bool = true;
						//std::cout << elset_name << std::endl;
					}
				}

				if (!gen_keyword_bool) {

					std::vector <int> nset_temp;

					while (true) {
						std::string tmp2;
						std::streamoff previous = infile.tellg();
						getline(infile, tmp2);

						if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
							infile.seekg(previous, std::ios::beg);
							break;
						}
						else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
							continue;
						}
						else {
							ReadLine(tmp2, row);
							nset_temp.insert(nset_temp.end(), row.begin(), row.end());
						}
					}
					Model->Node_Set[nset_name] = nset_temp;
				}
				else {
					std::vector <int> nset_temp;
					while (true) {
						std::string tmp2;
						std::streamoff previous = infile.tellg();
						getline(infile, tmp2);

						if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
							infile.seekg(previous, std::ios::beg);
							break;
						}
						else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
							continue;
						}
						else {
							ReadLine(tmp2, row);
							for (auto m = row[0]; m <= row[1]; m += row[2]) {
								nset_temp.push_back(m);
							}
						}
					}
					Model->Node_Set[nset_name] = nset_temp;
				}
				//PrintScreen_Vector(nset_name, Mesh.Node_Set[nset_name], 1);
			}
			//*******************************************************

			//**************************Section**************************
			else if (MatchKeyword(Terms[0].c_str(), "*SOLIDSECTION")) {
				
				if (Terms.size() < 3) {
					std::cerr << "Error: Section information is not defined accurately.\n";
					system("pause");
					exit(-1);
				}

				std::string elset;
				std::string material;
				if (UpperString(Terms[1]).find("ELSET=") != std::string::npos) {
					elset = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("ELSET=") + 6);
				}
				else {
					std::cerr << "Error: No 'elset' defined for a solid section.\n";
					exit(-1);
				}

				if (UpperString(Terms[2]).find("MATERIAL=") != std::string::npos) {
					material = UpperString(Terms[2]).substr(UpperString(Terms[2]).find("MATERIAL=") + 9);
				}
				else
				{
					std::cerr << "Error: No 'material' defined for a solid section.\n";
					exit(-1);
				}

				//std::cout << elset << ", " << material << std::endl;
				std::vector <std::string> Material_Type;
				Material_Type.push_back(material);
				Material_Type.push_back("SolidSection");
				Model->SectionSet_Material_Type[elset] = Material_Type;
			}
			//***********************************************************

			//**************************Surface*****************************
			else if (MatchKeyword(Terms[0].c_str(), "*SURFACE")) {
				if (Terms.size() < 3) {
					std::cerr << "Error: surface information is not defined accurately.\n";
					system("pause");
					exit(-1);
				}

				std::string Surface_Name;

				if (UpperString(Terms[2]).find("NAME=") != std::string::npos) {
					Surface_Name = UpperString(Terms[2]).substr(UpperString(Terms[2]).find("NAME=") + 5);
				}
				else {
					std::cerr << "Error: the surface name is not found\n";
					system("pause");
					exit(-1);
				}
				//std::cout << "Test Surface Name: " << Surface_Name << std::endl;

				std::string Surface_Type;
				if (UpperString(Terms[1]).find("TYPE=") != std::string::npos) {
					Surface_Type = UpperString(Terms[1]).substr(UpperString(Terms[1]).find("TYPE=") + 5);
				}
				else {
					std::cerr << "Error: the surface type is not found.\n";
					system("pause");
					exit(-1);
				}

				Model->SurfaceName_Type[Surface_Name] = Surface_Type;

				while (true) {
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					
					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						DivideTerms(tmp2, &Terms);
						std::map <std::string, std::string> Set_Surface;
						Set_Surface[Terms[0]] = Terms[1];
						Model->SurfaceName_Sets[Surface_Name] = Set_Surface;
					}
				}
			}
			//**************************************************************

			//**************************Step************************
			else if (MatchKeyword(Terms[0].c_str(), "*STEP")) {

				if (UpperString(Terms[1]).find("NAME=") != std::string::npos) {
					Step_Name = UpperString(Terms[2]).substr(UpperString(Terms[1]).find("NAME=") + 5);
					Model->Number_of_Steps++;
				}
				else {
					std::cerr << "Error: the step name is not found\n";
					system("pause");
					exit(-1);
				}
			}
			else if (MatchKeyword(Terms[0].c_str(), "*ENDSTEP")) {
				Step_Name = std::string();
			}
			//******************************************************

			//**************************Load*****************************
			else if (MatchKeyword(Terms[0].c_str(), "*DSLOAD")) {
				while (true) {
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					DivideTerms(tmp2, &Terms);

					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						std::string Surface = Terms[0];
						std::string LoadType = Terms[1];
						double Magnitude = atof(Terms[2].c_str());

						Model->DsloadSurface_LoadType[Surface] = LoadType;
						Model->DsloadSurface_Magnitude[Surface] = Magnitude;

						if (!Step_Name.empty()) {
							Model->StepName_Dsloads[Step_Name].push_back(Surface);
						}
						else {
							std::cerr << "The *DSLOAD " << Surface << " does not belong to any step.\n";
							system("pause");
							exit(-1);
						}
					}
				}
			}
			//***********************************************************

			//**************************Boundary*****************************
			else if (MatchKeyword(Terms[0].c_str(), "*BOUNDARY")) {
				while (true) {
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						DivideTerms(tmp2, &Terms);

						std::string Set_or_Node_Index = Terms[0];
						int Direc1 = atoi(Terms[1].c_str());
						int Direc2 = atoi(Terms[2].c_str());

						double Magnitude = 0.0E0;
						if (Terms.size() >= 4) {
							Magnitude = atof(Terms[3].c_str());
						}

						std::vector <int> Directions;

						for (auto m = Direc1; m <= Direc2; m++) {
							Directions.push_back(m);
						}
						Model->BoundarySet_Directions[Terms[0]] = Directions;
						Model->BoundarySet_Magnitude[Terms[0]] = Magnitude;
					
						if (!Step_Name.empty()) {
							Model->StepName_Boundaries[Step_Name].push_back(Terms[0]);
						}
						else {
							std::cerr << "The *BOUNDARY " << Terms[0] << " does not belong to any step.\n";
							system("pause");
							exit(-1);
						}
					
					}
				}
			}
			//***************************************************************

			//**************************Static**************************
			else if (MatchKeyword(Terms[0].c_str(), "*STATIC")) {

				while (true) {
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if ((tmp2[0] == '*' && tmp2[1] != '*') ) {
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else if ((tmp2[0] == '*' && tmp2[1] == '*')) {
						continue;
					}
					else {
						DivideTerms(tmp2, &Terms);

						if (!Step_Name.empty()) {
							//Model->StepName_Boundaries[Step_Name].push_back(Terms[0]);
							for (auto mm = 0; mm < Terms.size(); mm++) {
								Model->StepName_Static_Parameters[Step_Name].push_back( atof(Terms[mm].c_str()) );
								//std::cout << "Test " << atof(Terms[mm].c_str()) << std::endl;
							}
						}
						else {
							std::cerr << "The *STATIC keyword does not belong to any step.\n";
							system("pause");
							exit(-1);
						}

					}
				}
			}
			//**********************************************************

		} //if (tmp[0] == '*' && tmp[1] != '*'){
	}

	infile.close();

	Model->ObtainConnectivitiesSurfaceUnitNormals();
	Model->ObtainElementsMaterialName();
	Model->Connectivities_Sig.resize(Model->Connectivities.size(), std::vector <double>(6));
	Model->Nodes_Disp.resize(Model->Nodes.size(), std::vector <double>(3));
	Model->ObtainElementsSurfacePressure();
	Model->ObtainNodesForce();
	Model->ObtainNodesAppliedBCs();
	Model->Connectivities_Amatrx.resize(Model->Connectivities.size());
	Model->Connectivities_Rhs.resize(Model->Connectivities.size());

	return 1;
}
