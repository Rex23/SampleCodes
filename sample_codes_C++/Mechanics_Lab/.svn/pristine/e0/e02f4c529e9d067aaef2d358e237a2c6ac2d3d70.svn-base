#include "Model.h"

Class_Model::Class_Model(Class_Mesh& Mesh, Class_Material& Material, struct_Section& Section, Class_Surface& Surface,
	Class_Load& Load, Class_Boundary& Boundary, Class_Step& Step):Mesh(Mesh), Material(Material), Section(Section),
	Surface(Surface), Load(Load), Boundary(Boundary), Step(Step)
{

}

int Class_Model::ObtainElementsMaterialName(Class_Mesh& Mesh, Class_Material& Material, struct_Section& Section)
{
	//std::map <std::string, std::vector <std::string> > SectionSet_Material_Type;

	Connectivities_Material_Name.resize(Mesh.Connectivities.size());

	for (std::map <std::string, std::vector <std::string> >::const_iterator iter1 = Section.SectionSet_Material_Type.begin();
		iter1 != Section.SectionSet_Material_Type.end(); iter1++) {
		std::string Set_Name = iter1->first;
		std::vector <std::string> Vector_Material_Type = iter1->second;
		std::string Material_Name = Vector_Material_Type[0];
		std::string  Material_Type = Material.Material_Name_Type.find(Material_Name)->second;
		std::vector <double> Material_Properties = Material.Material_Name_Properties.find(Material_Name) -> second;
		std::vector <int> Ele_Set = Mesh.Element_Set.find(Set_Name)->second;
		for (int m1 = 0; m1 < Ele_Set.size(); m1++) {
			int Ele_Index = Ele_Set[m1];
			Connectivities_Material_Name[Ele_Index - 1] = Material_Name;
		}
	}

	return 1;
}