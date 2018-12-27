#include "Material.h"

int Class_Material::Obtain_MaterialName_Properties(std::string Material_Name, std::vector <std::vector <double> >& CC)
{
	std::string  Material_Type = Material_Name_Type.find(Material_Name)->second;
	std::vector <double> Material_Properties = Material_Name_Properties.find(Material_Name)->second;

	if (Material_Type == "ELASTIC") {
		double EE, nu;
		EE = Material_Properties[0];
		nu = Material_Properties[1];
		ObtainIsotropicStiffness(EE, nu, CC);
	}

	return 1;
}