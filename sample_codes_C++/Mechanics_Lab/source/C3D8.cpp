#include "C3D8.h"

Class_C3D8::Class_C3D8(int Nnode_In, int NumofNaturalCoords_In, int NumberofIP_Normal_In)
{
	Nnode = Nnode_In;
	NumofNaturalCoords = NumofNaturalCoords_In;
	NumberofIP_Normal = NumberofIP_Normal_In;

	JJ.resize(NumofNaturalCoords, std::vector <double>(NumofNaturalCoords));
	InverseJJ.resize(NumofNaturalCoords, std::vector <double>(NumofNaturalCoords));
	NN.resize(Nnode);
	DNDX.resize(Nnode, std::vector <double>(3));
	BB.resize(6, std::vector <double>(3 * Nnode));
	Natural_Coords_IP.resize(NumberofIP_Normal, std::vector <double>(3));
	W0.resize(NumberofIP_Normal);
	DNN_DXi.resize(NumberofIP_Normal, std::vector <std::vector <double> >(Nnode, std::vector<double>(NumofNaturalCoords)));

	double* GaussIntegrationPoints_Temp = new double[6 * NumberofIP_Normal];

	HexahedronGaussIntegrationPoints(NumberofIP_Normal, GaussIntegrationPoints_Temp);
	
	for (auto m = 0; m<NumberofIP_Normal; m++) {
		Natural_Coords_IP[m][0] = GaussIntegrationPoints_Temp[0 * NumberofIP_Normal + m];
		Natural_Coords_IP[m][1] = GaussIntegrationPoints_Temp[1 * NumberofIP_Normal + m];
		Natural_Coords_IP[m][2] = GaussIntegrationPoints_Temp[2 * NumberofIP_Normal + m];
		W0[m] = GaussIntegrationPoints_Temp[3 * NumberofIP_Normal + m] * GaussIntegrationPoints_Temp[4 * NumberofIP_Normal + m] *
			GaussIntegrationPoints_Temp[5 * NumberofIP_Normal + m];
	}
	delete[] GaussIntegrationPoints_Temp;
	
	for (auto i1 = 0; i1 < NumberofIP_Normal; i1++) {
		ShapeFunctions_Natural(Natural_Coords_IP[i1], NN);
		ShapeFunctionsDerivatives_Natural(Natural_Coords_IP[i1], DNN_DXi[i1]);
	}
}

int Class_C3D8::Obtain_C3D8_Local_Normal_Assembly(const std::vector <std::vector <double> >& Coords, const std::vector <std::vector <double> >& CC,
	const std::vector <std::vector <double> >& U, std::vector <double>& Rhs, std::vector <double> Sig, std::vector <std::vector <double> >& Amatrx)
{
	int DOF = 3;
	int Nnode = Coords.size();

	for (auto i1 = 0; i1 < NumberofIP_Normal; i1++) {
		ObtainElementPoint_Jacobian(Coords, DNN_DXi[i1], JJ);
		ObtainElementPoint_VolAndInverseJacobian(JJ, Volume, InverseJJ);
		ObtainElementPoint_ShapeFunctionsDerivatives(DNN_DXi[i1], InverseJJ, DNDX);
		ObtainElementPoint_BMatrix(DNDX, BB);
		ObtainElement_pAmatrx(DOF, Nnode, BB, Volume, NN, W0[i1], CC, Amatrx, U, Rhs, Sig);
		
		//printf("Test NN: %f, %f, %f, %f, %f, %f, %f, %f\n", NN[0], NN[1], NN[2], NN[3], NN[4], NN[5], NN[6], NN[7]);
	}

	return 1;
}