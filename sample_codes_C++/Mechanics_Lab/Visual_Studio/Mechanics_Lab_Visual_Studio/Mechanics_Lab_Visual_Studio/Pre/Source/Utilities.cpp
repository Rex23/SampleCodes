#include "Utilities.h"

int New_Class_Utilities::Class_Utilities::ObtainIsotropicStiffness(const double& EE, const double& nu, std::vector <std::vector <double> >& CC)
{
	CC.resize(6, std::vector <double>(6));

	CC[0][0] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*(1.0 - nu);
	CC[1][1] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*(1.0 - nu);
	CC[2][2] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*(1.0 - nu);
	CC[0][1] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[0][2] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[1][0] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[1][2] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[2][0] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[2][1] = EE / (1.0 + nu) / (1.0 - 2.0*nu)*nu;
	CC[3][3] = EE / (1.0 + nu) / 2.0;
	CC[4][4] = EE / (1.0 + nu) / 2.0;
	CC[5][5] = EE / (1.0 + nu) / 2.0;

	return 1;
}

int New_Class_Utilities::Class_Utilities::ObtainElementPoint_BMatrix(vector <vector <double> >& DNDX, vector <vector <double> >& BB)
{
	if (BB[0].size() == 3 * 8) { //Hexahedron Element
		for (int i1 = 0; i1<DNDX.size(); i1++) {
			BB[1 - 1][i1 * 3 + 1 - 1] = DNDX[i1][0];
			BB[2 - 1][i1 * 3 + 2 - 1] = DNDX[i1][1];
			BB[3 - 1][i1 * 3 + 3 - 1] = DNDX[i1][2];

			BB[4 - 1][i1 * 3 + 1 - 1] = DNDX[i1][1];
			BB[4 - 1][i1 * 3 + 2 - 1] = DNDX[i1][0];

			BB[5 - 1][i1 * 3 + 1 - 1] = DNDX[i1][2];
			BB[5 - 1][i1 * 3 + 3 - 1] = DNDX[i1][0];

			BB[6 - 1][i1 * 3 + 2 - 1] = DNDX[i1][2];
			BB[6 - 1][i1 * 3 + 3 - 1] = DNDX[i1][1];
		}
	}
	else if (BB[0].size() == 3 * 10) { //Tetrahedron Element
		for (int i1 = 0; i1<DNDX.size(); i1++) {
			BB[1 - 1][i1 * 3 + 1 - 1] = DNDX[i1][0];
			BB[2 - 1][i1 * 3 + 2 - 1] = DNDX[i1][1];
			BB[3 - 1][i1 * 3 + 3 - 1] = DNDX[i1][2];

			BB[4 - 1][i1 * 3 + 1 - 1] = DNDX[i1][1];
			BB[4 - 1][i1 * 3 + 2 - 1] = DNDX[i1][0];

			BB[5 - 1][i1 * 3 + 1 - 1] = DNDX[i1][2];
			BB[5 - 1][i1 * 3 + 3 - 1] = DNDX[i1][0];

			BB[6 - 1][i1 * 3 + 2 - 1] = DNDX[i1][2];
			BB[6 - 1][i1 * 3 + 3 - 1] = DNDX[i1][1];
		}
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::ObtainElementPoint_ShapeFunctionsDerivatives(vector <vector <double> >& DNNDXi, vector <vector <double> >& InverseJJ, vector <vector <double> >& DNDX)
{
	Vector2Initialization(&DNDX, 0.0E0);

	if (DNNDXi.size() == 8) { //Hexahedron Element
		for (int i1 = 0; i1<DNNDXi.size(); i1++) {
			for (int mm = 0; mm<InverseJJ.size(); mm++) {
				DNDX[i1][0] = DNDX[i1][0] + InverseJJ[0][mm] * DNNDXi[i1][mm];
				DNDX[i1][1] = DNDX[i1][1] + InverseJJ[1][mm] * DNNDXi[i1][mm];
				DNDX[i1][2] = DNDX[i1][2] + InverseJJ[2][mm] * DNNDXi[i1][mm];
			}
		}
	}
	else if (DNNDXi.size() == 10) { //Tetrahedron Element
		for (int i1 = 0; i1<DNNDXi.size(); i1++) {
			for (int mm = 0; mm<InverseJJ.size(); mm++) {
				DNDX[i1][0] = DNDX[i1][0] + InverseJJ[mm][0] * DNNDXi[i1][mm];
				DNDX[i1][1] = DNDX[i1][1] + InverseJJ[mm][1] * DNNDXi[i1][mm];
				DNDX[i1][2] = DNDX[i1][2] + InverseJJ[mm][2] * DNNDXi[i1][mm];
			}
		}
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::ObtainElementPoint_VolAndInverseJacobian(vector <vector <double> >& JJ, double& Volume, vector <vector <double> >& InverseJJ)
{
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	if (JJ.size() == 3) { //Linear Hexahedron or Tetrahedron Element
		a11 = JJ[0][0]; a12 = JJ[0][1]; a13 = JJ[0][2];
		a21 = JJ[1][0]; a22 = JJ[1][1]; a23 = JJ[1][2];
		a31 = JJ[2][0]; a32 = JJ[2][1]; a33 = JJ[2][2];
		Volume = fabs(-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33);
		InverseJJ[0][0] = (-a23*a32 + a22*a33) / Volume;
		InverseJJ[0][1] = (a13*a32 - a12*a33) / Volume;
		InverseJJ[0][2] = (-a13*a22 + a12*a23) / Volume;
		InverseJJ[1][0] = (a23*a31 - a21*a33) / Volume;
		InverseJJ[1][1] = (-a13*a31 + a11*a33) / Volume;
		InverseJJ[1][2] = (a13*a21 - a11*a23) / Volume;
		InverseJJ[2][0] = (-a22*a31 + a21*a32) / Volume;
		InverseJJ[2][1] = (a12*a31 - a11*a32) / Volume;
		InverseJJ[2][2] = (-a12*a21 + a11*a22) / Volume;
	}
	else if (JJ.size() == 4) { //Tetrahedron Element
		double Jx[4], Jy[4], Jz[4];
		double Jxx[4 * 4], Jyy[4 * 4], Jzz[4 * 4];
		for (int i1 = 0; i1<4; i1++) {
			Jx[i1] = JJ[1][i1]; Jy[i1] = JJ[2][i1]; Jz[i1] = JJ[3][i1];
		}
		for (int i1 = 0; i1<4; i1++) {
			for (int i2 = 0; i2<4; i2++) {
				Jxx[i1 * 4 + i2] = Jx[i1] - Jx[i2];
				Jyy[i1 * 4 + i2] = Jy[i1] - Jy[i2];
				Jzz[i1 * 4 + i2] = Jz[i1] - Jz[i2];
			}
		}
		Volume = Jxx[(2 - 1) * 4 + (1 - 1)] * (Jyy[(2 - 1) * 4 + (3 - 1)] * Jzz[(3 - 1) * 4 + (4 - 1)] - Jyy[(3 - 1) * 4 + (4 - 1)] * Jzz[(2 - 1) * 4 + (3 - 1)])
			+ Jxx[(3 - 1) * 4 + (2 - 1)] * (Jyy[(3 - 1) * 4 + (4 - 1)] * Jzz[(1 - 1) * 4 + (2 - 1)] - Jyy[(1 - 1) * 4 + (2 - 1)] * Jzz[(3 - 1) * 4 + (4 - 1)])
			+ Jxx[(4 - 1) * 4 + (3 - 1)] * (Jyy[(1 - 1) * 4 + (2 - 1)] * Jzz[(2 - 1) * 4 + (3 - 1)] - Jyy[(2 - 1) * 4 + (3 - 1)] * Jzz[(1 - 1) * 4 + (2 - 1)]);

		InverseJJ[0][0] = 1.0 / Volume*(Jyy[(4 - 1) * 4 + (2 - 1)] * Jzz[(3 - 1) * 4 + (2 - 1)] - Jyy[(3 - 1) * 4 + (2 - 1)] * Jzz[(4 - 1) * 4 + (2 - 1)]);
		InverseJJ[0][1] = 1.0 / Volume*(Jxx[(3 - 1) * 4 + (2 - 1)] * Jzz[(4 - 1) * 4 + (2 - 1)] - Jxx[(4 - 1) * 4 + (2 - 1)] * Jzz[(3 - 1) * 4 + (2 - 1)]);
		InverseJJ[0][2] = 1.0 / Volume*(Jxx[(4 - 1) * 4 + (2 - 1)] * Jyy[(3 - 1) * 4 + (2 - 1)] - Jxx[(3 - 1) * 4 + (2 - 1)] * Jyy[(4 - 1) * 4 + (2 - 1)]);

		InverseJJ[1][0] = 1.0 / Volume*(Jyy[(3 - 1) * 4 + (1 - 1)] * Jzz[(4 - 1) * 4 + (3 - 1)] - Jyy[(3 - 1) * 4 + (4 - 1)] * Jzz[(1 - 1) * 4 + (3 - 1)]);
		InverseJJ[1][1] = 1.0 / Volume*(Jxx[(4 - 1) * 4 + (3 - 1)] * Jzz[(3 - 1) * 4 + (1 - 1)] - Jxx[(1 - 1) * 4 + (3 - 1)] * Jzz[(3 - 1) * 4 + (4 - 1)]);
		InverseJJ[1][2] = 1.0 / Volume*(Jxx[(3 - 1) * 4 + (1 - 1)] * Jyy[(4 - 1) * 4 + (3 - 1)] - Jxx[(3 - 1) * 4 + (4 - 1)] * Jyy[(1 - 1) * 4 + (3 - 1)]);

		InverseJJ[2][0] = 1.0 / Volume*(Jyy[(2 - 1) * 4 + (4 - 1)] * Jzz[(1 - 1) * 4 + (4 - 1)] - Jyy[(1 - 1) * 4 + (4 - 1)] * Jzz[(2 - 1) * 4 + (4 - 1)]);
		InverseJJ[2][1] = 1.0 / Volume*(Jxx[(1 - 1) * 4 + (4 - 1)] * Jzz[(2 - 1) * 4 + (4 - 1)] - Jxx[(2 - 1) * 4 + (4 - 1)] * Jzz[(1 - 1) * 4 + (4 - 1)]);
		InverseJJ[2][2] = 1.0 / Volume*(Jxx[(2 - 1) * 4 + (4 - 1)] * Jyy[(1 - 1) * 4 + (4 - 1)] - Jxx[(1 - 1) * 4 + (4 - 1)] * Jyy[(2 - 1) * 4 + (4 - 1)]);

		InverseJJ[3][0] = 1.0 / Volume*(Jyy[(1 - 1) * 4 + (3 - 1)] * Jzz[(2 - 1) * 4 + (1 - 1)] - Jyy[(1 - 1) * 4 + (2 - 1)] * Jzz[(3 - 1) * 4 + (1 - 1)]);
		InverseJJ[3][1] = 1.0 / Volume*(Jxx[(2 - 1) * 4 + (1 - 1)] * Jzz[(1 - 1) * 4 + (3 - 1)] - Jxx[(3 - 1) * 4 + (1 - 1)] * Jzz[(1 - 1) * 4 + (2 - 1)]);
		InverseJJ[3][2] = 1.0 / Volume*(Jxx[(1 - 1) * 4 + (3 - 1)] * Jyy[(2 - 1) * 4 + (1 - 1)] - Jxx[(1 - 1) * 4 + (2 - 1)] * Jyy[(3 - 1) * 4 + (1 - 1)]);
	}

	return 1;
}

template <typename ValueType>
inline void New_Class_Utilities::Class_Utilities::Vector2Initialization(vector< vector<ValueType> >* Vec, const ValueType& Initial_Value)
{
	for (int i1 = 0; i1<(*Vec).size(); i1++) {
		for (int i2 = 0; i2<(*Vec)[i1].size(); i2++) {
			(*Vec)[i1][i2] = Initial_Value;
		}
	}
}

int New_Class_Utilities::Class_Utilities::ObtainElementPoint_Jacobian(const vector <vector <double> >& Coords, vector <vector <double> >& DNNDXi, vector <vector <double> >& JJ)
{
	Vector2Initialization(&JJ, 0.0E0);

	if (Coords.size() == 8 || Coords.size() == 4) { //Linear Hexahedron or Linear Tetrahedron Element
		for (int i1 = 0; i1< Coords.size(); i1++) {
			JJ[0][0] = JJ[0][0] + Coords[i1][0] * DNNDXi[i1][0];
			JJ[0][1] = JJ[0][1] + Coords[i1][1] * DNNDXi[i1][0];
			JJ[0][2] = JJ[0][2] + Coords[i1][2] * DNNDXi[i1][0];

			JJ[1][0] = JJ[1][0] + Coords[i1][0] * DNNDXi[i1][1];
			JJ[1][1] = JJ[1][1] + Coords[i1][1] * DNNDXi[i1][1];
			JJ[1][2] = JJ[1][2] + Coords[i1][2] * DNNDXi[i1][1];

			JJ[2][0] = JJ[2][0] + Coords[i1][0] * DNNDXi[i1][2];
			JJ[2][1] = JJ[2][1] + Coords[i1][1] * DNNDXi[i1][2];
			JJ[2][2] = JJ[2][2] + Coords[i1][2] * DNNDXi[i1][2];
		}
	}
	else if (Coords.size() == 10) {
		JJ[0][0] = 1.0E0;  JJ[0][1] = 1.0E0; JJ[0][2] = 1.0E0; JJ[0][3] = 1.0E0;
		for (int i1 = 0; i1< Coords.size(); i1++) {
			JJ[1][0] = JJ[1][0] + Coords[i1][0] * DNNDXi[i1][0];
			JJ[1][1] = JJ[1][1] + Coords[i1][0] * DNNDXi[i1][1];
			JJ[1][2] = JJ[1][2] + Coords[i1][0] * DNNDXi[i1][2];
			JJ[1][3] = JJ[1][3] + Coords[i1][0] * DNNDXi[i1][3];

			JJ[2][0] = JJ[2][0] + Coords[i1][1] * DNNDXi[i1][0];
			JJ[2][1] = JJ[2][1] + Coords[i1][1] * DNNDXi[i1][1];
			JJ[2][2] = JJ[2][2] + Coords[i1][1] * DNNDXi[i1][2];
			JJ[2][3] = JJ[2][3] + Coords[i1][1] * DNNDXi[i1][3];

			JJ[3][0] = JJ[3][0] + Coords[i1][2] * DNNDXi[i1][0];
			JJ[3][1] = JJ[3][1] + Coords[i1][2] * DNNDXi[i1][1];
			JJ[3][2] = JJ[3][2] + Coords[i1][2] * DNNDXi[i1][2];
			JJ[3][3] = JJ[3][3] + Coords[i1][2] * DNNDXi[i1][3];
		}
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::ShapeFunctions_Natural(vector<double>& Natural_Coords, vector<double>& NN)
{
	if (NN.size() == 8) { //Hexahedron Element
		NN[0] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		NN[1] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		NN[2] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		NN[3] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		NN[4] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		NN[5] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		NN[6] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		NN[7] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[1])*(1.0 + Natural_Coords[2]);
	}
	else if (NN.size() == 10) { //Tetrahedron Element
		NN[0] = Natural_Coords[0] * (2.0*Natural_Coords[0] - 1.0);
		NN[1] = Natural_Coords[1] * (2.0*Natural_Coords[1] - 1.0);
		NN[2] = Natural_Coords[2] * (2.0*Natural_Coords[2] - 1.0);
		NN[3] = Natural_Coords[3] * (2.0*Natural_Coords[3] - 1.0);
		NN[4] = 4.0*Natural_Coords[0] * Natural_Coords[1];
		NN[5] = 4.0*Natural_Coords[1] * Natural_Coords[2];
		NN[6] = 4.0*Natural_Coords[2] * Natural_Coords[0];
		NN[7] = 4.0*Natural_Coords[0] * Natural_Coords[3];
		NN[8] = 4.0*Natural_Coords[1] * Natural_Coords[3];
		NN[9] = 4.0*Natural_Coords[2] * Natural_Coords[3];
	}
	else if (NN.size() == 4) {
		NN[0] = Natural_Coords[0];
		NN[1] = Natural_Coords[1];
		NN[2] = Natural_Coords[2];
		NN[3] = Natural_Coords[3];
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::HexahedronGaussIntegrationPoints(int NumberofIP, double* NaturalCoordinatesAndWeights)
{
	double* Xi1 = new double[NumberofIP];
	double* Eta1 = new double[NumberofIP];
	double* Mu1 = new double[NumberofIP];
	double* w1 = new double[NumberofIP];
	double* w2 = new double[NumberofIP];
	double* w3 = new double[NumberofIP];

	if (NumberofIP == 1) {
		double abscissa[1];
		double weights[1];

		abscissa[0] = 0.0000000000000000;
		weights[0] = 2.00000000000000000;

		int Temp = 0;

		for (int i1 = 0; i1<1; i1++) {
			for (int i2 = 0; i2<1; i2++) {
				for (int i3 = 0; i3<1; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 8) { //2*2*2
		for (int i = 0; i<NumberofIP; i++) {
			w1[i] = 1.0E0; w2[i] = 1.0E0; w3[i] = 1.0E0;
		}
		//Xi1[0] = 1.0E0/sqrt(3.0E0); Xi1[1] = -1.0E0/sqrt(3.0E0); Xi1[2] = -1.0E0/sqrt(3.0E0); Xi1[3] = 1.0E0/sqrt(3.0E0);
		//Xi1[4] = 1.0E0/sqrt(3.0E0); Xi1[5] = -1.0E0/sqrt(3.0E0); Xi1[6] = -1.0E0/sqrt(3.0E0); Xi1[7] = 1.0E0/sqrt(3.0E0);

		//Eta1[0] = 1.0E0/sqrt(3.0E0); Eta1[1] = 1.0E0/sqrt(3.0E0); Eta1[2] = -1.0E0/sqrt(3.0E0); Eta1[3] = -1.0E0/sqrt(3.0E0);
		//Eta1[4] = 1.0E0/sqrt(3.0E0); Eta1[5] = 1.0E0/sqrt(3.0E0); Eta1[6] = -1.0E0/sqrt(3.0E0); Eta1[7] = -1.0E0/sqrt(3.0E0);

		//Mu1[0] = 1.0E0/sqrt(3.0E0); Mu1[1] = 1.0E0/sqrt(3.0E0); Mu1[2] = 1.0E0/sqrt(3.0E0); Mu1[3] = 1.0E0/sqrt(3.0E0);
		//Mu1[4] = -1.0E0/sqrt(3.0E0); Mu1[5] = -1.0E0/sqrt(3.0E0); Mu1[6] = -1.0E0/sqrt(3.0E0); Mu1[7] = -1.0E0/sqrt(3.0E0);

		Xi1[0] = -1.0E0 / sqrt(3.0E0); Eta1[0] = -1.0E0 / sqrt(3.0E0); Mu1[0] = -1.0E0 / sqrt(3.0E0);
		Xi1[1] = 1.0E0 / sqrt(3.0E0); Eta1[1] = -1.0E0 / sqrt(3.0E0); Mu1[1] = -1.0E0 / sqrt(3.0E0);
		Xi1[2] = -1.0E0 / sqrt(3.0E0); Eta1[2] = 1.0E0 / sqrt(3.0E0); Mu1[2] = -1.0E0 / sqrt(3.0E0);
		Xi1[3] = 1.0E0 / sqrt(3.0E0); Eta1[3] = 1.0E0 / sqrt(3.0E0); Mu1[3] = -1.0E0 / sqrt(3.0E0);
		Xi1[4] = -1.0E0 / sqrt(3.0E0); Eta1[4] = -1.0E0 / sqrt(3.0E0); Mu1[4] = 1.0E0 / sqrt(3.0E0);
		Xi1[5] = 1.0E0 / sqrt(3.0E0); Eta1[5] = -1.0E0 / sqrt(3.0E0); Mu1[5] = 1.0E0 / sqrt(3.0E0);
		Xi1[6] = -1.0E0 / sqrt(3.0E0); Eta1[6] = 1.0E0 / sqrt(3.0E0); Mu1[6] = 1.0E0 / sqrt(3.0E0);
		Xi1[7] = 1.0E0 / sqrt(3.0E0); Eta1[7] = 1.0E0 / sqrt(3.0E0); Mu1[7] = 1.0E0 / sqrt(3.0E0);
	}
	else if (NumberofIP == 27) {
		double abscissa[3];
		double weights[3];

		abscissa[0] = 0.0000000000000000;
		abscissa[1] = -0.7745966692414834;
		abscissa[2] = 0.7745966692414834;

		weights[0] = 0.8888888888888888;
		weights[1] = 0.5555555555555556;
		weights[2] = 0.5555555555555556;

		int Temp = 0;

		for (int i1 = 0; i1<3; i1++) {
			for (int i2 = 0; i2<3; i2++) {
				for (int i3 = 0; i3<3; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 64) {
		double abscissa[4];
		double weights[4];

		abscissa[0] = -0.3399810435848563;
		abscissa[1] = 0.3399810435848563;
		abscissa[2] = -0.8611363115940526;
		abscissa[3] = 0.8611363115940526;

		weights[0] = 0.6521451548625461;
		weights[1] = 0.6521451548625461;
		weights[2] = 0.3478548451374538;
		weights[3] = 0.3478548451374538;

		int Temp = 0;

		for (int i1 = 0; i1<4; i1++) {
			for (int i2 = 0; i2<4; i2++) {
				for (int i3 = 0; i3<4; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 125) { //5*5*5
		double abscissa[5];
		double weights[5];

		abscissa[0] = 0.0000000000000000;
		abscissa[1] = -0.5384693101056831;
		abscissa[2] = 0.5384693101056831;
		abscissa[3] = -0.9061798459386640;
		abscissa[4] = 0.9061798459386640;

		weights[0] = 0.5688888888888889;
		weights[1] = 0.4786286704993665;
		weights[2] = 0.4786286704993665;
		weights[3] = 0.2369268850561891;
		weights[4] = 0.2369268850561891;

		int Temp = 0;

		for (int i1 = 0; i1<5; i1++) {
			for (int i2 = 0; i2<5; i2++) {
				for (int i3 = 0; i3<5; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 216) { //6*6*6
		double abscissa[6];
		double weights[6];

		abscissa[0] = 0.6612093864662645;
		abscissa[1] = -0.6612093864662645;
		abscissa[2] = -0.2386191860831969;
		abscissa[3] = 0.2386191860831969;
		abscissa[4] = -0.9324695142031521;
		abscissa[5] = 0.9324695142031521;

		weights[0] = 0.3607615730481386;
		weights[1] = 0.3607615730481386;
		weights[2] = 0.4679139345726910;
		weights[3] = 0.4679139345726910;
		weights[4] = 0.1713244923791704;
		weights[5] = 0.1713244923791704;

		int Temp = 0;

		for (int i1 = 0; i1<6; i1++) {
			for (int i2 = 0; i2<6; i2++) {
				for (int i3 = 0; i3<6; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 343) { //7*7*7
		double abscissa[7];
		double weights[7];

		abscissa[0] = 0.0000000000000000;
		abscissa[1] = 0.4058451513773972;
		abscissa[2] = -0.4058451513773972;
		abscissa[3] = -0.7415311855993945;
		abscissa[4] = 0.7415311855993945;
		abscissa[5] = -0.9491079123427585;
		abscissa[6] = 0.9491079123427585;

		weights[0] = 0.4179591836734694;
		weights[1] = 0.3818300505051189;
		weights[2] = 0.3818300505051189;
		weights[3] = 0.2797053914892766;
		weights[4] = 0.2797053914892766;
		weights[5] = 0.1294849661688697;
		weights[6] = 0.1294849661688697;

		int Temp = 0;

		for (int i1 = 0; i1<7; i1++) {
			for (int i2 = 0; i2<7; i2++) {
				for (int i3 = 0; i3<7; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 512) { //8*8*8
		double abscissa[8];
		double weights[8];

		abscissa[0] = -0.1834346424956498;
		abscissa[1] = 0.1834346424956498;
		abscissa[2] = -0.5255324099163290;
		abscissa[3] = 0.5255324099163290;
		abscissa[4] = -0.7966664774136267;
		abscissa[5] = 0.7966664774136267;
		abscissa[6] = -0.9602898564975363;
		abscissa[7] = 0.9602898564975363;

		weights[0] = 0.3626837833783620;
		weights[1] = 0.3626837833783620;
		weights[2] = 0.3137066458778873;
		weights[3] = 0.3137066458778873;
		weights[4] = 0.2223810344533745;
		weights[5] = 0.2223810344533745;
		weights[6] = 0.1012285362903763;
		weights[7] = 0.1012285362903763;

		int Temp = 0;

		for (int i1 = 0; i1<8; i1++) {
			for (int i2 = 0; i2<8; i2++) {
				for (int i3 = 0; i3<8; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 729) { //9*9*9
		double abscissa[9];
		double weights[9];

		abscissa[0] = 0.0E0;
		abscissa[1] = -0.8360311073266358;
		abscissa[2] = 0.8360311073266358;
		abscissa[3] = -0.9681602395076261;
		abscissa[4] = 0.9681602395076261;
		abscissa[5] = -0.3242534234038089;
		abscissa[6] = 0.3242534234038089;
		abscissa[7] = -0.6133714327005904;
		abscissa[8] = 0.6133714327005904;

		weights[0] = 0.3302393550012598;
		weights[1] = 0.1806481606948574;
		weights[2] = 0.1806481606948574;
		weights[3] = 0.0812743883615744;
		weights[4] = 0.0812743883615744;
		weights[5] = 0.3123470770400029;
		weights[6] = 0.3123470770400029;
		weights[7] = 0.2606106964029354;
		weights[8] = 0.2606106964029354;

		int Temp = 0;

		for (int i1 = 0; i1<9; i1++) {
			for (int i2 = 0; i2<9; i2++) {
				for (int i3 = 0; i3<9; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 1000) { //10*10*10
		double abscissa[10];
		double weights[10];

		abscissa[0] = -0.1488743389816312;
		abscissa[1] = 0.1488743389816312;
		abscissa[2] = -0.4333953941292472;
		abscissa[3] = 0.4333953941292472;
		abscissa[4] = -0.6794095682990244;
		abscissa[5] = 0.6794095682990244;
		abscissa[6] = -0.8650633666889845;
		abscissa[7] = 0.8650633666889845;
		abscissa[8] = -0.9739065285171717;
		abscissa[9] = 0.9739065285171717;

		weights[0] = 0.2955242247147529;
		weights[1] = 0.2955242247147529;
		weights[2] = 0.2692667193099963;
		weights[3] = 0.2692667193099963;
		weights[4] = 0.2190863625159820;
		weights[5] = 0.2190863625159820;
		weights[6] = 0.1494513491505806;
		weights[7] = 0.1494513491505806;
		weights[8] = 0.0666713443086881;
		weights[9] = 0.0666713443086881;

		int Temp = 0;

		for (int i1 = 0; i1<10; i1++) {
			for (int i2 = 0; i2<10; i2++) {
				for (int i3 = 0; i3<10; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 1331) { //11*11*11
		double abscissa[11];
		double weights[11];

		abscissa[0] = 0.0000000000000000;
		abscissa[1] = -0.2695431559523450;
		abscissa[2] = 0.2695431559523450;
		abscissa[3] = -0.5190961292068118;
		abscissa[4] = 0.5190961292068118;
		abscissa[5] = -0.7301520055740494;
		abscissa[6] = 0.7301520055740494;
		abscissa[7] = -0.8870625997680953;
		abscissa[8] = 0.8870625997680953;
		abscissa[9] = -0.9782286581460570;
		abscissa[10] = 0.9782286581460570;

		weights[0] = 0.2729250867779006;
		weights[1] = 0.2628045445102467;
		weights[2] = 0.2628045445102467;
		weights[3] = 0.2331937645919905;
		weights[4] = 0.2331937645919905;
		weights[5] = 0.1862902109277343;
		weights[6] = 0.1862902109277343;
		weights[7] = 0.1255803694649046;
		weights[8] = 0.1255803694649046;
		weights[9] = 0.0556685671161737;
		weights[10] = 0.0556685671161737;

		int Temp = 0;

		for (int i1 = 0; i1<11; i1++) {
			for (int i2 = 0; i2<11; i2++) {
				for (int i3 = 0; i3<11; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 1728) { //12*12*12
		double abscissa[12];
		double weights[12];

		abscissa[0] = -0.1252334085114689;
		abscissa[1] = 0.1252334085114689;
		abscissa[2] = -0.3678314989981802;
		abscissa[3] = 0.3678314989981802;
		abscissa[4] = -0.5873179542866175;
		abscissa[5] = 0.5873179542866175;
		abscissa[6] = -0.7699026741943047;
		abscissa[7] = 0.7699026741943047;
		abscissa[8] = -0.9041172563704749;
		abscissa[9] = 0.9041172563704749;
		abscissa[10] = -0.9815606342467192;
		abscissa[11] = 0.9815606342467192;

		weights[0] = 0.2491470458134028;
		weights[1] = 0.2491470458134028;
		weights[2] = 0.2334925365383548;
		weights[3] = 0.2334925365383548;
		weights[4] = 0.2031674267230659;
		weights[5] = 0.2031674267230659;
		weights[6] = 0.1600783285433462;
		weights[7] = 0.1600783285433462;
		weights[8] = 0.1069393259953184;
		weights[9] = 0.1069393259953184;
		weights[10] = 0.0471753363865118;
		weights[11] = 0.0471753363865118;

		int Temp = 0;

		for (int i1 = 0; i1<12; i1++) {
			for (int i2 = 0; i2<12; i2++) {
				for (int i3 = 0; i3<12; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 2197) { //13*13*13
		double abscissa[13];
		double weights[13];

		abscissa[0] = 0.0000000000000000;
		abscissa[1] = -0.2304583159551348;
		abscissa[2] = 0.2304583159551348;
		abscissa[3] = -0.4484927510364469;
		abscissa[4] = 0.4484927510364469;
		abscissa[5] = -0.6423493394403402;
		abscissa[6] = 0.6423493394403402;
		abscissa[7] = -0.8015780907333099;
		abscissa[8] = 0.8015780907333099;
		abscissa[9] = -0.9175983992229779;
		abscissa[10] = 0.9175983992229779;
		abscissa[11] = -0.9841830547185881;
		abscissa[12] = 0.9841830547185881;

		weights[0] = 0.2325515532308739;
		weights[1] = 0.2262831802628972;
		weights[2] = 0.2262831802628972;
		weights[3] = 0.2078160475368885;
		weights[4] = 0.2078160475368885;
		weights[5] = 0.1781459807619457;
		weights[6] = 0.1781459807619457;
		weights[7] = 0.1388735102197872;
		weights[8] = 0.1388735102197872;
		weights[9] = 0.0921214998377285;
		weights[10] = 0.0921214998377285;
		weights[11] = 0.0404840047653159;
		weights[12] = 0.0404840047653159;

		int Temp = 0;

		for (int i1 = 0; i1<13; i1++) {
			for (int i2 = 0; i2<13; i2++) {
				for (int i3 = 0; i3<13; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}
	else if (NumberofIP == 2744) { //14*14*14
		double abscissa[14];
		double weights[14];

		abscissa[0] = -0.1080549487073437;
		abscissa[1] = 0.1080549487073437;
		abscissa[2] = -0.3191123689278897;
		abscissa[3] = 0.3191123689278897;
		abscissa[4] = -0.5152486363581541;
		abscissa[5] = 0.5152486363581541;
		abscissa[6] = -0.6872929048116855;
		abscissa[7] = 0.6872929048116855;
		abscissa[8] = -0.8272013150697650;
		abscissa[9] = 0.8272013150697650;
		abscissa[10] = -0.9284348836635735;
		abscissa[11] = 0.9284348836635735;
		abscissa[12] = -0.9862838086968123;
		abscissa[13] = 0.9862838086968123;

		weights[0] = 0.2152638534631578;
		weights[1] = 0.2152638534631578;
		weights[2] = 0.2051984637212956;
		weights[3] = 0.2051984637212956;
		weights[4] = 0.1855383974779378;
		weights[5] = 0.1855383974779378;
		weights[6] = 0.1572031671581935;
		weights[7] = 0.1572031671581935;
		weights[8] = 0.1215185706879032;
		weights[9] = 0.1215185706879032;
		weights[10] = 0.0801580871597602;
		weights[11] = 0.0801580871597602;
		weights[12] = 0.0351194603317519;
		weights[13] = 0.0351194603317519;

		int Temp = 0;

		for (int i1 = 0; i1<14; i1++) {
			for (int i2 = 0; i2<14; i2++) {
				for (int i3 = 0; i3<14; i3++) {
					Temp = Temp + 1;
					Xi1[Temp - 1] = abscissa[i1];
					Eta1[Temp - 1] = abscissa[i2];
					Mu1[Temp - 1] = abscissa[i3];
					w1[Temp - 1] = weights[i1];
					w2[Temp - 1] = weights[i2];
					w3[Temp - 1] = weights[i3];
				}
			}
		}
	}

	for (int i = 0; i<NumberofIP; i++) {
		NaturalCoordinatesAndWeights[i] = Xi1[i];
		NaturalCoordinatesAndWeights[i + NumberofIP] = Eta1[i];
		NaturalCoordinatesAndWeights[i + 2 * NumberofIP] = Mu1[i];
		NaturalCoordinatesAndWeights[i + 3 * NumberofIP] = w1[i];
		NaturalCoordinatesAndWeights[i + 4 * NumberofIP] = w2[i];
		NaturalCoordinatesAndWeights[i + 5 * NumberofIP] = w3[i];
	}

	delete[] Xi1;
	delete[] Eta1;
	delete[] Mu1;
	delete[] w1;
	delete[] w2;
	delete[] w3;

	return 1;
}

int New_Class_Utilities::Class_Utilities::Cox_de_Boor_Recursion(const int& k, const std::vector <double>& Knot_Vector, std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	//std::vector <std::vector <std::vector <double> > > Coefficients;
	std::vector <std::vector <std::vector <double> > > Coefficients_New;

	//The initial one:

	//int Num_of_Splines = 0;
	//std::vector <int> Initial_Intervals;
	for (size_t m = 0; m < Knot_Vector.size() - 1; m++) {
		if (Knot_Vector[m] < Knot_Vector[m + 1]) {
			//Num_of_Splines = m + 1;
			Initial_Intervals.push_back(m + 1);
		}
	}

	Coefficients.resize(Knot_Vector.size() - 1); //The number of splines for k = 1
	for (size_t m = 0; m < Coefficients.size(); m++) {
		Coefficients[m].resize(Knot_Vector.size() - 1); //The number of intervals
		for (size_t m1 = 0; m1 < Coefficients[m].size(); m1++) {
			Coefficients[m][m1].resize(1);
			//Coefficients[m][m1][0] = 1.0;
		}
	}

	for (size_t m = 0; m < Initial_Intervals.size(); m++) {
		Coefficients[Initial_Intervals[m] - 1][Initial_Intervals[m] - 1][0] = 1.0;
	}

	for (int m = 0; m < k - 1; m++) {
		Cox_de_Boor_Recursion_Single(m + 1, Coefficients, Knot_Vector, Initial_Intervals, Coefficients_New);
		Coefficients = Coefficients_New;
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::Cox_de_Boor_Recursion_Single(const int& k, const std::vector <std::vector <std::vector <double> > >& Coefficients,
	const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients_New)
{
	Coefficients_New.resize(Coefficients.size() - 1);
	for (size_t m = 0; m < Coefficients_New.size(); m++) {
		Coefficients_New[m].resize(Knot_Vector.size() - 1);
		for (size_t m1 = 0; m1 < Coefficients_New[m].size(); m1++) {
			Coefficients_New[m][m1].resize(Coefficients[m][m1].size() + 1);
		}
	}

	for (size_t m = 0; m < Coefficients.size() - 1; m++) { // (xi <= t < xi+1 )

		for (size_t mm = 0; mm < Knot_Vector.size() - 1; mm++) {

			double x_i = Knot_Vector[m];
			double x_i_1;
			double x_i_k;
			double x_i_k_1;

			x_i_1 = Knot_Vector[m + 1];
			x_i_k = Knot_Vector[m + k];
			x_i_k_1 = Knot_Vector[m + k + 1];

			for (size_t m1 = 0; m1 < Coefficients[m][mm].size(); m1++) { // Orders 0 to k - 1

				double t1, t2;
				if (x_i_k - x_i == 0.0)
					t1 = 0.0;
				else
					t1 = (-x_i) / (x_i_k - x_i);

				if (x_i_k_1 - x_i_1 == 0.0)
					t2 = 0.0;
				else
					t2 = x_i_k_1 / (x_i_k_1 - x_i_1);

				Coefficients_New[m][mm][m1] = Coefficients[m][mm][m1] * t1 + Coefficients[m + 1][mm][m1] * t2;
			}
		}

		for (size_t mm = 0; mm < Knot_Vector.size() - 1; mm++) {

			double x_i = Knot_Vector[m];
			double x_i_1;
			double x_i_k;
			double x_i_k_1;

			x_i_1 = Knot_Vector[m + 1];
			x_i_k = Knot_Vector[m + k];
			x_i_k_1 = Knot_Vector[m + k + 1];

			double t1, t2;
			if (x_i_k - x_i == 0.0)
				t1 = 0.0;
			else
				t1 = 1.0 / (x_i_k - x_i);

			if (x_i_k_1 - x_i_1 == 0.0)
				t2 = 0.0;
			else
				t2 = (-1.0) / (x_i_k_1 - x_i_1);

			for (size_t m1 = 0; m1 < Coefficients[m][mm].size(); m1++) { // Orders 0 to k - 1
				Coefficients_New[m][mm][m1 + 1] += Coefficients[m][mm][m1] * t1 + Coefficients[m + 1][mm][m1] * t2;
			}
		}
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::PrintSplines(const int& Mode, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	if (Mode == 1) {
		for (size_t m = 0; m < Initial_Intervals.size(); m++) {
			if (m != 0) {
				std::cout << "**********************************************\n";
				std::cout << "\n" << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m] - 1 + 1] << ": " << std::endl;
			}
			else {
				std::cout << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m] - 1 + 1] << ": " << std::endl;
			}
			for (size_t m1 = 0; m1 < Coefficients.size(); m1++) {
				int Temp = 0;
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
						Temp++;
					}
				}
				if (Temp > 0) {
					std::cout << "\nSpline: " << m1 + 1 << std::endl;
				}
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (m2 != Coefficients[m1][Initial_Intervals[m] - 1].size() - 1) {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							std::cout << " + " << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
					else {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							std::cout << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
				}
			}
			std::cout << std::endl;
		}
	}
	else if (Mode == 2) {
		std::cout << Coefficients.size() << " Splines (Order " << Knot_Vector.size() - Coefficients.size() << ") ";
		std::cout << "for knot vector [";
		for (size_t m = 0; m < Knot_Vector.size(); m++) {
			if (m < Knot_Vector.size() - 1)
				std::cout << Knot_Vector[m] << ",";
			else
				std::cout << Knot_Vector[m];
		}
		std::cout << "]: ";
		std::cout << std::endl;
		for (size_t m1 = 0; m1 < Coefficients.size(); m1++) { //Number of Splines
			std::cout << "**********************************************\n";
			std::cout << "\nSpline: " << m1 + 1 << std::endl;
			for (size_t m = 0; m < Initial_Intervals.size(); m++) {
				std::cout << Knot_Vector[Initial_Intervals[m] - 1] << "<= t <" << Knot_Vector[Initial_Intervals[m]] << ": " << std::endl;
				bool Temp = false;
				for (int m2 = Coefficients[m1][Initial_Intervals[m] - 1].size() - 1; m2 >= 0; --m2) {
					if (m2 != Coefficients[m1][Initial_Intervals[m] - 1].size() - 1) {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							Temp = true;
							std::cout << " + " << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
					else {
						if (Coefficients[m1][Initial_Intervals[m] - 1][m2] != 0.0) {
							Temp = true;
							std::cout << Coefficients[m1][Initial_Intervals[m] - 1][m2] << "*t^" << m2;
						}
					}
				}
				if (Temp == false) {
					std::cout.precision(2);
					std::cout << 0.0;
				}
				std::cout << std::endl;
			}
		}
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::Paraview_Visualization(const std::string& FileName, const int& Intervals, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients)
{
	system("del *.vtk");
	size_t Num_Splines = Coefficients.size();
	for (size_t m = 0; m < Num_Splines; m++) {
		int NumberofPoints = Intervals + 2;
		std::ofstream FileOutput;

		size_t Location = FileName.find(".");
		std::string NewFileName = FileName.substr(0, Location) + "_" + std::to_string(m + 1) + ".vtk";
		FileOutput.open(NewFileName);

		FileOutput << "# vtk DataFile Version 4.0" << "\r\n";
		FileOutput << "XFEM Outputs" << "\r\n";
		FileOutput << "ASCII" << std::endl;
		FileOutput << "DATASET UNSTRUCTURED_GRID" << std::endl;
		FileOutput << std::endl;
		FileOutput << "POINTS" << " " << NumberofPoints << " " << "float" << std::endl;
		FileOutput << std::fixed;
		FileOutput << std::setprecision(8);

		std::vector <double> Point_Temp(3);
		for (int i2 = 0; i2 < NumberofPoints; i2++) {

			std::vector <double> Point(3);

			if (i2 <= NumberofPoints - 2) {
				double t = Knot_Vector.back() / (NumberofPoints - 1) * i2;
				Point[0] = t;

				int Index = 0;
				for (size_t mm = 0; mm < Initial_Intervals.size(); mm++) {
					double Inter_1 = Knot_Vector[Initial_Intervals[mm] - 1];
					double Inter_2 = Knot_Vector[Initial_Intervals[mm] - 1 + 1];

					if (t >= Inter_1 && t < Inter_2) {
						Index = Initial_Intervals[mm];
						break;
					}
				}

				double Spline_Value = 0.0;
				for (size_t mm = 0; mm < Coefficients[m][Index - 1].size(); mm++) {
					Spline_Value += Coefficients[m][Index - 1][mm] * pow(t, mm);
				}

				Point[2] = Spline_Value;

				if (i2 == NumberofPoints - 2) {
					Point_Temp = Point;
				}
			}
			else {
				Point = Point_Temp;
			}

			FileOutput << Point[0] << " " << Point[1] << " " << Point[2] << std::endl;
		}

		FileOutput << "CELLS" << " " << NumberofPoints << " " << 2 * NumberofPoints << std::endl;

		for (int i2 = 0; i2 < NumberofPoints; i2++) {
			FileOutput << 1 << " " << i2 << std::endl;
		}

		FileOutput << "CELL_TYPES" << " " << NumberofPoints << std::endl;
		for (int i2 = 0; i2 < NumberofPoints; i2++) {
			FileOutput << 1 << std::endl;
		}
		FileOutput.close();
	}

	return 1;
}

char* New_Class_Utilities::Class_Utilities::StringToUpper(char* str)
{
	if (!no(str)) {
		unsigned int i = 0;
		while (str[i]) {
			str[i] = toupper(str[i]);
			i++;
		}
	}
	return str;
}

int New_Class_Utilities::Class_Utilities::no(char* param)
{
	if (int(param) != 0 && strlen(param)>0) {
		return 0;
	}
	return 1;
}

int New_Class_Utilities::Class_Utilities::MatchKeyword(const char* str, const char* key2)
{
	char key1[250];
	strcpy(key1, str);
	int ismatch = 0;

	if (strcmp(str, ","))
	{
		char *KEY = StringToUpper(strtok(key1, ","));
		if (!no(KEY) && strcmp(KEY, key2) == 0)
			ismatch = 1;
	}
	return ismatch;
}

int New_Class_Utilities::Class_Utilities::MatchString(const char* str, const char* key2) {
	char key1[250];
	strcpy(key1, str);
	char* sptr = strstr(StringToUpper(key1), key2);
	if (sptr != NULL)
		return sptr - key1 + 1;
	return 0;
}

std::string New_Class_Utilities::Class_Utilities::UpperString(const std::string& Term)
{
	std::string Upper_Term(Term);

	std::locale loc;
	for (std::string::size_type i = 0; i<Term.length(); ++i)
		Upper_Term[i] = std::toupper(Term[i], loc);

	return Upper_Term;
}

void New_Class_Utilities::Class_Utilities::DivideTerms(std::string tmp, std::vector<std::string>* Terms)
{
	(*Terms).clear();

	vector<string> Terms2;

	Terms2.clear();

	int pre_pos = 0;
	//int RowIndex = 0;
	for (size_t i = 0; i < tmp.size(); i++)
	{
		string tmp_str;
		if (tmp[i] == ',')
		{
			//stringstream r;
			tmp_str = tmp.substr(pre_pos, i - pre_pos);
			//r << tmp_str;
			//int tt;
			//r >> tt;

			Terms2.push_back(tmp_str);

			//if (RowIndex == WidthofInputLines -1) break;
			//++RowIndex;

			pre_pos = i + 1;
		}
		else if (i == tmp.size() - 1)
		{
			//stringstream r;
			tmp_str = tmp.substr(pre_pos, i + 1 - pre_pos);
			//r << tmp_str;
			//double tt;
			//r >> tt;

			if (tmp_str.back() == '\r') {
				
				//tmp_str.erase(tmp_str.size() - 1);
				
				tmp_str = tmp_str.substr(0, tmp_str.size() - 1);
				//std::cout << "Found one!\n";
				//int Pause; std::cin >> Pause;
			}

			Terms2.push_back(tmp_str);

			break;
		}
	}

	//*******************************Get rid of the spaces*******************************
	for (size_t i = 0; i<Terms2.size(); i++) {
		string Temp = "";

		for (size_t j = 0; j<Terms2[i].size(); j++) {
			if ((Terms2[i])[j] != ' ') {
				Temp = Temp + (Terms2[i])[j];
			}
		}
		(*Terms).push_back(Temp);
		//printf("Test Temp, %s\n", Temp.c_str());
	}
	//***********************************************************************************
}