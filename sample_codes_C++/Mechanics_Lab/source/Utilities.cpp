#include "Utilities.h"

double New_Class_Utilities::Class_Utilities::Length(const double* Point1, const double* Point2)
{
	return (sqrt(pow(Point1[0] - Point2[0], 2) + pow(Point1[1] - Point2[1], 2) + pow(Point1[2] - Point2[2], 2)));
}

double New_Class_Utilities::Class_Utilities::Length(const vector<double>& Point1, const vector<double>& Point2)
{
	double Temp = 0.0E0;

	for (int i1 = 0; i1<Point1.size(); i1++) {
		Temp += pow(Point1[i1] - Point2[i1], 2);
	}

	Temp = sqrt(Temp);

	return (Temp);
}

double New_Class_Utilities::Class_Utilities::ObtainPolygonArea(const std::vector <std::vector <double> >& Points)
{
	if (Points.size() == 4) {
		double Area1, Area2;
		double s1, s2;
		double a, b, c, d, e;

		a = Length(Points[0], Points[1]);
		b = Length(Points[1], Points[2]);
		c = Length(Points[2], Points[3]);
		d = Length(Points[3], Points[0]);
		e = Length(Points[0], Points[2]);

		s1 = (a + b + e) / 2.0E0;
		s2 = (e + c + d) / 2.0E0;
		Area1 = sqrt(fabs(s1*(s1 - a)*(s1 - b)*(s1 - e)));
		Area2 = sqrt(fabs(s2*(s2 - e)*(s2 - c)*(s2 - d)));
		return (Area1 + Area2);
	}
	else if (Points.size() == 3) {
		double Area;
		double s, a, b, c;

		a = Length(Points[0], Points[1]);
		b = Length(Points[1], Points[2]);
		c = Length(Points[2], Points[0]);
		s = (a + b + c) / 2.0E0;
		Area = sqrt(fabs(s*(s - a)*(s - b)*(s - c)));
		return Area;
	}
}

std::vector <double> New_Class_Utilities::Class_Utilities::RotationOfVector(const std::vector <std::vector <double> >& RotationMatrix, const std::vector <double> Vector)
{
	std::vector <double> Rotated_Vector(3);

	for (auto m = 0; m < RotationMatrix.size(); m++) {
		for (auto m1 = 0; m1 < RotationMatrix[m].size(); m1++) {
			Rotated_Vector[m] += RotationMatrix[m][m1] * Vector[m1];
		}
	}

	return Rotated_Vector;
}

double New_Class_Utilities::Class_Utilities::DotProduct(const double* Vector1, const double* Vector2)
{
	return (Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2]);
}

double New_Class_Utilities::Class_Utilities::DotProduct(const vector<double>& Vector1, const vector<double>& Vector2)
{
	return (Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2]);
}

double New_Class_Utilities::Class_Utilities::Dot_Product(const vector<double> Vector1, const vector<double> Vector2)
{
	if (Vector1.size() != Vector2.size()) {
		fprintf(stdout, "Error in Vector Dimensions!\n");
		exit(-46);
	}

	double Temp = 0;
	for (int m1 = 0; m1<Vector1.size(); m1++) {
		Temp += Vector1[m1] * Vector2[m1];
	}
	return Temp;
}

int New_Class_Utilities::Class_Utilities::ObtainLocalToGlobalRotationMatrix(const std::vector<std::vector <double> >& Local_Directions, std::vector <std::vector <double> >& Rotation_Matrix)
{
	Rotation_Matrix.resize(3, std::vector <double>(3));
	std::vector <double> E1(3), E2(3), E3(3), E10(3), E20(3), E30(3);

	E1 = Local_Directions[0];
	E2 = Local_Directions[1];
	E3 = Local_Directions[2];

	E10[0] = 1.0; E20[1] = 1.0E0; E30[2] = 1.0E0;

	Rotation_Matrix[0][0] = DotProduct(E10, E1); Rotation_Matrix[0][1] = DotProduct(E10, E2); Rotation_Matrix[0][2] = DotProduct(E10, E3);
	Rotation_Matrix[1][0] = DotProduct(E20, E1); Rotation_Matrix[1][1] = DotProduct(E20, E2); Rotation_Matrix[1][2] = DotProduct(E20, E3);
	Rotation_Matrix[2][0] = DotProduct(E30, E1); Rotation_Matrix[2][1] = DotProduct(E30, E2); Rotation_Matrix[2][2] = DotProduct(E30, E3);

	return 1;
}

int New_Class_Utilities::Class_Utilities::ObtainGlobalToLocalRotationMatrix(const std::vector<std::vector <double> >& Local_Directions, std::vector <std::vector <double> >& Rotation_Matrix)
{
	std::vector <double> E1(3), E2(3), E3(3), E10(3), E20(3), E30(3);

	E1 = Local_Directions[0];
	E2 = Local_Directions[1];
	E3 = Local_Directions[2];

	E10[0] = 1.0; E20[1] = 1.0E0; E30[2] = 1.0E0;

	Rotation_Matrix[0][0] = DotProduct(E1, E10); Rotation_Matrix[0][1] = DotProduct(E1, E20); Rotation_Matrix[0][2] = DotProduct(E1, E30);
	Rotation_Matrix[1][0] = DotProduct(E2, E10); Rotation_Matrix[1][1] = DotProduct(E2, E20); Rotation_Matrix[1][2] = DotProduct(E2, E30);
	Rotation_Matrix[2][0] = DotProduct(E3, E10); Rotation_Matrix[2][1] = DotProduct(E3, E20); Rotation_Matrix[2][2] = DotProduct(E3, E30);

	return 1;
}

void New_Class_Utilities::Class_Utilities::VectorNormalization(double *A_Vector)
{
	double A_Vector2[3];

	A_Vector2[0] = A_Vector[0];
	A_Vector2[1] = A_Vector[1];
	A_Vector2[2] = A_Vector[2];

	A_Vector[0] = A_Vector2[0] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[1] = A_Vector2[1] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[2] = A_Vector2[2] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
}

int New_Class_Utilities::Class_Utilities::VectorNormalization2(vector<double>& A_Vector)
{
	vector<double> A_Vector2(3);

	A_Vector2[0] = A_Vector[0];
	A_Vector2[1] = A_Vector[1];
	A_Vector2[2] = A_Vector[2];

	A_Vector[0] = A_Vector2[0] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[1] = A_Vector2[1] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[2] = A_Vector2[2] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));

	if (A_Vector[0] != A_Vector[0] || A_Vector[1] != A_Vector[1] || A_Vector[2] != A_Vector[2]) {
		return 0;
	}
	else {
		return 1;
	}
}

int New_Class_Utilities::Class_Utilities::VectorNormalization2(double* A_Vector)
{
	vector<double> A_Vector2(3);

	A_Vector2[0] = A_Vector[0];
	A_Vector2[1] = A_Vector[1];
	A_Vector2[2] = A_Vector[2];

	A_Vector[0] = A_Vector2[0] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[1] = A_Vector2[1] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[2] = A_Vector2[2] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));

	if (A_Vector[0] != A_Vector[0] || A_Vector[1] != A_Vector[1] || A_Vector[2] != A_Vector[2]) {
		return 0;
	}
	else {
		return 1;
	}
}

void New_Class_Utilities::Class_Utilities::VectorNormalization(vector<double>& A_Vector)
{
	vector<double> A_Vector2(3);

	A_Vector2[0] = A_Vector[0];
	A_Vector2[1] = A_Vector[1];
	A_Vector2[2] = A_Vector[2];

	A_Vector[0] = A_Vector2[0] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[1] = A_Vector2[1] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
	A_Vector[2] = A_Vector2[2] / sqrt(pow(A_Vector2[0], 2) + pow(A_Vector2[1], 2) + pow(A_Vector2[2], 2));
}

void New_Class_Utilities::Class_Utilities::CrossProduct(const vector<double>& First_Vector, const vector<double>& Second_Vector, vector<double>& CrossProduct)
{
	CrossProduct[0] = First_Vector[1] * Second_Vector[2] - First_Vector[2] * Second_Vector[1];
	CrossProduct[1] = First_Vector[2] * Second_Vector[0] - First_Vector[0] * Second_Vector[2];
	CrossProduct[2] = First_Vector[0] * Second_Vector[1] - First_Vector[1] * Second_Vector[0];
}

void New_Class_Utilities::Class_Utilities::CrossProduct(const double *First_Vector, const double *Second_Vector, double* CrossProduct) {
	CrossProduct[0] = First_Vector[1] * Second_Vector[2] - First_Vector[2] * Second_Vector[1];
	CrossProduct[1] = First_Vector[2] * Second_Vector[0] - First_Vector[0] * Second_Vector[2];
	CrossProduct[2] = First_Vector[0] * Second_Vector[1] - First_Vector[1] * Second_Vector[0];
}

int New_Class_Utilities::Class_Utilities::PlaneUnitNormalVector2(const std::vector<double>& First_Node, const std::vector<double>& Second_Node, const std::vector<double>& Third_Node, std::vector<double>& Unit_Normal)
{
	Unit_Normal.resize(3);

	double First_Vector[3], Second_Vector[3];
	double vector_Temp[3];
	double CrossProduct_Temp[3];

	First_Vector[0] = Second_Node[0] - First_Node[0];
	First_Vector[1] = Second_Node[1] - First_Node[1];
	First_Vector[2] = Second_Node[2] - First_Node[2];

	Second_Vector[0] = Third_Node[0] - Second_Node[0];
	Second_Vector[1] = Third_Node[1] - Second_Node[1];
	Second_Vector[2] = Third_Node[2] - Second_Node[2];

	CrossProduct(First_Vector, Second_Vector, CrossProduct_Temp);
	vector_Temp[0] = CrossProduct_Temp[0];
	vector_Temp[1] = CrossProduct_Temp[1];
	vector_Temp[2] = CrossProduct_Temp[2];
	int Return_Temp = VectorNormalization2(vector_Temp);

	Unit_Normal[0] = vector_Temp[0];
	Unit_Normal[1] = vector_Temp[1];
	Unit_Normal[2] = vector_Temp[2];

	return Return_Temp;
}

int New_Class_Utilities::Class_Utilities::ShapeFunctionsDerivatives_Natural(const std::vector <double>& Natural_Coords, std::vector <std::vector <double> >& DNNDXi)
{
	if (DNNDXi.size() == 8) { //Hexahedron Element

		DNNDXi[0][0] = -1.0 / 8.0*(1.0 - Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		DNNDXi[0][1] = -1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[2]);
		DNNDXi[0][2] = -1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[1]);

		DNNDXi[1][0] = 1.0 / 8.0*(1.0 - Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		DNNDXi[1][1] = -1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[2]);
		DNNDXi[1][2] = -1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[1]);

		DNNDXi[2][0] = 1.0 / 8.0*(1.0 + Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		DNNDXi[2][1] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[2]);
		DNNDXi[2][2] = -1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[1]);

		DNNDXi[3][0] = -1.0 / 8.0*(1.0 + Natural_Coords[1])*(1.0 - Natural_Coords[2]);
		DNNDXi[3][1] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[2]);
		DNNDXi[3][2] = -1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[1]);

		DNNDXi[4][0] = -1.0 / 8.0*(1.0 - Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		DNNDXi[4][1] = -1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[2]);
		DNNDXi[4][2] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 - Natural_Coords[1]);

		DNNDXi[5][0] = 1.0 / 8.0*(1.0 - Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		DNNDXi[5][1] = -1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[2]);
		DNNDXi[5][2] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 - Natural_Coords[1]);

		DNNDXi[6][0] = 1.0 / 8.0*(1.0 + Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		DNNDXi[6][1] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[2]);
		DNNDXi[6][2] = 1.0 / 8.0*(1.0 + Natural_Coords[0])*(1.0 + Natural_Coords[1]);

		DNNDXi[7][0] = -1.0 / 8.0*(1.0 + Natural_Coords[1])*(1.0 + Natural_Coords[2]);
		DNNDXi[7][1] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[2]);
		DNNDXi[7][2] = 1.0 / 8.0*(1.0 - Natural_Coords[0])*(1.0 + Natural_Coords[1]);
	}
	else if (DNNDXi.size() == 10) { //Quadratic Tetrahedron Element

		DNNDXi[0][0] = 4.0E0*Natural_Coords[0] - 1.0E0;
		DNNDXi[0][1] = 0.0E0;
		DNNDXi[0][2] = 0.0E0;
		DNNDXi[0][3] = 0.0E0;

		DNNDXi[1][0] = 0.0E0;
		DNNDXi[1][1] = 4.0E0*Natural_Coords[1] - 1.0E0;
		DNNDXi[1][2] = 0.0E0;
		DNNDXi[1][3] = 0.0E0;

		DNNDXi[2][0] = 0.0E0;
		DNNDXi[2][1] = 0.0E0;
		DNNDXi[2][2] = 4.0E0*Natural_Coords[2] - 1.0E0;
		DNNDXi[2][3] = 0.0E0;

		DNNDXi[3][0] = 0.0E0;
		DNNDXi[3][1] = 0.0E0;
		DNNDXi[3][2] = 0.0E0;
		DNNDXi[3][3] = 4.0E0*Natural_Coords[3] - 1.0E0;

		DNNDXi[4][0] = 4.0E0*Natural_Coords[1];
		DNNDXi[4][1] = 4.0E0*Natural_Coords[0];
		DNNDXi[4][2] = 0.0E0;
		DNNDXi[4][3] = 0.0E0;

		DNNDXi[5][0] = 0.0E0;
		DNNDXi[5][1] = 4.0E0*Natural_Coords[2];
		DNNDXi[5][2] = 4.0E0*Natural_Coords[1];
		DNNDXi[5][3] = 0.0E0;

		DNNDXi[6][0] = 4.0E0*Natural_Coords[2];
		DNNDXi[6][1] = 0.0E0;
		DNNDXi[6][2] = 4.0E0*Natural_Coords[0];
		DNNDXi[6][3] = 0.0E0;

		DNNDXi[7][0] = 4.0E0*Natural_Coords[3];
		DNNDXi[7][1] = 0.0E0;
		DNNDXi[7][2] = 0.0E0;
		DNNDXi[7][3] = 4.0E0*Natural_Coords[0];

		DNNDXi[8][0] = 0.0E0;
		DNNDXi[8][1] = 4.0E0*Natural_Coords[3];
		DNNDXi[8][2] = 0.0E0;
		DNNDXi[8][3] = 4.0E0*Natural_Coords[1];

		DNNDXi[9][0] = 0.0E0;
		DNNDXi[9][1] = 0.0E0;
		DNNDXi[9][2] = 4.0E0*Natural_Coords[3];
		DNNDXi[9][3] = 4.0E0*Natural_Coords[2];
	}
	else if (DNNDXi.size() == 4) { //Linear Tetrahedron Element

		DNNDXi[0][0] = 1.0E0;
		DNNDXi[0][1] = 0.0E0;
		DNNDXi[0][2] = 0.0E0;

		DNNDXi[1][0] = 0.0E0;
		DNNDXi[1][1] = 1.0E0;
		DNNDXi[1][2] = 0.0E0;

		DNNDXi[2][0] = 0.0E0;
		DNNDXi[2][1] = 0.0E0;
		DNNDXi[2][2] = 1.0E0;

		DNNDXi[3][0] = -1.0E0;
		DNNDXi[3][1] = -1.0E0;
		DNNDXi[3][2] = -1.0E0;
	}

	return 1;
}

int New_Class_Utilities::Class_Utilities::ObtainElement_pAmatrx(const int& DOF, const int& Nnode, const std::vector <std::vector <double> >& BB, 
	const double& Volume, const std::vector <double>& NN, const double& W0, const std::vector <std::vector <double> >& CC, 
	std::vector <std::vector <double> >& Amatrx, const std::vector <std::vector <double> >& U, std::vector <double>& Rhs, std::vector <double>& Sig)
{
	//for (int mm = 0; mm < BB.size(); mm++) {
	//	for (int mm1 = 0; mm1 < BB[mm].size(); mm1++) {
	//		printf("Test BB: %d, %d, %f\n", mm + 1, mm1 + 1, BB[mm][mm1]);
	//	}
	//}

	//int Pause; cin >> Pause;

	double Temp;
	int i1_Pseudo, i4_Pseudo;
	div_t divresult;

	vector <vector <double> >* BB_Compress = new vector <vector <double> >(BB.size());
	vector <vector <unsigned int> >* BB_Compress_Index = new vector <vector <unsigned int> >(BB.size());
	for (unsigned int i1 = 0; i1<BB.size(); i1++) {
		for (unsigned int i2 = 0; i2<BB[i1].size(); i2++) {
			if (BB[i1][i2] != 0.0E0) {
				(*BB_Compress)[i1].push_back(BB[i1][i2]);
				(*BB_Compress_Index)[i1].push_back(i2);
			}
		}
	}

	for (unsigned int i2 = 0; i2<(*BB_Compress).size(); i2++) {
		for (unsigned int i3 = 0; i3<(*BB_Compress).size(); i3++) {
			if (CC[i2][i3] != 0.0E0) {
				for (unsigned int i1 = 0; i1<(*BB_Compress)[i2].size(); i1++) {
					for (unsigned int i4 = 0; i4<(*BB_Compress)[i3].size(); i4++) {
						if ((*BB_Compress_Index)[i2][i1] <= (*BB_Compress_Index)[i3][i4]) {
							Temp = (*BB_Compress)[i2][i1] * CC[i2][i3] * (*BB_Compress)[i3][i4] * Volume*W0;
							//printf("Test a few: %f, %f, %f, %f\n", (*BB_Compress)[i2][i1], CC[i2][i3], Volume, W0);
							divresult = div((*BB_Compress_Index)[i2][i1], 3);
							i1_Pseudo = divresult.quot*DOF + divresult.rem;
							divresult = div((*BB_Compress_Index)[i3][i4], 3);
							i4_Pseudo = divresult.quot*DOF + divresult.rem;
							Amatrx[i1_Pseudo][i4_Pseudo] += Temp;
						}
					}
				}
			}
		}
	}

	//PrintScreen_Vector2("Test Amatrx", Amatrx, 1);

	delete BB_Compress;
	delete BB_Compress_Index;

	//for (unsigned int i3 = 0; i3<6; i3++) {
	//	for (unsigned int i4 = 0; i4<6; i4++) {
	//		for (unsigned int i2 = 0; i2<BB[0].size(); i2++) {
	//			int ii1 = i2 / DOF;
	//			int ii2 = i2 % DOF;

	//			if (CC[i3][i4] == 0.0E0 || BB[i4][i2] == 0.0E0 || U[ii1][ii2] == 0.0E0) {
	//			}
	//			else {
	//				Sig[i3] += CC[i3][i4] * BB[i4][i2] * U[ii1][ii2];
	//			}
	//		}
	//	}
	//}

	for (int i1 = 0; i1<BB[0].size(); i1++) {
		for (int i2 = 0; i2<6; i2++) {
			if (BB[i2][i1] == 0.0E0 || Sig[i2] == 0.0E0) {
			}
			else {
				divresult = div(i1, 3);
				i1_Pseudo = divresult.quot*DOF + divresult.rem;
				Temp = -BB[i2][i1] * Sig[i2] * Volume*W0;
				Rhs[i1_Pseudo] += Temp;
			}
		}
	}

	return 1;
}

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
void New_Class_Utilities::Class_Utilities::Vector2Initialization(vector< vector<ValueType> >* Vec, const ValueType& Initial_Value)
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