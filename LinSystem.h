#pragma once
#include "SquareMatrix.h"
#include "Vector.h"
#include "InterfaceLinSystem.h"
class LinSystem : public InterfaceLinSystem {
private:
	int sysDim;
	double err;
	bool is_certain;//is_determined;
	SqrMatrix sys_matrix;
	Vector fill_koef, vdecision; //free_koeffs_vector
public:
	std::string to_string();
	Vector decision() const; //solution() const; �������
	double error() const;
	bool certain() const; //determined() const; ������������
	int dimens() const;

	LinSystem(SqrMatrix, Vector);
	LinSystem(double **, int);
	LinSystem(double **, double *, int);
	
	Vector fill_vector();//free_vector();

	SqrMatrix matrix();
	void matrix_method();

	Matrix extended_matrix();//����������� �������
	Matrix sysswap();//����������� ������� � ���������� ��������
	void optimal();//����� ������������ �������
	//void Cramer_method();

	void Gauss();
	void LU();
	void square_root();
	bool simple_iteration();
	bool Zeidel();
	
	~LinSystem();
};

