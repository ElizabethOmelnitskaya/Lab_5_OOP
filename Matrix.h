#pragma once
#include "InterfaceMatrix.h"
#include "Vector.h"

class Matrix; class SqrMatrix;

Matrix operator*(double, const Matrix&);

class Matrix : public InterfaceMatrix {
private:
	const int CONST_ROWS = 3, CONST_COLS = 3;
	int row, col;
protected:
	double **M;
public:
	Matrix();
	Matrix(int, int);
	Matrix(const Vector &);
	Matrix(const Matrix &);
	Matrix(double, int, int);
	Matrix(double **, int, int);

	std::string to_string() const;
	int rows() const;
	int cols() const;

	bool operator==(const Matrix &);
	bool operator!=(const Matrix &);

	Matrix& operator=(const Vector &);
	Matrix& operator=(const Matrix &other);
	Matrix operator+(const Matrix &);
	Matrix operator-(const Matrix &);
	Matrix operator*(const Matrix &);
	Matrix operator*(double);
	Matrix operator/(double);

	Vector operator*(const Vector &);
	Matrix & operator+=(const Matrix &);
	Matrix & operator-=(const Matrix &);
	Matrix & operator*=(const Matrix &);
	Matrix & operator*=(double);
	Matrix & operator/=(double);

	
	double * operator[](int);
	double min_norm() const;

	double cube_norm() const; // ���������� ����� max(abs(xi))
	double octo_norm() const; //�������������� ����� abs(x1)+abs(x2)+...+abs(xn)
	double euclid_norm() const; // sqrt((x1^2)+x2^2+...+xn^2) ��������� �����(�����������)
	
	Matrix create_transposed() const;
	void transpose();
	
	operator SqrMatrix();
	operator Vector();
	
	friend class Vector;
	friend Matrix operator*(double, const Matrix&);
	
	virtual ~Matrix();
};

