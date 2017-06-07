#pragma once
#include "InterfaceMatrix.h"
#include "Vector.h"

class Matrix; class SquareMatrix;

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

	std::string to_string() const;
	int rows() const;
	int cols() const;

	bool operator==(const Matrix &);
	Matrix& operator=(const Vector &);
	Matrix operator+(const Matrix &);
	Matrix operator-(const Matrix &);
	Matrix operator*(const Matrix &);
	Matrix operator*(double);
	Matrix operator/(double);
	
	double * operator[](int);
	double min_norm() const;

	double cube_norm() const; // кубическая норма max(abs(xi))
	double octo_norm() const; //октоэдрическая норма abs(x1)+abs(x2)+...+abs(xn)
	double euclid_norm() const; // sqrt((x1^2)+x2^2+...+xn^2) евклидова норма(сферическая)
	
	Matrix create_transposed() const;
	void transpose();
	
	operator SquareMatrix();
	operator Vector();
	
	friend class Vector;
	friend Matrix operator*(double, const Matrix&);
	
	virtual ~Matrix();
};

