#pragma once
#include "Matrix.h"
class SqrMatrix : public Matrix {
private:
	int m_dim;
public:
	SqrMatrix();
	SqrMatrix(int);
	SqrMatrix(double, int);
	SqrMatrix(double **, int);
	SqrMatrix(const Matrix&);
	SqrMatrix(const SqrMatrix&);

	SqrMatrix& operator=(const Matrix&);

	bool symmetry() const;
	bool minors() const;
	bool diagonal_dominating() const;
	int dimens() const;
	SqrMatrix cr_revers() const; // create

	double determ() const;
	~SqrMatrix();
};

