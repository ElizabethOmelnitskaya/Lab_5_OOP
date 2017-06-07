#pragma once

#include "InterfaceVector.h"
#include <iostream>
#include <vector>

class Vector; class Matrix;

Vector operator*(double, const Vector &);

class Vector : public InterfaceVector {
private:
	const int CONST_DIM = 5;
	double *vect;
	int dim;
public:
	Vector();
	Vector(int);
	Vector(double*, int);
	Vector(double, int);
	template<class list>
	Vector(list);
	Vector(const Vector &);
	Vector(const Matrix &);

	bool operator==(const Vector &);
	bool operator!=(const Vector &);

	Vector & operator=(const Vector &);
	Vector & operator=(const Matrix &);
	double operator*(const Vector &);
	Vector operator*(double);
	Vector operator*(const Matrix &);
	Vector operator+(const Vector &);
	Vector operator-(const Vector &);
	Vector operator/(double);

	Vector& operator+=(const Vector &);
	Vector& operator-=(const Vector &);
	Vector& operator*=(double);
	Vector& operator/=(double);

	double cube_norm() const; // кубическая норма max(abs(xi))
	double octo_norm() const; //октоэдрическая норма abs(x1)+abs(x2)+...+abs(xn)
	double euclid_norm() const; // sqrt((x1^2)+x2^2+...+xn^2) евклидова норма(сферическая)

	std::string to_string() const;
	int dimens() const;
	operator Matrix();
	double operator[](int) const;
	double &operator[](int);

	friend class Matrix;
	friend Vector operator*(double, const Vector &);

	~Vector();
};

template<class list>
inline Vector::Vector(list lst) {
	dim = lst.size();
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = lst.at(i);
}
