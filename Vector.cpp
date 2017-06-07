#include "Vector.h"
#include "Matrix.h"
#include <cmath>
#include <string>

Vector :: Vector(){
	dim = CONST_DIM;
	vect = new double[dim];
	for (int i = 0; i < dim; i++) vect[i] = 0;
}

Vector :: Vector(int dim) {
	if (dim <= 0) throw 1;
	this->dim = dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++) vect[i] = 0;
}

Vector::Vector(double *other_vect, int other_dim) {
	if (other_dim <= 0) throw 1;
	dim = other_dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++) vect[i] = other_vect[i];
}

Vector::Vector(double def_val, int dim) {
	if (dim <= 0) throw 1;
	this->dim = dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++) vect[i] = def_val;
}

Vector::Vector(const Matrix &M) {
	if (M.row != 1 && M.col != 1) throw 1;
	else if (M.row == 1) {
		dim = M.col;
		vect = new double[M.col];
		for (int i = 0; i < M.col; i++)
			vect[i] = M.M[0][i];
	}
	else {
		dim = M.row;
		vect = new double[M.row];
		for (int i = 0; i < M.row; i++)
			vect[i] = M.M[i][0];
	}
}

Vector::Vector(const Vector &other) {
	dim = other.dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = other.vect[i];
}

bool Vector::operator==(const Vector &other) {
	if (dim != other.dim) return false;
	for (int i = 0; i < dim; i++) {
		if (fabs(vect[i] - other.vect[i]) >= EPS) return false;
	}
	return true;
}

Vector & Vector::operator=(const Vector &other) {
	delete[] vect;
	dim = other.dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++) vect[i] = other.vect[i];
	return *this;
}

Vector & Vector::operator=(const Matrix &M) {
	delete[] vect;
	if (M.row != 1 && M.col != 1) throw 1;
	else if (M.col == 1) {
		dim = M.row;
		vect = new double[M.row];
		for (int i = 0; i < M.row; i++) vect[i] = M.M[i][0];
	}
	else {
		dim = M.col;
		vect = new double[M.col];
		for (int i = 0; i < M.col; i++) vect[i] = M.M[0][i];
	}
	return *this;
}

Vector operator*(double koeff, const Vector &right) {
	Vector res(right);
	for (int i = 0; i < res.dim; i++) res.vect[i] *= koeff;
	return res;
}

bool Vector::operator!=(const Vector &other) {
	if (dim != other.dim) return true;
	for (int i = 0; i < dim; i++)
		if (fabs(vect[i] - other.vect[i]) >= EPS)
			return true;
	return false;
}

Vector::operator Matrix() { return Matrix(&vect, 1, dim); }

double Vector::operator*(const Vector &other) {
	double comp = 0;
	for (int i = 0; i < dim; i++)
		comp += vect[i] * other.vect[i];
	return comp;
}

Vector operator*(double koeff, const Vector &right) {
	Vector result(right);
	for (int i = 0; i < result.dim; i++)
		result.vect[i] *= koeff;
	return result;
}

Vector Vector::operator*(double koeff) {
	Vector res(*this);
	for (int i = 0; i < dim; i++)
		res.vect[i] *= koeff;
	return res;
}

Vector Vector::operator+(const Vector &other) {
	if (dim != other.dim) throw 1;
	Vector res(dim);
	for (int i = 0; i < dim; i++) res.vect[i] = vect[i] + other.vect[i];
	return res;
}

Vector Vector::operator-(const Vector &other) {
	if (dim != other.dim) throw 1;
	Vector res(dim);
	for (int i = 0; i < dim; i++) res.vect[i] = vect[i] - other.vect[i];
	return res;
}

Vector Vector::operator*(const Matrix &right) {
	if (dim != right.row) throw 1;
	Vector res(right.col);
	for (int ivect = 0; ivect < res.dim; ivect++){
		for (int sum_i = 0; sum_i < right.row; sum_i++) {
			res.vect[ivect] += vect[sum_i] * right.M[sum_i][ivect];
		}
	}
	return res;
}

Vector Vector::operator/(double koeff) {
	if (fabs(koeff) < EPS) throw 1;
	Vector result(*this);
	for (int i = 0; i < dim; i++) result.vect[i] /= koeff;
	return result;
}

double Vector::cube_norm() const {
	double norm = fabs(vect[0]);
	for (int i = 1; i < dim; i++) {
		if (fabs(vect[i]) > norm) norm = fabs(vect[i]);
	}
	return norm;
}

double Vector::octo_norm() const {
	double norm = 0;
	for (int i = 0; i < dim; i++) norm += fabs(vect[i]);
	return norm;
}

double Vector::euclid_norm() const {
	double norm = 0;
	for (int i = 0; i < dim; i++)
		norm += vect[i] * vect[i];
	return sqrt(norm);
}

std::string Vector::to_string() const {
	std::string result = "->";
	int to_ins_i_val;
	double to_ins_d_val;
	for (int i = 0; i < dim - 1; i++) {
		to_ins_i_val = round(vect[i]);
		to_ins_d_val = vect[i];
		result += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + "; ";
	}
	to_ins_i_val = round(vect[dim - 1]);
	to_ins_d_val = vect[dim - 1];
	result += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + "<-";
	return result;
}

int Vector::dimens() const { return dim; }

double Vector::operator[](int index) const {
	if (index < 0 && index >= dim) throw 1;
	else return vect[index];
}

double & Vector::operator[](int index) {
	if (index < 0 && index >= dim) throw 1;
	else return vect[index];
}

Vector::~Vector() {
	delete[] vect;
}

Vector & Vector::operator+=(const Vector &other) {
	if (dim != other.dim) throw 1;
	for (int i = 0; i < dim; i++)
		vect[i] += other.vect[i];
	return *this;
}

Vector & Vector::operator-=(const Vector &other) {
	if (dim != other.dim)
		throw 1;
	for (int i = 0; i < dim; i++)
		vect[i] -= other.vect[i];
	return *this;
}

Vector & Vector::operator*=(double koeff) {
	for (int i = 0; i < dim; i++) vect[i] *= koeff;
	return *this;
}

Vector & Vector::operator/=(double koeff) {
	if (fabs(koeff) < EPS) throw 1;
	for (int i = 0; i < dim; i++) vect[i] /= koeff;
	return *this;
}
