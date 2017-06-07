#include "Matrix.h"
#include "SquareMatrix.h"
#include <string>

int Matrix::rows() const { return row; }

int Matrix::cols() const { return col; }

Matrix::Matrix(){
	row = CONST_ROWS;
	col = CONST_COLS;
	M = new double*[row];
	for (int i = 0; i < row; i++) M[i] = new double[col];
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){ 
			M[i][j] = 0; 
		}
	}
}

Matrix::Matrix(int n, int m){
	if (n <= 0 || m <= 0) throw 1;
	this->row = n;
	this->col = m;
	M = new double*[n];
	for (int i = 0; i < n; i++) M[i] = new double[m];
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			M[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix &other) {
	row = other.row;
	col = other.col;
	M = new double*[row];
	for (int i = 0; i < row; i++)
		M[i] = new double[col];
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] = other.M[i][j];
}

Matrix::Matrix(const Vector &vector) {
	row = 1;
	col = vector.dim;
	M = new double*[1];
	M[0] = new double[col];
	for (int i = 0; i < col; i++){
		M[0][i] = vector.vect[i];
	}
}

Matrix::Matrix(double def_val, int n, int m) {
	if (n <= 0 || m <= 0) throw 1;
	this->row = n;
	this->col = m;
	M = new double*[n];
	for (int i = 0; i < n; i++) M[i] = new double[m];
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++) {
			M[i][j] = def_val;
		}
	}
}

Matrix::Matrix(double **M, int n, int m) {
	if (n <= 0 || m <= 0) throw 1;
	this->row = n;
	this->col = m;
	this->M = new double*[n];
	for (int i = 0; i < n; i++)
		this->M[i] = new double[m];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			this->M[i][j] = M[i][j];
}

Matrix::Matrix(double default_val, int rows_num, int cols_num) {
	if (rows_num <= 0 || cols_num <= 0) throw 1;
	this->row = rows_num;
	this->col = cols_num;
	M = new double*[rows_num];
	for (int i = 0; i < rows_num; i++) M[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			M[i][j] = default_val;
}

Matrix::Matrix(double **mat, int rows_num, int cols_num) {
	if (rows_num <= 0 || cols_num <= 0) throw 1;
	this->row = rows_num;
	this->col = cols_num;
	this->M = new double*[rows_num];
	for (int i = 0; i < rows_num; i++) this->M[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			this->M[i][j] = mat[i][j];
}

double * Matrix::operator[](int index) {
	if (index < 0 || index >= row) throw 1;
	return M[index];
}

bool Matrix::operator==(const Matrix &other) {
	if (row != other.row || col != other.col) return false;
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			if (fabs(M[i][j] - other.M[i][j]) >= 0) return false;
		}
	}
	return true;
}

bool Matrix::operator!=(const Matrix &other) {
	if (row != other.row || col != other.col) return true;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			if (fabs(M[i][j] - other.M[i][j]) >= 0)
				return true;
	return false;
}

Matrix & Matrix::operator=(const Vector &vector){
	for (int i = 0; i < row; i++) delete[] M[i];
	delete[] M;
	row = 1;
	col = vector.dim;
	M = new double*[1];
	M[0] = new double[col];
	for (int i = 0; i < col; i++) M[0][i] = vector.vect[i];
	return *this;
}

Matrix & Matrix::operator=(const Matrix &other) {
	for (int i = 0; i < row; i++) delete[] M[i];
	delete[] M;
	row = other.row;
	col = other.col;
	M = new double*[row];
	for (int i = 0; i < row; i++) M[i] = new double[col];
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] = other.M[i][j];
	return *this;
}

Matrix Matrix::operator+(const Matrix &other) {
	if (row != other.row || col != other.col) throw 1;
	Matrix res(row, col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			res.M[i][j] = M[i][j] + other.M[i][j];
	return res;
}

Matrix Matrix::operator-(const Matrix &other) {
	if (row != other.row || col != other.col) throw 1;
	Matrix res(row, col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			res.M[i][j] = M[i][j] + other.M[i][j];
	return res;
}

Matrix Matrix::operator*(const Matrix &other) {
	if (row != other.col || col!= other.row) throw 1;
	Matrix res(row, other.col);
	for (int i = 0; i < res.row; i++)
		for (int j = 0; j < res.col; j++)
			for (int sum_index = 0; sum_index < col; sum_index++)
				res.M[i][j] += M[i][sum_index] * other.M[sum_index][j];
	return res;
}

Matrix Matrix::operator*(double koeff) {
	Matrix res(*this);
	for (int i = 0; i < res.row; i++){
		for (int j = 0; j < res.col; j++) {
			res.M[i][j] *= koeff;
		}
	}
	return res;
}

Matrix Matrix::operator/(double koeff) {
	Matrix res(*this);
	for (int i = 0; i < res.row; i++){
		for (int j = 0; j < res.col; j++){
			res.M[i][j] /= koeff;
		}
	}
	return res;
}

Matrix operator*(double koeff, const Matrix &right) {
	Matrix res(right);
	for (int i = 0; i < res.row; i++)
		for (int j = 0; j < res.col; j++)
			res.M[i][j] *= koeff;
	return res;
}

Vector Matrix::operator*(const Vector &right) {
	if (col != right.dimens()) throw 1;
	Vector res(row);
	for (int ivect = 0; ivect < row; ivect++)
		for (int isum = 0; isum < col; isum++)
			res.vect[ivect] += M[ivect][isum] * right.vect[isum];
	return res;
}

Matrix & Matrix::operator+=(const Matrix &other) {
	if (row != other.row || col != other.col) throw 1;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] += other.M[i][j];
	return *this;
}

Matrix & Matrix::operator-=(const Matrix &other) {
	if (row != other.row || col != other.col) throw 1;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] -= other.M[i][j];
	return *this;
}

Matrix & Matrix::operator*=(const Matrix &other) {
	*this = *this * other;
	return *this;
}

Matrix & Matrix::operator*=(double koeff) {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] *= koeff;
	return *this;
}

Matrix & Matrix::operator/=(double koeff) {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			M[i][j] /= koeff;
	return *this;
}

std::string Matrix::to_string() const {
	std::string res = row == 1 ? "->" : "<-";
	int to_ins_i_val;
	double to_ins_d_val;
	for (int i = 0; i < row; i++) {
		if (i != 0 && i != row - 1) res += "*";
		else if (i == row - 1) res += "#";
		for (int j = 0; j < col; j++) {
			to_ins_i_val = round(M[i][j]);
			to_ins_d_val = M[i][j];
			res += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + (j != col - 1 ? "\t\t" : "\t");
		}
		if (i != 0 && i != row - 1) res += "|#\n";
		else if (i == 0) res += "#|\n";
	}
	res += row == 1 ? "->\n\n" : "<-\n\n";
	return res;
}

Matrix Matrix::create_transposed() const {
	Matrix res(col, row);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			res.M[j][i] = M[i][j];
		}		
	}
	return res;
}

void Matrix::transpose() {
	double **trans_M = new double*[col];
	for (int i = 0; i < col; i++) trans_M[i] = new double[row];
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			trans_M[j][i] = M[i][j];
		}
	}	
	for (int i = 0; i < row; i++) delete[] M[i];
	delete[] M;
	M = trans_M;
}

Matrix::operator SqrMatrix() {
	if (row != col) throw 1;
	return SqrMatrix(M, row);
}

Matrix::operator Vector() { return Vector(*this); }

double Matrix::min_norm() const {
	if (cube_norm() <= octo_norm() && cube_norm() <= euclid_norm()) return cube_norm();
	else if (octo_norm() <= cube_norm() && octo_norm() <= euclid_norm()) return octo_norm();
	else return euclid_norm();
}

double Matrix::cube_norm() const {
	double norm = 0, hyp_max = 0;
	for (int j = 0; j < col; j++) norm += fabs(M[0][j]);
	for (int i = 1; i < row; i++) {
		hyp_max = 0;
		for (int j = 0; j < col; j++) hyp_max += fabs(M[i][j]);
		if (hyp_max > norm) norm = hyp_max;
	}
	return norm;
}

double Matrix::octo_norm() const {
	double norm = 0, hyp_max = 0;
	for (int i = 0; i < row; i++)
		norm += fabs(M[i][0]);
	for (int j = 1; j < col; j++) {
		hyp_max = 0;
		for (int i = 0; i < row; i++)
			hyp_max += fabs(M[i][j]);
		if (hyp_max > norm)
			norm = hyp_max;
	}
	return norm;
}

double Matrix::euclid_norm() const {
	double norm = 0;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			norm += M[i][j] * M[i][j];
	return sqrt(norm);
}

Matrix::~Matrix() {
	for (int i = 0; i < row; i++) delete[] M[i];
	delete[] M;
}


