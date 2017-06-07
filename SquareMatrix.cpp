#include "SquareMatrix.h"

SqrMatrix::SqrMatrix() : Matrix() {}
SqrMatrix::SqrMatrix(int matrix_dim) : Matrix(matrix_dim, matrix_dim) { this->m_dim = matrix_dim; }
SqrMatrix::SqrMatrix(double **M, int m_dim) : Matrix(M, m_dim, m_dim) { this->m_dim = m_dim; }
SqrMatrix::SqrMatrix(double def_val, int m_dim) : Matrix(def_val, m_dim, m_dim) { this->m_dim = m_dim; }

SqrMatrix::SqrMatrix(const Matrix &other) : Matrix(other) {
	if (other.cols() != other.rows()) throw 1;
	m_dim = other.rows();
}
SqrMatrix::SqrMatrix(const SqrMatrix &other) : Matrix(other) { m_dim = other.m_dim; }

bool SqrMatrix::symmetry() const {
	for (int i = 0; i < m_dim; i++){
		for (int j = i + 1; j < m_dim; j++){
			if (fabs(M[i][j] - M[j][i]) >= EPS) return false;
		}
	}
	return true;
}

SqrMatrix & SqrMatrix::operator=(const Matrix &other) {
	if (other.cols() != other.rows()) throw 1;
	Matrix tmp(other);
	for (int i = 0; i < m_dim; i++) delete[] M[i];
	delete[] M;
	m_dim = tmp.rows();
	M = new double*[m_dim];
	for (int i = 0; i < m_dim; i++)
		M[i] = new double[m_dim];
	for (int i = 0; i < m_dim; i++)
		for (int j = 0; j < m_dim; j++)
			M[i][j] = tmp[i][j];
	return *this;
}

int SqrMatrix::dimens() const { return m_dim; }

bool SqrMatrix::minors() const {
	SqrMatrix M_tmp(M, 1);
	for (int i = 1; i < m_dim; i++, M_tmp = SqrMatrix(M, i)){
		if (fabs(M_tmp.determ()) < EPS) return false;
	}
	return true;
}

bool SqrMatrix::diagonal_dominating() const {
	double max_row = 0; // max in row abs 
	for (int i = 0; i < m_dim; i++) {
		max_row = fabs(M[i][0]);
		for (int j = 1; j < m_dim; j++)
			if (max_row < fabs(M[i][j])) max_row = fabs(M[i][j]);
		if (fabs(max_row - M[i][i]) > EPS) return false;
	}
	return true;
}

double SqrMatrix::determ() const {
	double det = 1, tmp = 0;
	int trans = 0;
	SqrMatrix copy(*this);
	for (int idiag = 0, max_irow; idiag < m_dim; idiag++) {
		if (fabs(copy.M[idiag][idiag]) < EPS) {
			max_irow = idiag;
			for (int j = idiag + 1; j < m_dim; j++) {
				if (fabs(copy.M[idiag][j]) > fabs(copy.M[idiag][idiag])) max_irow = j;
			}
			if (fabs(copy.M[idiag][max_irow]) < EPS) return 0;
			trans++;
			for (int i = 0; i < m_dim; i++) {
				tmp = copy.M[i][max_irow];
				copy.M[i][max_irow] = copy.M[i][idiag];
				copy.M[i][idiag] = tmp;
			}
		}
		for (int i = idiag + 1; i < m_dim; i++) {
			tmp = copy.M[i][idiag];
			for (int j = idiag; j < m_dim; j++) copy.M[i][j] -= copy.M[idiag][j] / copy.M[idiag][idiag] * tmp;
		}
	}
	for (int idiag = 0; idiag < m_dim; idiag++) det *= copy.M[idiag][idiag];
	det *= pow(-1, trans);
	return det;
}

SqrMatrix SqrMatrix::cr_revers() const {
	double tmp = 0;
	SqrMatrix res(m_dim), copy(*this);

	for (int i = 0; i < m_dim; i++) res.M[i][i] = 1;

	for (int idiag = 0, max_icol; idiag < m_dim; idiag++) {
		if (fabs(copy.M[idiag][idiag]) < EPS) {
			max_icol = idiag;
			for (int i = idiag + 1; i < m_dim; i++) {
				if (fabs(copy.M[i][idiag]) > fabs(copy.M[max_icol][idiag])) max_icol = i;
			}
			if (fabs(copy.M[max_icol][idiag]) < EPS) throw 1;
			for (int j = 0; j < m_dim; j++) {
				tmp = copy.M[max_icol][j];
				copy.M[max_icol][j] = copy.M[idiag][j];
				copy.M[idiag][j] = tmp;
				tmp = res.M[max_icol][j];
				res.M[max_icol][j] = res.M[idiag][j];
				res.M[idiag][j] = tmp;
			}
		}
		for (int i = idiag + 1; i < m_dim; i++) {
			tmp = copy.M[i][idiag];
			for (int j = idiag; j < m_dim; j++)
				copy.M[i][j] -= copy.M[idiag][j] / copy.M[idiag][idiag] * tmp;
			for (int j = 0; j < m_dim; j++)
				res.M[i][j] -= res.M[idiag][j] / copy.M[idiag][idiag] * tmp;
		}
	}

	for (int idiag = m_dim - 1; idiag >= 0; idiag--) {
		for (int i = idiag - 1; i >= 0; i--) {
			tmp = copy.M[i][idiag];
			for (int j = 0; j < m_dim; j++) {
				copy.M[i][j] -= copy.M[idiag][j] / copy.M[idiag][idiag] * tmp;
				res.M[i][j] -= res.M[idiag][j] / copy.M[idiag][idiag] * tmp;
			}
		}
	}

	for (int i = 0; i < m_dim; i++) {
		tmp = copy.M[i][i];
		for (int j = 0; j < m_dim; j++) {
			copy.M[i][j] /= tmp;
			res.M[i][j] /= tmp;
		}
	}
	return res;
}

SqrMatrix::~SqrMatrix() {}
