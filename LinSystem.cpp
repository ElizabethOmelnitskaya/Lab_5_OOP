#include "LinSystem.h"
#include <string>

LinSystem::LinSystem(SqrMatrix sys_matrix, Vector free) {
	if (sys_matrix.dimens() != free.dimens()) throw 1;
	this->sys_matrix = sys_matrix;
	this->fill_koef = free;
	sysDim = sys_matrix.dimens();
	is_certain = fabs(sys_matrix.determ()) >= EPS;
	if (is_certain) optimal();
}

LinSystem::LinSystem(double **sys_extended, int sys_dim) { //extended_system_matrix
	this->sysDim = sys_dim;
	sys_matrix = SqrMatrix(sys_extended, sys_dim);
	fill_koef = Vector(sys_dim);
	for (int i = 0; i < sys_dim; i++) fill_koef[i] = sys_extended[i][sys_dim];
	is_certain = fabs(sys_matrix.determ()) >= EPS;
	if (is_certain) optimal();
}

LinSystem::LinSystem(double **system_matrix, double *free_koeffs_vector, int sys_dim) {
	this->sysDim = sys_dim;
	this->sys_matrix = SqrMatrix(system_matrix, sys_dim);
	this->fill_koef = Vector(free_koeffs_vector, sys_dim);
	is_certain = fabs(this->sys_matrix.determ()) >= EPS;
	if (is_certain) optimal();
}

SqrMatrix LinSystem::matrix() { return sys_matrix; }

void LinSystem::matrix_method() {
	if (is_certain) vdecision = sys_matrix.cr_revers() * fill_koef;
	else throw 1;
}

Vector LinSystem::decision() const {
	if (!is_certain) throw 1;
	return vdecision;
}

double LinSystem::error() const { return err; }

bool LinSystem::certain() const { return is_certain; }

int LinSystem::dimens() const { return sysDim; }

Vector LinSystem::fill_vector() { return fill_koef; }

std::string LinSystem::to_string(){
	std::string res = "SLAU:\n";
	int i_val;
	double d_val;

	for (int i = 0, first_i; i < sysDim; i++) {
		first_i = -1;
		for (int j = 0; j < sysDim; j++)
			if (fabs(sys_matrix[i][j]) > EPS) {
				first_i = j;
				break;
			}
		for (int j = 0; j < sysDim && first_i >= 0; j++) {
			i_val = round(sys_matrix[i][j]);
			d_val = sys_matrix[i][j];
			if (j == first_i) res += (fabs(sys_matrix[i][j] - 1) < EPS ? "" : fabs(sys_matrix[i][j] + 1) < EPS ? "-" : (fabs(d_val - i_val) < EPS ? std::to_string(i_val) : std::to_string(d_val) + " * ")) + "x" + std::to_string(j + 1);
			else if (sys_matrix[i][j] > EPS) res += " + " + (fabs(d_val - i_val) < EPS ? std::to_string(i_val) : std::to_string(d_val) + " & ") + "x" + std::to_string(j + 1);
			else if (sys_matrix[i][j] < -EPS) res += " - " + (fabs(d_val - i_val) < EPS ? std::to_string(-i_val) : std::to_string(-d_val) + " & ") + "x" + std::to_string(j + 1);
			else if (fabs(sys_matrix[i][j] - 1) < EPS) res += " + x" + std::to_string(j + 1);
			else if (fabs(sys_matrix[i][j] + 1) < EPS) res += " - x" + std::to_string(j + 1);
		}
		d_val = fill_koef[i];
		i_val = round(fill_koef[i]);
		res += " = " + (fabs(d_val - i_val) < EPS ? std::to_string(i_val) : std::to_string(d_val)) + "\n\n";
	}
	res += "Result:";
	if (!is_certain) res += " Once res or NOT res!";
	else
		for (int i = 0; i < sysDim; i++) {
			d_val = vdecision[i];
			i_val = round(vdecision[i]);
			res += "\n\nx" + std::to_string(i + 1) + " = " + (fabs(d_val - i_val) < EPS ? std::to_string(i_val) : std::to_string(d_val)) + ";";
		}
	return res += "\n\n EPS < " + std::to_string(err) + "\n\n"; //погрешность не превышает
}

Matrix LinSystem::extended_matrix() { 
	Matrix res(sysDim, sysDim + 1);
	for (int i = 0; i < sysDim; i++)
		for (int j = 0; j < sysDim; j++)
			res[i][j] = sys_matrix[i][j];
	for (int i = 0; i < sysDim; i++)
		res[i][sysDim] = fill_koef[i];
	return res;
}

Matrix LinSystem::sysswap() {
	double tmp = 0;
	Matrix result(extended_matrix());
	for (int main_diag_i = 0, abs_lower_max_i; main_diag_i < result.rows(); main_diag_i++) {
		abs_lower_max_i = main_diag_i;
		for (int i = main_diag_i + 1; i < result.rows(); i++)
			if (fabs(result[i][main_diag_i]) > fabs(result[abs_lower_max_i][main_diag_i]))
				abs_lower_max_i = i;
		if (main_diag_i != abs_lower_max_i)
			for (int j = 0; j < result.cols(); j++) {
				tmp = result[main_diag_i][j];
				result[main_diag_i][j] = result[abs_lower_max_i][j];
				result[abs_lower_max_i][j] = tmp;
			}

	}
	return result;
}

void LinSystem::optimal() {
	if (!Zeidel()) {
		Gauss();
		err = 0;
	}
	if (!is_certain) throw 1;
}

/*void LinearSquareSystem::Cramer_method() {
	if (!is_determined)
		throw 1;
	solution_vector = MathVector(sys_dim);
	double tmp_det = system_matrix.determinant();
	SquareMatrix tmp_system_matrix(system_matrix);
	for (int j = 0; j < sys_dim; j++) {
		tmp_system_matrix = system_matrix;
		for (int i = 0; i < sys_dim; i++)
			tmp_system_matrix[i][j] = free_koeffs_vector[i];
		solution_vector[j] = tmp_system_matrix.determinant() / tmp_det;
	}
}*/

void LinSystem::Gauss() {
	if (!is_certain)
		throw 1;
	double tmp = 0;
	Matrix tmp_extend(extended_matrix()); // tmp_extend - текущая расширенная матрица
	for (int idiag = 0, low_imax; idiag < sysDim; idiag++, low_imax = idiag) { //low_imax - нижний по модулю макс
		if (fabs(tmp_extend[idiag][idiag]) < EPS) {
			low_imax = idiag;
			for (int i = idiag + 1; i < sysDim; i++){
				if (fabs(tmp_extend[i][idiag]) > fabs(tmp_extend[low_imax][idiag])) low_imax = i;
			}
			if (idiag != low_imax)
				for (int j = 0; j < sysDim + 1; j++) {
					tmp = tmp_extend[idiag][j];
					tmp_extend[idiag][j] = tmp_extend[low_imax][j];
					tmp_extend[low_imax][j] = tmp;
				}
		}
		if (fabs(tmp_extend[idiag][idiag]) < EPS)
			throw 1;
		for (int i = idiag + 1; i < sysDim; i++) {
			tmp = tmp_extend[i][idiag];
			for (int j = idiag; j < sysDim + 1; j++)
				tmp_extend[i][j] -= tmp_extend[idiag][j] / tmp_extend[idiag][idiag] * tmp;
		}
	}

	for (int sol_i = sysDim - 1; sol_i >= 0; sol_i--) {
		vdecision[sol_i] = tmp_extend[sol_i][sysDim];
		for (int i = sysDim - 1; i > sol_i; i--)
			vdecision[sol_i] -= vdecision[i] * tmp_extend[sol_i][i];
		vdecision[sol_i] /= tmp_extend[sol_i][sol_i];
	}
}

void LinSystem::LU() {
	if (!sys_matrix.minors() || !is_certain) throw 1;
	Vector tmp_vector(sysDim);
	SqrMatrix upp(sysDim), low(sysDim);
	
	for (int i = 0; i < sysDim; i++) {
		tmp_vector[i] = 0;
		for (int j = 0; j < sysDim; j++) {
			low[i][j] = 0;
			if (i != j) upp[i][j] = 0;
			else upp[i][j] = 1;
		}
	}
	for (int i = 0; i < sysDim; i++) {
		for (int j = i + 1; j < sysDim; j++) {
			upp[i][j] = sys_matrix[i][j];
			for (int isum = 0; isum < i; isum++)
				upp[i][j] -= low[i][isum] * upp[isum][j];
			upp[i][j] /= low[i][i];
		}
		for (int j = i; j < sysDim; j++) {
			low[j][i] = sys_matrix[j][i];
			for (int isum = 0; isum < i; isum++)
				low[j][i] -= low[j][isum] * upp[isum][i];
		}
	}
	for (int i = 0; i < sysDim; i++) {
		tmp_vector[i] = fill_koef[i];
		for (int j = 0; j < i; j++) tmp_vector[i] -= tmp_vector[j] * low[i][j];
		tmp_vector[i] /= low[i][i];
	}
	for (int i = sysDim - 1; i >= 0; i--) {
		vdecision[i] = tmp_vector[i];
		for (int j = sysDim - 1; j > i; j--)
			vdecision[i] -= vdecision[j] * upp[i][j];
	}
}

void LinSystem::square_root() {
	if (!sys_matrix.minors() || !is_certain || !sys_matrix.symmetry()) throw 1;
	SqrMatrix sqrM(sysDim);
	Vector Vtmp(sysDim);
	int *sgn = new int[sysDim]; // вектор
	for (int i = 0; i < sysDim; i++){
		for (int j = 0; j < sysDim; j++){
			sqrM[i][j] = 0;
		}
	}
	for (int j = 0; j < sysDim; j++) {
		for (int i = 0; i < j; i++) {
			sqrM[i][j] = sys_matrix[i][j];
			for (int k = 0; k < i; k++) sqrM[i][j] -= sqrM[k][i] * sqrM[k][j] * sgn[k];
			sqrM[i][j] /= sqrM[i][i] * sgn[i];
		}
		sqrM[j][j] = sys_matrix[j][j];
		for (int k = 0; k < j; k++) sqrM[j][j] -= sgn[k] * sqrM[k][j] * sqrM[k][j];
		sgn[j] = sqrM[j][j];
		sqrM[j][j] = sqrt(fabs(sqrM[j][j]));
		sgn[j] = sgn[j] > 0 ? 1 : sgn[j] < 0 ? -1 : 0;
	}
	for (int i = 0; i < sysDim; i++) {
		Vtmp[i] = fill_koef[i];
		for (int j = 0; j < i; j++) Vtmp[i] -= Vtmp[j] * sqrM[j][i];
		Vtmp[i] /= sqrM[i][i];
	}
	for (int i = sysDim - 1; i >= 0; i--) {
		vdecision[i] = Vtmp[i];
		for (int j = sysDim - 1; j > i; j--) vdecision[i] -= vdecision[j] * sqrM[i][j] * sgn[i];
		vdecision[i] /= sqrM[i][i] * sgn[i];
	}

	delete[] sgn;
}

bool LinSystem::simple_iteration() {
	if (!is_certain) throw 1;
	bool stop = false;
	unsigned short int norm_num = 0;
	double tmp = 0;
	Vector Vbetta(sysDim), tmp_vector(sysDim);
	SqrMatrix Malpha(sysDim);
	Matrix c_matrix = sys_matrix.diagonal_dominating() ? extended_matrix() : sysswap();
	for (int i = 0; i < sysDim; i++) {
		for (int j = 0; j < sysDim; j++) {
			if (i == j) Malpha[i][j] = 0;
			else Malpha[i][j] = -c_matrix[i][j] / c_matrix[i][i];
		}
		Vbetta[i] = fill_koef[i] / c_matrix[i][i];
	}
	if (Malpha.min_norm() >= 1) return false;
	for (int i = 0; i < sysDim; i++) vdecision[i] = Vbetta[i];
	if (Malpha.cube_norm() <= Malpha.octo_norm() && Malpha.cube_norm() <= Malpha.euclid_norm()) norm_num = 0;
	else if (Malpha.octo_norm() <= Malpha.cube_norm() && Malpha.octo_norm() <= Malpha.euclid_norm()) norm_num = 1;
	else norm_num = 2;
	do {
		for (int i = 0; i < sysDim; i++) {
			tmp_vector[i] = Vbetta[i];
			for (int j = 0; j < sysDim; j++) tmp_vector[i] += Malpha[i][j] * vdecision[j];
		}
		for (int i = 0; i < sysDim; i++) vdecision[i] -= tmp_vector[i];
		err = Malpha.min_norm() / (1 - Malpha.min_norm());
		err *= !norm_num ? vdecision.cube_norm() : norm_num == 1 ? vdecision.octo_norm() : vdecision.euclid_norm();
		if (err < EPS) stop = true;
		for (int i = 0; i < sysDim; i++) vdecision[i] = tmp_vector[i];
	} while (!stop);
	return true;
}

bool LinSystem::Zeidel() {
	bool stop = false;
	unsigned short int norm_num = 0;
	double tmp = 0;
	Vector Vbetta(sysDim), Vtmp(sysDim);
	SqrMatrix Malpha(sysDim);

	if (!is_certain) throw 1;
	
	Matrix c_matrix = sys_matrix.diagonal_dominating() ? extended_matrix() : sysswap();
	for (int i = 0; i < sysDim; i++) {
		for (int j = 0; j < sysDim; j++) {
			if (i == j) Malpha[i][j] = 0;
			else Malpha[i][j] = -c_matrix[i][j] / c_matrix[i][i];
		}
		Vbetta[i] = fill_koef[i] / c_matrix[i][i];
	}
	if (Malpha.min_norm() >= 1) return false;
	for (int i = 0; i < sysDim; i++) vdecision[i] = Vbetta[i];
	if (Malpha.cube_norm() <= Malpha.octo_norm() && Malpha.cube_norm() <= Malpha.euclid_norm()) norm_num = 0;
	else if (Malpha.octo_norm() <= Malpha.cube_norm() && Malpha.octo_norm() <= Malpha.euclid_norm()) norm_num = 1;
	else norm_num = 2;
	do {
		for (int i = 0; i < sysDim; i++) {
			Vtmp[i] = Vbetta[i];
			for (int j = i; j < sysDim; j++)
				Vtmp[i] += Malpha[i][j] * vdecision[j];
			for (int j = 0; j < i; j++)
				Vtmp[i] += Malpha[i][j] * Vtmp[j];
		}
		for (int i = 0; i < sysDim; i++) vdecision[i] -= Vtmp[i];
		err = Malpha.min_norm() / (1 - Malpha.min_norm());
		err *= !norm_num ? vdecision.cube_norm() : norm_num == 1 ? vdecision.octo_norm() : vdecision.euclid_norm();
		if (err < EPS) stop = true;
		for (int i = 0; i < sysDim; i++)
			vdecision[i] = Vtmp[i];

	} while (!stop);
	return true;
}

LinSystem::~LinSystem() {}
