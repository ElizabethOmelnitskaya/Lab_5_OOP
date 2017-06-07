#include "Function.h"
#include <vector>

Function::Function(func suitable) { this->suitable = suitable; }

Function::Function(Vector Vpol) {
	this->Vpol = Vpol;
	suitable = NULL;
}

Function::Function(const Function& other) { suitable = other.suitable; }

double Function::operator()(double x) {
	if (suitable) return suitable(x);
	double res = 0, tmp_comp = 1;
	for (int i = Vpol.dimens() - 1; i >= 0; i--) {
		res += Vpol[i] * tmp_comp;
		tmp_comp *= x;
	}
	return res;
}

Function & Function::operator=(const Function& other) {
	suitable = other.suitable;
	return *this;
}

double Function::deriv_n(double x, unsigned int order) {
	double eps = order <= 4 ? eps * pow(10, order - 1) : 0.1;
	std::vector<double> func_derivs;
	if (order == 0) return (*this)(x);
	
	for (int i = 0; i <= order; i++) { func_derivs.push_back((*this)(x + ((int)order - 2 * i) * eps)); }

	for (int order_index = 0; order_index < order; order_index++) {
		for (int i = 0; i < order - order_index; i++){
			func_derivs[i] = (func_derivs[i] - func_derivs[i + 1]) / (2 * eps);
		}
		func_derivs.pop_back();
	}
	return func_derivs[0];
}

double Function::min(double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if ((*this)(x1) < (*this)(x2)) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else  a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
} 

double Function::max(double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if ((*this)(x1) > (*this)(x2)) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good; 
	}
	return (a + b) / 2;
}

double Function::n_deriv_min(unsigned int order, double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (deriv_n(x1, order) < deriv_n(x2, order)) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

double Function::n_deriv_max(unsigned int order, double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (deriv_n(x1, order) > deriv_n(x2, order)) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

double Function::abs_min(double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (fabs((*this)(x1)) < fabs((*this)(x2))) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

double Function::abs_max(double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (fabs((*this)(x1)) > fabs((*this)(x2))) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

double Function::abs_n_deriv_min(unsigned int order, double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (fabs(deriv_n(x1, order)) < fabs(deriv_n(x2, order))) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

double Function::abs_n_deriv_max(unsigned int order, double a, double b) {
	double x1 = 0, x2 = 0;
	if (a > b) {
		double tmp = a;
		a = b;
		b = tmp;
	}
	while (fabs(b - a) > CONST_EPS) {
		x1 = b - (b - a) / good;
		x2 = a + (b - a) / good;
		if (fabs(deriv_n(x1, order)) > fabs(deriv_n(x2, order))) b = x2, x2 = x1, x1 = b - (b - a) / good;
		else  a = x1, x1 = x2, x2 = a + (b - a) / good;
	}
	return (a + b) / 2;
}

short int Function::signum(double a, double b, double step = 0.01) {
	bool pos = false, neg = false;
	for (double x = a; x <= b; x += step) {
		if ((*this)(x) < -CONST_EPS) neg = true;
		else if ((*this)(x) > CONST_EPS) pos = true;
	}
	return (pos && neg) ? 2 : (pos ? 1 : (neg ? -1 : 0));
}

short int Function::n_deriv_signum(unsigned int order, double a, double b, double step = 0.01) {
	bool pos = false, neg = false;
	for (double x = a; x <= b; x += step) {
		if (deriv_n(x, order) < -CONST_EPS) neg = true;
		else if (deriv_n(x, order) > CONST_EPS) pos = true;
	}
	return (pos && neg) ? 2 : (pos ? 1 : (neg ? -1 : 0));
}

void Function::half_devision(double func_val, double a, double b, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (((*this)(a)-func_val) * ((*this)(b)-func_val) >= 0) throw 1;
	solve = (a + b) / 2; iters = 0;
	while (iters < iters_lim_num && (b - a) / 2 > eps) {
		if (((*this)(solve)-func_val) * ((*this)(a)-func_val) < 0) b = (a + b) / 2;
		else a = (a + b) / 2;

		solve = (a + b) / 2;
		iters++;
	}
	abs_err = (a - b) / 2;
}

void Function::iteration(double func_val, double a, double b, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (((*this)(a)-func_val) * ((*this)(b)-func_val) > 0) throw 1;
	else if (n_deriv_signum(1, a, b) == 2) throw 2;
	double iter_koeff = 2 * n_deriv_signum(1, a, b) / (fabs(deriv_n(abs_n_deriv_max(1, a, b), 1)) + fabs(deriv_n(abs_n_deriv_min(1, a, b), 1)));
	double abs_err_koeff = fabs(deriv_n(abs_n_deriv_max(1, a, b), 1)) * fabs(iter_koeff) - 1;
	solve = (a + b) / 2;
	do {
		abs_err = fabs(((*this)(solve)-func_val) * iter_koeff * abs_err_koeff / (1 - abs_err_koeff));
		solve = solve - ((*this)(solve)-func_val) * iter_koeff;
		iters++;
	} while (fabs((*this)(solve)-func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

void Function::chords(double func_val, double a, double b, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (n_deriv_signum(1, a, b) == 2 || n_deriv_signum(2, a, b) == 2) throw 1;
	iters = 0;
	solve = ((*this)(b)-func_val) * n_deriv_signum(2, a, b) > 0 ? a : b;
	if (fabs(solve - a) < eps)
		do
		{
			solve = b - ((*this)(b)-func_val) * (b - solve) / ((*this)(b)-(*this)(solve));
			abs_err = fabs(b - ((*this)(b)-func_val) * (b - solve) / ((*this)(b)-(*this)(solve)) - solve);
			iters++;
		} while (fabs((*this)(solve)-func_val) > eps && iters < iters_lim_num && abs_err > eps);
	else
		do
		{
			solve = a - ((*this)(a)-func_val) * (solve - a) / ((*this)(solve)-(*this)(a));
			abs_err = fabs(a - ((*this)(a)-func_val) * (solve - a) / ((*this)(solve)-(*this)(a)) - solve);
			iters++;
		} while (fabs((*this)(solve)-func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

void Function:: tangents(double func_val, double a, double b, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (n_deriv_signum(1, a, b) == 2 || n_deriv_signum(2, a, b) == 2) throw 1;
	solve = ((*this)(b)-func_val) * deriv_n(b, 2) > 0 ? b : a;
	iters = 0;
	do {
		solve = solve - ((*this)(solve)-func_val) / deriv_n(solve, 1);
		abs_err = fabs(((*this)(solve)-func_val) / deriv_n(solve, 1));
		iters++;
	} while (fabs((*this)(solve)-func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

Function::~Function() {}
