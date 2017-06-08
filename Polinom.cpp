#include "Polinom.h"
#include <string>
#include <vector>

std::ostream & operator<<(std::ostream &os, Polinom &pol) {
	return os << pol.to_string();
}

Polinom::Polinom(int pol_pow) : Function(Vector(1, pol_pow + 1)) {
	this->pol_pow = pol_pow;
	pol_koef = Vector(1, pol_pow + 1);
}

Polinom::Polinom(Vector pol_koef) : Function(pol_koef){
	this->pol_koef = pol_koef;
	pol_pow = pol_koef.dimens() - 1;
}

Polinom::Polinom(const Polinom &other) : Function(other.pol_koef) {
	pol_koef = other.pol_koef;
	pol_pow = other.pol_pow;
}

Polinom::Polinom(double *pol_koef, int pol_pow) : Function(Vector(pol_koef, pol_pow + 1)){
	this->pol_pow = pol_pow;
	this->pol_koef = Vector(pol_koef, pol_pow + 1);
}

int Polinom::pow() const { return pol_pow; }

std::string Polinom::to_string() const {
	std::string res = "";
	int ival = round(pol_koef[0]);
	double dval = pol_koef[0];

	switch (pol_pow) {
	case 0:
		res = fabs(dval - ival) < EPS ? std::to_string(ival) : std::to_string(dval);
		break;
	case 1:
		res = (fabs(pol_koef[0] - 1) < EPS ? "" : fabs(pol_koef[0] + 1) < EPS ? "-" : (fabs(dval - ival) < EPS ? std::to_string(ival) : std::to_string(dval) + "*")) + "x";
		break;
	default:
		res = (fabs(pol_koef[0] - 1) < EPS ? "" : fabs(pol_koef[0] + 1) < EPS ? "-" : (fabs(dval - ival) < EPS ? std::to_string(ival) : std::to_string(dval) + "*")) + "x^" + std::to_string(pol_pow);
		break;
	}
	for (int i = 1; i <= pol_pow; i++) {
		ival = round(pol_koef[i]);
		dval = pol_koef[i];
		if (fabs(pol_koef[i]) < EPS) continue;
		else if (pol_koef[i] > EPS) res += "+" + (fabs(dval - ival) < EPS ? std::to_string(ival) : std::to_string(dval));
		else if (fabs(pol_koef[i] + 1) < EPS) {
			res += "-";
			if (pol_pow == i) res += "1";
		}
		else if (fabs(pol_koef[i] - 1) < EPS) {
			res += "+";
			if (pol_pow == i) res += "1";
		}
		else res += "-" + (fabs(dval - ival) < EPS ? std::to_string(-ival) : std::to_string(-dval));
		if (fabs(dval - ival) > EPS && i != pol_pow) res += "*";
		if (i != pol_pow) res += "x^" + std::to_string(pol_pow - i);
		else if (pol_pow - i == 1) res += "x";
	}
	return res;
}

double Polinom::operator()(double x) {
	double res = 0, tmp = 1;
	for (int i = pol_pow; i >= 0; i--) {
		res += pol_koef[i] * tmp;
		tmp *= x;
	}
	return res;
}

Polinom & Polinom::operator=(const Polinom &other) {
	pol_pow = other.pol_pow;
	pol_koef = other.pol_koef;
	return *this;
}

bool Polinom::operator==(const Polinom &other) { return pol_koef == other.pol_koef; }

bool Polinom::operator!=(const Polinom &other) { return pol_koef != other.pol_koef; }

Polinom Polinom::operator+(const Polinom &other) {
	int max_pow = pol_pow >= other.pol_pow ? pol_pow : other.pol_pow,
		this_index = pol_pow, other_index = other.pol_pow;
	Polinom result(max_pow);
	for (; other_index >= 0 && this_index >= 0; other_index--, this_index--)
		result.pol_koef[(this_index >= other_index ? this_index : other_index)] = pol_koef[this_index] + other.pol_koef[other_index];
	for (int i = (this_index >= 0 ? this_index : other_index); i >= 0; i--)
		result.pol_koef[i] = (this_index >= 0 ? pol_koef[i] : other.pol_koef[i]);
	return result;
}

Polinom Polinom::operator-(const Polinom &other) {
	int max_pow = pol_pow >= other.pol_pow ? pol_pow : other.pol_pow,
		this_index = pol_pow, other_index = other.pol_pow;
	Polinom res(max_pow);
	for (; other_index >= 0 && this_index >= 0; other_index--, this_index--)
		res.pol_koef[(this_index >= other_index ? this_index : other_index)] = pol_koef[this_index] - other.pol_koef[other_index];
	for (int i = (this_index >= 0 ? this_index : other_index); i >= 0; i--)
		res.pol_koef[i] = (this_index >= 0 ? pol_koef[i] : -other.pol_koef[i]);
	return res;
}

Polinom Polinom::operator*(const Polinom &other) {
	Polinom res(pol_pow + other.pol_pow);
	res.pol_koef = Vector((double)0, res.pol_pow + 1);
	for (int i = res.pol_pow; i >= 0; i--)
		for (int j = pol_pow; j >= 0; j--)
			for (int k = other.pol_pow; k >= 0; k--)
				if (j + k == i) res.pol_koef[i] += pol_koef[j] * other.pol_koef[k];
	return res;
}

Polinom Polinom::operator/(const Polinom &other) {
	if (other.pol_pow > pol_pow)
		return Polinom(Vector(1));
	std::vector<double> vect;
	Vector Vres(pol_pow - other.pol_pow + 1), Vtmp(pol_pow - other.pol_pow + 1);
	Polinom tmp(*this);
	for (int i = 0; i <= pol_pow - other.pol_pow; i++) {
		Vtmp = Vector(pol_pow - other.pol_pow + 1 - i);
		Vres[i] = tmp.pol_koef[0] / other.pol_koef[0];
		Vtmp[0] = Vres[i];
		tmp = tmp - Polinom(Vtmp) * other;
		vect.clear();
		for (int i = 0; i <= tmp.pol_pow; i++)
			vect.push_back(tmp.pol_koef[i]);
		vect.erase(vect.begin());
		tmp = Polinom(Vector(vect));
	}
	return Polinom(Vres);
}

Polinom Polinom::operator%(const Polinom &other) {
	if (other.pol_pow > pol_pow)
		return Polinom(Vector(1));
	std::vector<double> vect;
	Vector Vres(pol_pow - other.pol_pow + 1), Vtmp(pol_pow - other.pol_pow + 1);
	Polinom tmp(*this);
	for (int i = 0; i <= pol_pow - other.pol_pow; i++) {
		Vtmp = Vector(pol_pow - other.pol_pow + 1 - i);
		Vres[i] = tmp.pol_koef[0] / other.pol_koef[0];
		Vtmp[0] = Vres[i];
		tmp = tmp - Polinom(Vtmp) * other;
		if (tmp.pol_pow) {
			vect.clear();
			for (int j = 0; j <= tmp.pol_pow; j++)
				vect.push_back(tmp.pol_koef[j]);
			vect.erase(vect.begin());
			tmp = Polinom(Vector(vect));
		}
	}
	if (tmp.pol_pow) {
		vect.clear();
		for (int j = 0; j <= tmp.pol_pow; j++)
			vect.push_back(tmp.pol_koef[j]);
		while (!*vect.begin() && vect.begin() < vect.end() - 1)
			vect.erase(vect.begin());
		tmp = Polinom(Vector(vect));
	}
	return tmp;
}

Polinom Polinom::operator-() {
	Polinom res(*this);
	res.pol_koef *= -1;
	return res;
}

Polinom & Polinom::operator+=(const Polinom &other) {
	*this = *this + other;
	return *this;
}

Polinom & Polinom::operator-=(const Polinom &other) {
	*this = *this - other;
	return *this;
}

Polinom & Polinom::operator*=(const Polinom &other) {
	*this = *this * other;
	return *this;
}

Polinom & Polinom::operator/=(const Polinom &other) {
	*this = *this / other;
	return *this;
}

int Polinom::Sturm(double a, double b) { // a - левая граница b - правая граница
	int a_res = 0, b_res = 0;
	std::vector<Polinom> St_sys;
	St_sys.push_back(*this);
	St_sys.push_back(derivative(1));
	for (int i = 0; St_sys[i] % St_sys[i + 1] != Polinom(Vector(1)); i++)
		St_sys.push_back(-(St_sys[i] % St_sys[i + 1]));
	for (int i = 0; i < St_sys.size() - 1; i++)
		if (St_sys[i](b) * St_sys[i + 1](b) < 0) b_res++;
	for (int i = 0; i < St_sys.size() - 1; i++)
		if (St_sys[i](a) * St_sys[i + 1](a) < 0) a_res++;
	
	return a_res - b_res;
}

Polinom Polinom::derivative(int order) {
	if (order < 0) throw 1;
	else if (order > pol_pow) return Polinom(Vector(1));
	Polinom res(pol_pow - order);
	for (int i = 0; i <= res.pol_pow; i++) {
		res.pol_koef[i] = pol_koef[i];
		for (int j = i; j < order + i; j++)
			res.pol_koef[i] *= pol_pow - j;
	}
	return res;
}

double Polinom::deriv_n(double x, unsigned int order) { return derivative(order)(x); }

Polinom::~Polinom() {}

