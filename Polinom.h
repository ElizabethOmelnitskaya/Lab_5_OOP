#pragma once
#include "Function.h"
#include "Vector.h"

class Polinom : public Function {
private:
	int pol_pow;
	Vector pol_koef;
public:
	Polinom(int);
	Polinom(Vector);
	Polinom(const Polinom &);
	Polinom(double *, int);
	
	template<class list>
	Polinom(list);

	int pow() const;
	std::string to_string() const;

	double operator()(double);
	Polinom & operator=(const Polinom &);
	
	bool operator==(const Polinom &);
	bool operator!=(const Polinom &);
	Polinom operator+(const Polinom &);
	Polinom operator-(const Polinom &);
	Polinom operator*(const Polinom &);
	Polinom operator/(const Polinom &);
	Polinom operator%(const Polinom &);
	Polinom operator-();
	Polinom & operator+=(const Polinom &);
	Polinom & operator-=(const Polinom &);
	Polinom & operator*=(const Polinom &);
	Polinom & operator/=(const Polinom &);

	int Sturm(double, double);

	Polinom derivative(int); // производная n-го порядка  
	double deriv_n(double, unsigned int) override;
	
	~Polinom();
};

template<class init_list>
inline Polinom::Polinom(init_list) : Function(Vector(init_list)){
	pol_koef = Vector(init_list);
	pol_pow = pol_koef.dimens() - 1;
}

std::ostream& operator<<(std::ostream &, Polinom &);