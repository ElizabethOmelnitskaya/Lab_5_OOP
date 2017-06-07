#pragma once
#include "InterfaceFanction.h"
#include "Vector.h"

class Function : public InterfaceFunction {
public:
	typedef double(*func)(double);			//������ �������������� ������� ����� ����������
	typedef double(*n_func)(double, unsigned int);	//������ ������������������ �������������� ������� ����� ����������
	const double pi = 3.14, e = 2.7, CONST_EPS = 0.001, good = (1 + sqrt(5)) / 2, H_INTEGR_CONST = 0.01;
private:
	Vector Vpol;
	func suitable; // ���������� ������� 
protected:
	Function(Vector);
public:
	Function(func);
	Function(const Function &);
	double operator()(double);
	Function & operator=(const Function &);
	
	virtual double deriv_n(double, unsigned int); // ����������� n-�� ������� 
	double min(double, double);
	double max(double, double);
	
	double abs_min(double, double);
	double abs_max(double, double);
	
	double n_deriv_min(unsigned int, double, double);
	double n_deriv_max(unsigned int, double, double);

	double abs_n_deriv_min(unsigned int, double, double);
	double abs_n_deriv_max(unsigned int, double, double);

	short int signum(double, double, double);
	short int n_deriv_signum(unsigned int, double, double, double);

	void half_devision(double, double, double, double &, int &, double &, int = 100, double = 0.001); //����� ����������� �������
	void iteration(double, double, double, double &, int &, double &, int = 100, double = 0.001); // ����� ������� ��������
	void chords(double, double, double, double &, int &, double &, int = 100, double = 0.001); // ����� ����
	void tangents(double, double, double, double &, int &, double &, int = 100, double = 0.001); // ����� �����������

	~Function();
};