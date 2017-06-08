#define _CRT_SECURE_NO_WARNINGS
#include "Polinom.h"
#include "Function.h"
#include <cmath>
#include <ctime>
#include "SquareMatrix.h"
#include "LinSystem.h"
#include <string>

double func(double x) { return (x - 2) * (x - 1); }
using namespace std;

void main() {
	srand(time(0));
	setlocale(LC_ALL, "ru");

	int iters = 0, pow;
	double sol = 0, delt = 0;
	cout << "Input power: ";
	cin >> pow;
	cout << endl;
	Polinom pol(pow);
	Vector Vtmp(pow + 1);
	for (int i = 0; i <= pow; i++) {
		cout << "Input koef for x " << i << ": ";
		cin >> Vtmp[i];
	}
	cout << endl;
	pol = Polinom(Vtmp);
	cout << "The introduced polynomial: " << pol << endl;
	cout << endl;
	cout << endl;
	for (int i = 0; i <= pow; i++)
		cout << "The derivative is of the order of: " << i << endl << pol.derivative(i) << "\n";
	cout << endl;
	cout << "The number of decisions in a given interval : " << pol.Sturm(-10, 200) << "\n" << endl;
	
	system("pause");
}