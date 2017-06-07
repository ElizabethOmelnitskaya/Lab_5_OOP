#pragma once
#include "Vector.h"

class InterfaceLinSystem {
public:
	virtual std::string to_string() = 0;
	virtual int dimens() const = 0;
	virtual double error() const = 0;
	virtual Vector decision() const = 0; // решение(solution)
	
};

std::ostream& operator<<(std::ostream &, InterfaceLinSystem &);