#pragma once
#define EPS 0.0001
#include <iostream>

class InterfaceVector {
public:
	virtual double cube_norm() const = 0;
	virtual double octo_norm() const = 0;
	virtual double euclid_norm() const = 0;
	virtual std::string to_string() const = 0;
};

std::ostream& operator<<(std::ostream&, const InterfaceVector&);
