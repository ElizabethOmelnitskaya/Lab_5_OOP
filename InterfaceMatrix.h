#pragma once
#include "InterfaceVector.h"

class InterfaceMatrix : public InterfaceVector {
public:
	virtual int rows() const = 0;
	virtual int cols() const = 0;
	virtual void transpose() = 0;
};
