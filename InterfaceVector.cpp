#include "InterfaceVector.h"

std::ostream & operator<<(std::ostream &os, const InterfaceVector &vect) {
	os << vect.to_string().c_str();
	return os;
}
