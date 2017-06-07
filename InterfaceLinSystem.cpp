#include "InterfaceLinSystem.h"

std::ostream & operator<<(std::ostream &os, InterfaceLinSystem &sys) {
	return os << sys.to_string().c_str();
}
