#include <analysis/modules/DecayFinder/Descriptor.h>
using namespace Belle2;

#include <string>

// A || B

bool OrDescriptor::match(int index) const {
  return left->match(index) || right->match(index);
}

std::string OrDescriptor::repr() const {
  return "(" + left->repr() + " || " + right->repr() + ")";
}

// A && B

bool AndDescriptor::match(int index) const {
  return left->match(index) && right->match(index);
}

std::string AndDescriptor::repr() const {
  return "(" + left->repr() + " && " + right->repr() + ")";
}
