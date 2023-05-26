#include "Exception.h"

Exception::Exception(char const *m) : message(new stringstream()) {
  *message << m;
}

Exception::Exception(const Exception &ex) : message(new stringstream()) {
  *message << ex.GetMessage();
}

Exception::~Exception() { delete message; }
const char *Exception::GetMessage() const { return message->str().c_str(); }
