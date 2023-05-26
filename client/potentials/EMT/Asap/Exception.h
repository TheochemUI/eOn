#ifndef _EXCEPTION_H
#define _EXCEPTION_H

#include <sstream>
using std::stringstream;

// this Exception class can be used as follows:
//
// throw Exception("Error");
//
// throw Exception("Error ") << n + 7 << " while opening file: " << filename;

class Exception {
public:
  Exception(char const *m);
  Exception(const Exception &ex);
  ~Exception();
  template <class T> Exception &operator<<(const T &x) {
    *message << x;
    return *this;
  }
  const char *GetMessage() const;

private:
  stringstream *message;
};

#endif // _EXCEPTION_H
