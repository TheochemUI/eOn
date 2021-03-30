#ifndef EON_LOG_H
#define EON_LOG_H

#include "Parameters.h"

void log_init(Parameters *p, char *filename);

void log_close();

void log(const char* format, ...);

void log_file(const char* format, ...);

#endif
