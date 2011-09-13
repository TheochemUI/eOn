#include <stdio.h>
#include <stdarg.h>

static FILE *logfile=NULL;

void log(const char* format, ...) {
    va_list args;

    if (logfile == NULL) {
        logfile = fopen("client.log", "w");
        setvbuf(logfile, NULL, _IOLBF, 0);
    }

    va_start(args, format);
    vfprintf(logfile, format, args);
    vfprintf(stdout,  format, args);
    va_end(args);

    fprintf(stdout, "\n");
}
