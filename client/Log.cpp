#include <stdio.h>
#include <stdarg.h>
#include "Log.h"

static FILE *logfile=NULL;
static Parameters *params = NULL;

void log_init(Parameters *p, char *filename) {
    params = p;
    if (logfile == NULL) {
        if (params->checkpoint) {
            logfile = fopen(filename, "a");
        }else{
            logfile = fopen(filename, "w");
        }
        //setvbuf(logfile, NULL, _IOLBF, 0);
    }
}

void log_close() {
    if (logfile != NULL) {
        fclose(logfile);
        logfile = NULL;
    }
}

void log(const char* format, ...) {
    if (logfile == NULL) {
        fprintf(stderr, "error: log() called before log_init\n");
        return;
    }
    va_list args;
    if (params->quiet == false) {
        va_start(args, format);
        vfprintf(stdout, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        //fflush(logfile);
    }

    if (params->writeLog == true) {
        va_start(args, format);
        vfprintf(logfile, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        //fflush(logfile);
    }
}

void log_file(const char* format, ...) {
    if (logfile == NULL) {
        fprintf(stderr, "error: log() called before log_init\n");
        return;
    }
    if (params->writeLog == true) {
        va_list args;
        va_start(args, format);
        vfprintf(logfile, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        //fflush(logfile);
    }
}
