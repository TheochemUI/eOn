#ifdef STACKTRACE

#ifdef __cplusplus
extern "C" {
#endif

typedef void abortfunc();
void setSignalHandlers(int node, abortfunc *ab);

#ifdef __cplusplus
}
#endif

#endif /* STACKTRACE */
