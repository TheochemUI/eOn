Replace per-site `scoped_interpreter` guards with lazy singleton `eonc::ensure_interpreter()` in `PyGuard.h`. Python interpreter is only started when a Python-based potential is actually used.
