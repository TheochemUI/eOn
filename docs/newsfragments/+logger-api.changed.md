Migrate to modern C++20 logging API (EonLogger.h).

Replaced verbose quill logger initialization throughout codebase with new `eonc::log::Scoped` RAII helper and `eonc::log::get_file()` utility. Eliminates 14-line boilerplate for file loggers and manual initialization for default loggers. Net reduction of 74 lines while improving code clarity and safety.
