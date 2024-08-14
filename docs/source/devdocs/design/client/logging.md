# Logging

```{versionchanged} 3.x
```

`spdlog` is used for managing the logs and for writing some of the output files.
This is handled through an RAII compatible interface, since the logger is
effectively a singleton.


## Singleton `LogManager` Class

The `LogManager` class is designed to be a singleton, meaning that only one
instance of it exists during the application's runtime. This instance is
responsible for setting up the loggers, managing log sinks, and ensuring that
all log messages are appropriately handled.

```{code-block} cpp
class LogManager {
public:
    static LogManager *getInstance();

private:
    void setup_logger();
    void cleanup_logger();
    LogManager() { setup_logger(); }
    ~LogManager() { cleanup_logger(); }
    LogManager(const LogManager &) = delete;
    LogManager &operator=(const LogManager &) = delete;

    static std::atomic<LogManager *> instance;
    static std::mutex lmMutex;
};
```

### Implementation Details

- **Thread-Safe Singleton Initialization:** The LogManager class uses a
  combination of `std::atomic` and `std::mutex` to ensure that the singleton
  instance is created in a thread-safe manner. This is crucial in multi-threaded
  applications where multiple threads might try to initialize the logger
  simultaneously.

- **Setup and Cleanup:** The `setup_logger` method configures the loggers,
  setting up various sinks such as console and file sinks, and registers these
  loggers with spdlog. The `cleanup_logger` method ensures that all loggers are
  properly flushed and shutdown, releasing any resources they might be holding.

- **RAII Compliance:** By constructing the logger in the LogManager constructor
  and cleaning it up in the destructor, the class follows RAII principles. This
  means that the logger setup and teardown are automatically handled as the
  LogManager instance is created and destroyed.

### Thread-Safe Singleton Pattern

The `getInstance` method ensures that only one instance of LogManager is
created, even in a multi-threaded environment:

```{code-block} cpp
LogManager *LogManager::getInstance() {
    LogManager *lm = instance.load(std::memory_order_acquire);
    if (!lm) {
        std::lock_guard<std::mutex> myLock(lmMutex);
        lm = instance.load(std::memory_order_relaxed);
        if (!lm) {
            lm = new LogManager();
            instance.store(lm, std::memory_order_release);
        }
    }
    [[maybe_unused]] volatile bool dummy{}; // Prevent being optimized away
    return lm;
```

- **Memory Ordering:** The use of `std::memory_order_acquire` and
  `std::memory_order_release` ensures that the initialization of the singleton is
  visible to all threads once it's complete.

- **Double-Checked Locking:** The method employs a double-checked locking
  pattern, which first checks if the instance is already created (without
  locking), and only locks if the instance is nullptr. This minimizes the
  locking overhead.

## Logger Setup

The `setup_logger` function configures the logger with the following details:

- **Sinks:** Two sinks are used—one for console output (`stdout_color_sink_mt`) and another for file output (`basic_file_sink_mt`). The file sink is configured to overwrite existing logs.

- **Log Patterns:** The logger is set up with specific patterns to format the log output. For example, the traceback logger includes detailed information such as log level, source file, line number, and function name.

- **Log Levels:** The log level is set to trace to capture all levels of log
  messages.

```{code-block} cpp
// --- Start Logging setup
eonc::LogManager::getInstance();
auto logger = spdlog::get("combi");
// End logging setup
```

- **Accessing the Logger:** After the LogManager is initialized, you can obtain the logger using `spdlog::get("combi")`. This logger can then be used to log messages throughout the application.

- **Automatic Cleanup:** When the application terminates, the  destructor is
  called, automatically cleaning up all loggers and flushing any remaining log
  messages.

### Task-Specific Loggers
- **Base Logger:** The "combi" logger is the default and is used for general
  application logs.

- **Task Loggers:** Tasks like PointJob and ConjugateGradients create their own
  loggers to write results to specific files. For example:
  
#### Examples

- `PointJob` creates a "point" logger that logs energy and force data to
`results.dat`. 
- `ConjugateGradients` creates a "cg" logger for logging
optimization progress to `_cg.log`.

This setup allows each task to log its results independently while maintaining a
shared logging infrastructure.
