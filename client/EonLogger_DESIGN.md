# EonLogger.h - Modern C++20 Logging Made Elegant

## Problem Statement (Graeme's Complaint)

The old logging setup was verbose, repetitive, and annoying:

**Before (spdlog/old quill setup):**
```cpp
class ConjugateGradients {
  quill::Logger *m_log{nullptr};
  
  ConjugateGradients(...) {
    // 14 lines of boilerplate!
    m_log = quill::Frontend::create_or_get_logger(
        "cg",
        quill::Frontend::create_or_get_sink<quill::FileSink>(
            "_cg.log",
            []() {
              quill::FileSinkConfig cfg;
              cfg.set_open_mode('w');
              return cfg;
            }(),
            quill::FileEventNotifier{}),
        quill::PatternFormatterOptions{
            quill::PatternFormatterOptions{"%(message)"}},
        quill::ClockSourceType::System);
  }
};
```

For simple classes using the default logger:
```cpp
class Matter {
  quill::Logger *m_log{nullptr};
  
  Matter(...) {
    m_log = quill::Frontend::get_logger("combi");  // Still manual
  }
};
```

## Solution: C++20 Elegant Logging

**After (with EonLogger.h):**

### 1. Simple Case - Default Logger
```cpp
#include "EonLogger.h"

class Matter {
  eonc::log::Scoped m_log;  // That's it! Auto-initializes to "combi"
  
  // No constructor code needed!
  
  void someMethod() {
    LOG_INFO(m_log, "positions updated");
  }
};
```

### 2. File Logger Case
```cpp
#include "EonLogger.h"

class ConjugateGradients {
  eonc::log::FileScoped m_log{"cg", "_cg.log"};  // One line!
  
  // No constructor code needed!
  
  void step() {
    LOG_DEBUG(m_log, "iteration {}", m_cg_i);
  }
};
```

### 3. One-Shot Logging (No Member Variable)
```cpp
#include "EonLogger.h"

void quick_function() {
  EONC_LOG_INFO("this is a quick log");
  EONC_LOG_DEBUG("value: {}", some_value);
  // No logger variable needed!
}
```

## Key Features

### 1. **Zero Boilerplate for Default Logger**
```cpp
// Old way (2 lines + constructor init)
quill::Logger *m_log{nullptr};
// In constructor: m_log = quill::Frontend::get_logger("combi");

// New way (1 line, auto-init)
eonc::log::Scoped m_log;
```

### 2. **Implicit Conversion**
The `Scoped` and `FileScoped` types implicitly convert to `quill::Logger*`, so all existing `LOG_*` macros work without changes:

```cpp
eonc::log::Scoped m_log;
LOG_INFO(m_log, "message");  // Works! No .get() or ->ptr() needed
```

### 3. **C++20 Features Used**
- `[[nodiscard]]` - Prevent accidental logger discard
- `std::string_view` - Zero-copy string parameters
- `std::source_location` - Ready for future source-aware logging
- `noexcept` - Clear exception guarantees
- Implicit conversions with `operator Logger*()`

### 4. **Type Safety**
```cpp
// These both work:
eonc::log::Scoped log1;
auto* log2 = eonc::log::get();

// But Scoped is safer (RAII, clear ownership)
```

## Migration Guide

### Step 1: Replace Manual Initialization

**Before:**
```cpp
class MyClass {
  quill::Logger *m_log{nullptr};
public:
  MyClass() {
    m_log = quill::Frontend::get_logger("combi");
  }
};
```

**After:**
```cpp
#include "EonLogger.h"

class MyClass {
  eonc::log::Scoped m_log;  // Remove constructor line!
};
```

### Step 2: Replace File Logger Setup

**Before:**
```cpp
quill::Logger *m_log{nullptr};

Constructor() {
  m_log = quill::Frontend::create_or_get_logger(
      "name",
      quill::Frontend::create_or_get_sink<quill::FileSink>(
          "file.log",
          []() {
            quill::FileSinkConfig cfg;
            cfg.set_open_mode('w');
            return cfg;
          }(),
          quill::FileEventNotifier{}),
      quill::PatternFormatterOptions{"%(message)"},
      quill::ClockSourceType::System);
}
```

**After:**
```cpp
#include "EonLogger.h"

eonc::log::FileScoped m_log{"name", "file.log"};
// Remove entire constructor block!
```

### Step 3: One-Shot Logs

**Before:**
```cpp
void helper_function() {
  auto* log = quill::Frontend::get_logger("combi");
  LOG_INFO(log, "message");
}
```

**After:**
```cpp
#include "EonLogger.h"

void helper_function() {
  EONC_LOG_INFO("message");  // No logger variable!
}
```

## Benefits

1. **Less Code**: 90% reduction in logger setup boilerplate
2. **Clearer Intent**: `eonc::log::Scoped m_log;` says exactly what it does
3. **Safer**: RAII, no manual nullptr checks, no forgotten initialization
4. **Modern**: Uses C++20 features throughout
5. **Compatible**: Drop-in replacement, all existing LOG_* macros work
6. **Documented**: Each function/class has clear usage examples

## Example: Full Class Refactor

**Before (98 lines with logging setup):**
```cpp
class ConjugateGradients : public Optimizer {
  quill::Logger *m_log{nullptr};
  
public:
  ConjugateGradients(std::shared_ptr<ObjectiveFunction> a_objf,
                     const Parameters &a_params)
      : Optimizer(a_objf, OptType::CG, a_params),
        m_directionOld{(a_objf->getPositions()).setZero()},
        m_forceOld{(a_objf->getPositions()).setZero()},
        m_cg_i{0} {
    m_log = quill::Frontend::create_or_get_logger(
        "cg",
        quill::Frontend::create_or_get_sink<quill::FileSink>(
            "_cg.log",
            []() {
              quill::FileSinkConfig cfg;
              cfg.set_open_mode('w');
              return cfg;
            }(),
            quill::FileEventNotifier{}),
        quill::PatternFormatterOptions{
            quill::PatternFormatterOptions{"%(message)"}},
        quill::ClockSourceType::System);
  }
  
  void step() {
    LOG_DEBUG(m_log, "CG step {}", m_cg_i);
  }
};
```

**After (85 lines, logger setup is 1 line):**
```cpp
#include "EonLogger.h"

class ConjugateGradients : public Optimizer {
  eonc::log::FileScoped m_log{"cg", "_cg.log"};  // ONE LINE!
  
public:
  ConjugateGradients(std::shared_ptr<ObjectiveFunction> a_objf,
                     const Parameters &a_params)
      : Optimizer(a_objf, OptType::CG, a_params),
        m_directionOld{(a_objf->getPositions()).setZero()},
        m_forceOld{(a_objf->getPositions()).setZero()},
        m_cg_i{0} {
    // No logger setup code!
  }
  
  void step() {
    LOG_DEBUG(m_log, "CG step {}", m_cg_i);
  }
};
```

**Result: 13 lines removed, zero functionality lost.**

## Graeme Was Right (Sort Of)

The old setup *was* ugly and provided no value **for the boilerplate**.  
The logs themselves are valuable (you use them all the time).  
We've now made the setup as elegant as the logs are useful.

One line. No manual init. No nullptr. Just works.
