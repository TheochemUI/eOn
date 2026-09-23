/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "eon/fpe_handler.h"

#include <cfenv>
#include <csignal>
#include <cstdint>
#include <cstdio>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <float.h>
#include <windows.h>
#endif

#if defined(__linux__)
#include <ucontext.h>
#endif
#if defined(__APPLE__)
#include <sys/ucontext.h>
#endif

#if defined(__APPLE__) && defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace eonc {

#ifdef _WIN32
// Report each exception class once. Clearing the status alone is not enough
// for a true continue: the faulting op re-executes and re-traps forever.
// After the first report, unmask-trapping is demoted for that class so the
// instruction completes with the IEEE default (Inf/NaN) and the process
// proceeds.
static LONG WINAPI windowsFPEHandler(EXCEPTION_POINTERS *info) {
  DWORD code = info->ExceptionRecord->ExceptionCode;
  static bool reported_div = false;
  static bool reported_inv = false;
  static bool reported_ovf = false;
  static bool reported_other = false;
  switch (code) {
  case EXCEPTION_FLT_DIVIDE_BY_ZERO:
    if (!reported_div) {
      reported_div = true;
      fprintf(stderr, "FPE (continuing, masking further): division by zero\n");
    }
    break;
  case EXCEPTION_FLT_INVALID_OPERATION:
    if (!reported_inv) {
      reported_inv = true;
      fprintf(stderr, "FPE (continuing, masking further): invalid operation\n");
    }
    break;
  case EXCEPTION_FLT_OVERFLOW:
    if (!reported_ovf) {
      reported_ovf = true;
      fprintf(stderr, "FPE (continuing, masking further): overflow\n");
    }
    break;
  case EXCEPTION_FLT_UNDERFLOW:
  case EXCEPTION_FLT_INEXACT_RESULT:
  case EXCEPTION_FLT_DENORMAL_OPERAND:
  case EXCEPTION_FLT_STACK_CHECK:
    if (!reported_other) {
      reported_other = true;
      fprintf(stderr, "FPE (continuing, masking further): other float fault\n");
    }
    break;
  default:
    return EXCEPTION_CONTINUE_SEARCH;
  }
  // Re-mask every class we care about so CONTINUE_EXECUTION does not re-trap.
  _clearfp();
  unsigned int control = 0;
  _controlfp_s(&control, _MCW_EM, _MCW_EM);
  return EXCEPTION_CONTINUE_EXECUTION;
}
#else
// MXCSR exception-mask bits (Intel SDM): bit7=IM, bit8=DM, bit9=ZM,
// bit10=OM, bit11=UM, bit12=PM. Sticky status flags are bits 0-5.
static constexpr unsigned MXCSR_MASK_IM = 1u << 7;
static constexpr unsigned MXCSR_MASK_ZM = 1u << 9;
static constexpr unsigned MXCSR_MASK_OM = 1u << 10;

// One DIV/IDIV. length 0 means the bytes at IP are not that instruction.
// width is the operand size in bytes (1, 2, 4, or 8).
struct DecodedDiv {
  int length;
  int width;
};

// #DE is a fault: the saved IP points at the divide. Only DIV and IDIV
// raise it. Prefixes and a ModRM (plus SIB/displacement) are enough to
// measure that one instruction; anything else is left untouched.
static DecodedDiv decode_div_insn(const unsigned char *code, bool long_mode) {
  DecodedDiv out{0, 0};
  const unsigned char *p = code;
  const unsigned char *limit = code + 15;
  bool operand16 = false;
  bool rex_w = false;
  bool addr16 = false;
  while (p < limit) {
    unsigned char c = *p;
    if (c == 0x66) {
      operand16 = true;
      ++p;
      continue;
    }
    if (c == 0x67) {
      if (!long_mode) {
        addr16 = true;
      }
      ++p;
      continue;
    }
    if (c == 0xF0 || c == 0xF2 || c == 0xF3 || c == 0x26 || c == 0x2E ||
        c == 0x36 || c == 0x3E || c == 0x64 || c == 0x65) {
      ++p;
      continue;
    }
    if (long_mode && c >= 0x40 && c <= 0x4F) {
      rex_w = (c & 0x08) != 0;
      ++p;
      continue;
    }
    break;
  }
  if (addr16 || p >= limit) {
    return out;
  }
  unsigned char opcode = *p++;
  if (opcode != 0xF6 && opcode != 0xF7) {
    return out;
  }
  if (p >= limit) {
    return out;
  }
  unsigned char modrm = *p++;
  unsigned mod = modrm >> 6;
  unsigned reg = (modrm >> 3) & 7u;
  unsigned rm = modrm & 7u;
  if (reg != 6u && reg != 7u) {
    return out;
  }
  if (mod != 3u) {
    bool have_sib = rm == 4u;
    unsigned sib_base = 0;
    if (have_sib) {
      if (p >= limit) {
        return out;
      }
      sib_base = static_cast<unsigned>(*p++ & 7u);
    }
    if ((mod == 0u && rm == 5u) || (have_sib && mod == 0u && sib_base == 5u) ||
        mod == 2u) {
      if (p + 4 > limit) {
        return out;
      }
      p += 4;
    } else if (mod == 1u) {
      if (p + 1 > limit) {
        return out;
      }
      ++p;
    }
  }
  int length = static_cast<int>(p - code);
  if (length < 2 || length > 15) {
    return out;
  }
  out.length = length;
  if (opcode == 0xF6) {
    out.width = 1;
  } else if (rex_w) {
    out.width = 8;
  } else if (operand16) {
    out.width = 2;
  } else {
    out.width = 4;
  }
  return out;
}

struct X86DivRegs {
  uintptr_t *ip;
  uintptr_t *ax;
  uintptr_t *dx;
};

static X86DivRegs x86_div_regs(void *scp) {
  X86DivRegs regs{nullptr, nullptr, nullptr};
  if (scp == nullptr) {
    return regs;
  }
  auto *ctx = static_cast<ucontext_t *>(scp);
#if defined(__linux__) && defined(__x86_64__)
  regs.ip = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_RIP]);
  regs.ax = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_RAX]);
  regs.dx = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_RDX]);
#elif defined(__linux__) && defined(__i386__)
  regs.ip = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_EIP]);
  regs.ax = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_EAX]);
  regs.dx = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext.gregs[REG_EDX]);
#elif defined(__APPLE__) && defined(__x86_64__)
  if (ctx->uc_mcontext != nullptr) {
    regs.ip = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__rip);
    regs.ax = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__rax);
    regs.dx = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__rdx);
  }
#elif defined(__APPLE__) && defined(__i386__)
  if (ctx->uc_mcontext != nullptr) {
    regs.ip = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__eip);
    regs.ax = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__eax);
    regs.dx = reinterpret_cast<uintptr_t *>(&ctx->uc_mcontext->__ss.__edx);
  }
#else
  (void)ctx;
#endif
  return regs;
}

// Advance past a faulting DIV/IDIV and clear its quotient. Returns false
// when IP does not point at a divide (INTO's trap IP is already next).
static bool advance_integer_div(void *scp) {
  X86DivRegs regs = x86_div_regs(scp);
  if (regs.ip == nullptr || regs.ax == nullptr || regs.dx == nullptr) {
    return false;
  }
  const auto *code = reinterpret_cast<const unsigned char *>(*regs.ip);
  DecodedDiv div = decode_div_insn(code, sizeof(void *) == 8);
  if (div.length == 0) {
    return false;
  }
  if (div.width <= 2) {
    *regs.ax &= ~uintptr_t{0xFFFFu};
    if (div.width == 2) {
      *regs.dx &= ~uintptr_t{0xFFFFu};
    }
  } else {
    *regs.ax = 0;
    *regs.dx = 0;
  }
  *regs.ip += static_cast<uintptr_t>(div.length);
  return true;
}

static void write_fault_rip(void *scp) {
#if defined(__linux__) && defined(__x86_64__)
  auto *ctx_log = static_cast<ucontext_t *>(scp);
  unsigned long rip =
      static_cast<unsigned long>(ctx_log->uc_mcontext.gregs[REG_RIP]);
  char hex[] = "FPE rip=0x0000000000000000\n";
  for (int i = 0; i < 16; ++i) {
    unsigned nibble = static_cast<unsigned>((rip >> (4 * (15 - i))) & 0xFu);
    hex[10 + i] =
        static_cast<char>(nibble < 10 ? '0' + nibble : 'a' + (nibble - 10));
  }
  write(STDERR_FILENO, hex, sizeof(hex) - 1);
#else
  (void)scp;
#endif
}

static void fpe_signal_handler(int sig, siginfo_t *sip, void *scp) {
  // Async-signal-safe only: write(2) and sig_atomic_t. No iostream, malloc,
  // backtrace, or fenv helpers (fedisableexcept / feclearexcept are not
  // async-signal-safe). All continue-state is written into the saved ucontext
  // so it is restored on sigreturn.
  //
  // x86 cannot "continue" past a trapped FP op by clearing sticky flags:
  // flags are bits 0-5 of MXCSR/swd, but the exception MASK bits live at
  // MXCSR 7-12. Clearing 0x3F leaves trapping armed, so the faulting
  // instruction re-executes on the same operands and re-raises forever
  // (report, sigreturn, refault) -- multi-GB identical stderr lines and a
  // client stuck at ~100% CPU. Mask the class in the restored MXCSR so
  // re-execution produces the IEEE default (Inf/NaN) and proceeds.
  static volatile sig_atomic_t reported_div = 0;
  static volatile sig_atomic_t reported_inv = 0;
  static volatile sig_atomic_t reported_ovf = 0;
  static volatile sig_atomic_t reported_unk = 0;
  static volatile sig_atomic_t reported_int = 0;

  static constexpr char prefix[] = "FPE (continuing, masking further): ";
  static constexpr char msg_div[] = "division by zero\n";
  static constexpr char msg_inv[] = "invalid operation\n";
  static constexpr char msg_ovf[] = "overflow\n";
  static constexpr char msg_unk[] = "unknown\n";

  // Default: mask all three classes we enable at startup, so an unknown
  // si_code cannot leave trapping armed and re-storm.
  unsigned mxcsr_mask_bits = MXCSR_MASK_IM | MXCSR_MASK_ZM | MXCSR_MASK_OM;
  volatile sig_atomic_t *reported = &reported_unk;
  const char *msg = msg_unk;
  size_t msg_len = sizeof(msg_unk) - 1;

  // Integer #DE has no mask bit. Returning here re-executes the divide.
  if (sip->si_code == FPE_INTDIV || sip->si_code == FPE_INTOVF) {
    if (reported_int == 0) {
      reported_int = 1;
      static constexpr char iprefix[] = "FPE (continuing, skipping divide): ";
      static constexpr char msg_idiv[] = "integer divide\n";
      static constexpr char msg_iovf[] = "integer overflow\n";
      const char *imsg = sip->si_code == FPE_INTDIV ? msg_idiv : msg_iovf;
      size_t imsg_len = sizeof(msg_iovf) - 1;
      if (sip->si_code == FPE_INTDIV) {
        imsg_len = sizeof(msg_idiv) - 1;
      }
      write(STDERR_FILENO, iprefix, sizeof(iprefix) - 1);
      write(STDERR_FILENO, imsg, imsg_len);
      write_fault_rip(scp);
    }
    if (!advance_integer_div(scp) && sip->si_code == FPE_INTDIV) {
      static constexpr char stuck[] =
          "FPE integer divide: could not skip faulting instruction\n";
      write(STDERR_FILENO, stuck, sizeof(stuck) - 1);
      _exit(128 + SIGFPE);
    }
    return;
  }

  switch (sip->si_code) {
  case FPE_FLTDIV:
    reported = &reported_div;
    msg = msg_div;
    msg_len = sizeof(msg_div) - 1;
    mxcsr_mask_bits = MXCSR_MASK_ZM;
    break;
  case FPE_FLTINV:
    reported = &reported_inv;
    msg = msg_inv;
    msg_len = sizeof(msg_inv) - 1;
    mxcsr_mask_bits = MXCSR_MASK_IM;
    break;
  case FPE_FLTOVF:
    reported = &reported_ovf;
    msg = msg_ovf;
    msg_len = sizeof(msg_ovf) - 1;
    mxcsr_mask_bits = MXCSR_MASK_OM;
    break;
  default:
    break;
  }

  if (*reported == 0) {
    *reported = 1;
    write(STDERR_FILENO, prefix, sizeof(prefix) - 1);
    write(STDERR_FILENO, msg, msg_len);
    // First-fault RIP for post-mortem addr2line / offline diagnosis.
    write_fault_rip(scp);
  }

#if defined(__linux__) && (defined(__x86_64__) || defined(__i386__))
  ucontext_t *ctx = static_cast<ucontext_t *>(scp);
  if (ctx->uc_mcontext.fpregs) {
    // Clear sticky exception FLAGS (bits 0-5) and arm the MASK bit(s) for
    // the fault class (bits 7-12). Mask sticks after sigreturn because the
    // restored MXCSR becomes the live CPU state.
    ctx->uc_mcontext.fpregs->swd &= ~0x3Fu;
    ctx->uc_mcontext.fpregs->mxcsr &= ~0x3Fu;
    ctx->uc_mcontext.fpregs->mxcsr |= mxcsr_mask_bits;
    // x87 control word: mask bits are 0-5 of cwd (IM, DM, ZM, OM, UM, PM).
    // Set the matching masks so a legacy x87 fault cannot re-storm either.
    if (mxcsr_mask_bits & MXCSR_MASK_ZM) {
      ctx->uc_mcontext.fpregs->cwd |= (1u << 2); // x87 ZM
    }
    if (mxcsr_mask_bits & MXCSR_MASK_IM) {
      ctx->uc_mcontext.fpregs->cwd |= (1u << 0); // x87 IM
    }
    if (mxcsr_mask_bits & MXCSR_MASK_OM) {
      ctx->uc_mcontext.fpregs->cwd |= (1u << 3); // x87 OM
    }
  }
#elif defined(__linux__) && defined(__aarch64__)
  // ARM polarity is the opposite of MXCSR: FPCR trap-enable bits SET mean
  // trap. Clearing sticky FPSR flags alone re-executes with trapping still
  // armed (eOn-qkuw / eOn-g7fl on aarch64). Clear the matching FPCR enables.
  constexpr uint32_t kFpsimdMagic = 0x46508001u;
  constexpr unsigned kFpcrIoe = 1u << 8;
  constexpr unsigned kFpcrDze = 1u << 9;
  constexpr unsigned kFpcrOfe = 1u << 10;
  unsigned fpcr_clear = kFpcrIoe | kFpcrDze | kFpcrOfe;
  if (mxcsr_mask_bits == MXCSR_MASK_ZM) {
    fpcr_clear = kFpcrDze;
  } else if (mxcsr_mask_bits == MXCSR_MASK_IM) {
    fpcr_clear = kFpcrIoe;
  } else if (mxcsr_mask_bits == MXCSR_MASK_OM) {
    fpcr_clear = kFpcrOfe;
  }
  auto *ctx = static_cast<ucontext_t *>(scp);
  unsigned char *p =
      reinterpret_cast<unsigned char *>(ctx->uc_mcontext.__reserved);
  unsigned char *end = p + sizeof(ctx->uc_mcontext.__reserved);
  struct A64Head {
    uint32_t magic;
    uint32_t size;
  };
  struct Fpsimd {
    A64Head head;
    uint32_t fpsr;
    uint32_t fpcr;
  };
  while (p + sizeof(A64Head) <= end) {
    auto *h = reinterpret_cast<A64Head *>(p);
    if (h->magic == 0 || h->size < sizeof(A64Head)) {
      break;
    }
    if (h->magic == kFpsimdMagic && h->size >= sizeof(Fpsimd) &&
        p + h->size <= end) {
      auto *f = reinterpret_cast<Fpsimd *>(p);
      f->fpsr &= ~0x1Fu;
      f->fpcr &= ~fpcr_clear;
      break;
    }
    if (h->size == 0) {
      break;
    }
    p += h->size;
  }
#elif defined(__APPLE__) && defined(__x86_64__)
  // Darwin restores SSE state from uc_mcontext->__fs. Mask bits have the
  // same polarity as Linux MXCSR: set means the class does not trap. Without
  // that update the faulting instruction re-executes and re-raises forever.
  auto *ctx = static_cast<ucontext_t *>(scp);
  if (ctx->uc_mcontext) {
    auto &fs = ctx->uc_mcontext->__fs;
    fs.__fpu_mxcsr &= ~0x3Fu;
    fs.__fpu_mxcsr |= mxcsr_mask_bits;
    if (mxcsr_mask_bits & MXCSR_MASK_ZM) {
      fs.__fpu_fcw.__zdiv = 1;
    }
    if (mxcsr_mask_bits & MXCSR_MASK_IM) {
      fs.__fpu_fcw.__invalid = 1;
    }
    if (mxcsr_mask_bits & MXCSR_MASK_OM) {
      fs.__fpu_fcw.__ovrfl = 1;
    }
    fs.__fpu_fsw.__invalid = 0;
    fs.__fpu_fsw.__denorm = 0;
    fs.__fpu_fsw.__zdiv = 0;
    fs.__fpu_fsw.__ovrfl = 0;
    fs.__fpu_fsw.__undfl = 0;
    fs.__fpu_fsw.__precis = 0;
  }
#elif defined(__APPLE__) && defined(__aarch64__)
  constexpr unsigned kFpcrIoe = 1u << 8;
  constexpr unsigned kFpcrDze = 1u << 9;
  constexpr unsigned kFpcrOfe = 1u << 10;
  unsigned fpcr_clear = kFpcrIoe | kFpcrDze | kFpcrOfe;
  if (mxcsr_mask_bits == MXCSR_MASK_ZM) {
    fpcr_clear = kFpcrDze;
  } else if (mxcsr_mask_bits == MXCSR_MASK_IM) {
    fpcr_clear = kFpcrIoe;
  } else if (mxcsr_mask_bits == MXCSR_MASK_OM) {
    fpcr_clear = kFpcrOfe;
  }
  auto *ctx = static_cast<ucontext_t *>(scp);
  if (ctx->uc_mcontext) {
    ctx->uc_mcontext->__ns.__fpsr &= ~0x1Fu;
    ctx->uc_mcontext->__ns.__fpcr &= ~fpcr_clear;
  }
#endif
  (void)sig;
}
#endif

void enableFPE() {
#ifdef _WIN32
  // Register Windows SEH handler for FPE reporting
  SetUnhandledExceptionFilter(windowsFPEHandler);
  // Enable floating-point exceptions on Windows
  _controlfp_s(nullptr, 0, _MCW_EM);
  _controlfp_s(nullptr, ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW),
               _MCW_EM);
#elif defined(__unix__)
  // Enable floating-point exceptions on Unix
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif defined(__APPLE__) && defined(__aarch64__)
  // ARM: trap-enable bits live in FPCR (IOE/DZE/OFE), not FPSR flags.
  fenv_t env;
  fegetenv(&env);
  env.__fpcr |= ((1u << 8) | (1u << 9) | (1u << 10));
  fesetenv(&env);
#elif defined(__APPLE__) && defined(__x86_64__)
  // Enable floating-point exceptions on Intel macOS
  _MM_SET_EXCEPTION_MASK(
      _MM_MASK_MASK &
      ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW));
#else
  fprintf(stderr, "FPE trapping not supported on this platform.\n");
#endif

#ifndef _WIN32
  // Register POSIX signal handler
  struct sigaction act;
  act.sa_sigaction = fpe_signal_handler;
  sigemptyset(&act.sa_mask);
  act.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &act, nullptr);
#endif
}

void disableFPE() {
#ifdef _WIN32
  // Mask all floating-point exceptions (restore default behavior)
  unsigned int control;
  _controlfp_s(&control, _MCW_EM, _MCW_EM);
#elif defined(__unix__)
  fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif defined(__APPLE__)
  fenv_t env;
  fegetenv(&env);
#if defined(__aarch64__)
  env.__fpcr &= ~((1u << 8) | (1u << 9) | (1u << 10));
#elif defined(__x86_64__)
  // enableFPE clears MXCSR IM/ZM/OM. Restoring that environment unchanged
  // leaves the traps armed.
  env.__mxcsr |= (MXCSR_MASK_IM | MXCSR_MASK_ZM | MXCSR_MASK_OM);
  env.__mxcsr &= ~0x3Fu;
  env.__control = static_cast<unsigned short>(env.__control | (1u << 0) |
                                              (1u << 2) | (1u << 3));
  env.__status = static_cast<unsigned short>(env.__status & ~0x3Fu);
#endif
  fesetenv(&env);
#if defined(__x86_64__)
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | _MM_MASK_INVALID |
                         _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW);
#endif
#endif
}

void FPEHandler::eat_fpe() {
  std::lock_guard<std::mutex> lock(mutex_);
  feholdexcept(&orig_feenv);
}

void FPEHandler::restore_fpe() {
  std::lock_guard<std::mutex> lock(mutex_);
  fesetenv(&orig_feenv);
}

} // namespace eonc
