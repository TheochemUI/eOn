#!/usr/bin/env bash
# Install native deps for pyeonclient manylinux / sdist builds (quill, capnp SHARED+PIC).
# Cap'n Proto must be SHARED + PIC: static .a without -fPIC cannot link into libptlrpc.so.
set -euo pipefail
PREFIX="${PYEONCLIENT_DEPS_PREFIX:-/usr/local}"
# manylinux images run as root; GHA ubuntu runners need sudo for /usr/local.
if [[ "$(id -u)" -eq 0 ]]; then
  SUDO=()
  LDCONFIG=(ldconfig)
else
  SUDO=(sudo)
  LDCONFIG=(sudo ldconfig)
fi
export PATH="${PREFIX}/bin:${PATH:-}"
export LD_LIBRARY_PATH="${PREFIX}/lib64:${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib64/pkgconfig:${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CMAKE_PREFIX_PATH="${PREFIX}${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"

if command -v yum >/dev/null 2>&1; then
  yum install -y eigen3-devel gcc-gfortran patchelf cmake zlib-devel 2>/dev/null || true
elif command -v apt-get >/dev/null 2>&1; then
  "${SUDO[@]}" apt-get update -qq
  "${SUDO[@]}" apt-get install -y --no-install-recommends \
    libeigen3-dev gfortran ninja-build cmake pkg-config zlib1g-dev patchelf \
    build-essential 2>/dev/null || true
fi

export CARGO_HOME="${CARGO_HOME:-${HOME}/.cargo}"
export PATH="${CARGO_HOME}/bin:${PATH}"
if ! command -v rustc >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
  # shellcheck disable=SC1091
  source "${CARGO_HOME}/env"
fi
command -v cbindgen >/dev/null 2>&1 || cargo install cbindgen

# capnproto.org has dropped the TLS handshake mid-build (curl exit 35).
# Retry, then take the same tag from GitHub. The GitHub archive wraps the
# C++ tree in c++/; the release tarball is already that tree.
fetch_tarball() {
  local dest="$1"
  shift
  local url attempt
  for url in "$@"; do
    for attempt in 1 2 3; do
      if curl -fsSL --retry 5 --retry-delay 2 --retry-connrefused -o "${dest}" "${url}"; then
        return 0
      fi
      echo "download failed (${url}, attempt ${attempt})" >&2
      rm -f "${dest}"
      sleep $((attempt * 2))
    done
  done
  echo "could not download ${dest}" >&2
  return 1
}

if ! pkg-config --exists capnp 2>/dev/null; then
  rm -rf /tmp/capnp-src /tmp/capnp-build
  fetch_tarball /tmp/capnp.tgz \
    https://github.com/capnproto/capnproto/archive/refs/tags/v1.0.2.tar.gz \
    https://capnproto.org/capnproto-c++-1.0.2.tar.gz
  mkdir -p /tmp/capnp-src
  tar -xzf /tmp/capnp.tgz -C /tmp/capnp-src --strip-components=1
  capnp_src=/tmp/capnp-src
  if [[ -f /tmp/capnp-src/c++/CMakeLists.txt ]]; then
    capnp_src=/tmp/capnp-src/c++
  fi
  cmake -S "${capnp_src}" -B /tmp/capnp-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_TESTING=OFF
  cmake --build /tmp/capnp-build -j"$(nproc 2>/dev/null || echo 2)"
  "${SUDO[@]}" cmake --install /tmp/capnp-build
  "${LDCONFIG[@]}" 2>/dev/null || true
fi
# Ensure capnp CLI can load shared libs (libcapnpc.so)
export LD_LIBRARY_PATH="${PREFIX}/lib64:${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
"${LDCONFIG[@]}" 2>/dev/null || true
command -v capnp
capnp --version || true
pkg-config --modversion capnp

if ! pkg-config --exists quill 2>/dev/null; then
  rm -rf /tmp/quill-src /tmp/quill-build
  fetch_tarball /tmp/quill.tgz \
    https://github.com/odygrd/quill/archive/refs/tags/v11.0.2.tar.gz
  mkdir -p /tmp/quill-src
  tar -xzf /tmp/quill.tgz -C /tmp/quill-src --strip-components=1
  cmake -S /tmp/quill-src -B /tmp/quill-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DQUILL_BUILD_TESTS=OFF \
    -DQUILL_BUILD_EXAMPLES=OFF
  cmake --build /tmp/quill-build -j"$(nproc 2>/dev/null || echo 2)"
  "${SUDO[@]}" cmake --install /tmp/quill-build
  # Upstream quill often installs CMake config only; meson prefers pkg-config.
  if ! pkg-config --exists quill 2>/dev/null; then
    for d in "${PREFIX}/lib64/pkgconfig" "${PREFIX}/lib/pkgconfig"; do
      "${SUDO[@]}" mkdir -p "$d"
    done
    pc="${PREFIX}/lib/pkgconfig/quill.pc"
    if [[ ! -d "${PREFIX}/lib/pkgconfig" && -d "${PREFIX}/lib64/pkgconfig" ]]; then
      pc="${PREFIX}/lib64/pkgconfig/quill.pc"
    fi
    "${SUDO[@]}" tee "$pc" >/dev/null <<EOF
prefix=${PREFIX}
includedir=\${prefix}/include
Name: quill
Description: Asynchronous Low Latency C++ Logging Library
Version: 11.0.2
Cflags: -I\${includedir}
EOF
  fi
fi
export PKG_CONFIG_PATH="${PREFIX}/lib64/pkgconfig:${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
pkg-config --modversion quill
echo "pyeonclient_manylinux_deps: OK"
