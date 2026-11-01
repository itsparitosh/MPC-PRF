#!/bin/bash
# =====================================================================
# build_libsodium.sh -- builds libsodium inside this folder (no sudo needed)
#
# Uses the official release third_party/libsodium-1.0.20.tar.gz (from
# github.com/jedisct1/libsodium, release 1.0.20-RELEASE), checks its SHA-256,
# builds a static library and installs it into ./deps. The Makefile uses
# ./deps automatically when it exists. Takes one or two minutes.
# =====================================================================
set -e
cd "$(dirname "$0")/.."
ROOT=$(pwd)
TAR=third_party/libsodium-1.0.20.tar.gz
SHA=ebb65ef6ca439333c2bb41a0c1990587288da07f6c7fd07cb3a18cc18d30ce19
if command -v sha256sum >/dev/null; then GOT=$(sha256sum "$TAR" | cut -d' ' -f1); else GOT=$(shasum -a 256 "$TAR" | cut -d' ' -f1); fi
if [ "$GOT" != "$SHA" ]; then echo "SHA-256 of $TAR does not match the official release; stopping."; exit 1; fi
rm -rf deps/build && mkdir -p deps/build
tar -xzf "$TAR" -C deps/build
cd deps/build/libsodium-1.0.20
./configure --prefix="$ROOT/deps" --disable-shared --enable-static --with-pic > ../configure.log
make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)" > ../make.log
make install > ../install.log
cd "$ROOT" && rm -rf deps/build
echo "libsodium installed in $ROOT/deps; now run make"
