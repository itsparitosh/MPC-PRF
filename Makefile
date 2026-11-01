CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -march=native -Wall -Wextra -Wno-ignored-attributes -Wno-unused-result
INCLUDES  = -Iinclude
# libsodium: a local build in ./deps (tools/build_libsodium.sh, no sudo) if there is
# one, otherwise the system's (pkg-config also finds Homebrew installs on macOS).
DEPS := $(CURDIR)/deps
ifneq ($(wildcard $(DEPS)/lib/libsodium.a),)
SODIUM_CFLAGS := -I$(DEPS)/include
SODIUM_LIBS   := $(DEPS)/lib/libsodium.a
else
SODIUM_CFLAGS := $(shell pkg-config --cflags libsodium 2>/dev/null)
SODIUM_LIBS   := $(shell pkg-config --libs libsodium 2>/dev/null || echo -lsodium)
endif
INCLUDES += $(SODIUM_CFLAGS)
LIBS      = -lgmpxx -lgmp -lcrypto $(SODIUM_LIBS) -pthread

# A clear message instead of "sodium.h: No such file or directory"
ifeq ($(filter deps clean,$(MAKECMDGOALS)),)
HAVE_SODIUM := $(shell printf '\043include <sodium.h>\nint main(){return sodium_init();}' | $(CXX) -x c++ $(SODIUM_CFLAGS) - -x none -o /dev/null $(SODIUM_LIBS) 2>/dev/null && echo yes)
ifneq ($(HAVE_SODIUM),yes)
$(error libsodium is missing. Without sudo: run "make deps" (builds it inside this folder). With sudo: Ubuntu/Debian "sudo apt install libsodium-dev", Fedora "sudo dnf install libsodium-devel", macOS "brew install libsodium pkg-config")
endif
endif
HEADERS   = $(wildcard include/*/*.h)

all: mpcprf test_prf

mpcprf: src/mpcprf.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

test_prf: tests/test_prf.cpp tests/harness.h $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Itests $< -o $@ $(LIBS)

# unit tests (each test runs S0 and S1 as two processes over a socketpair)
test: test_prf
	./test_prf

# the real program: one honest run, then one cheating server per case
run-tests: mpcprf
	bash tests/run_tests.sh

clean:
	rm -f mpcprf test_prf
	rm -rf test_run

# builds libsodium into ./deps (no sudo needed)
deps:
	bash tools/build_libsodium.sh

.PHONY: deps all test run-tests clean
