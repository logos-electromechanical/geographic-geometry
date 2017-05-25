# This file uses GNU make extensions.

# Tweakable / overridable compilation options (debugging and optimization)
OPTS = -gdwarf-4 -g3 -Os -Wall
#OPTS = -gdwarf-4 -g3 -Wall

# Paths
GTEST_ROOT= submodules/googletest/googletest
GTEST_DIR= $(GTEST_ROOT)
GTEST_INC= $(GTEST_ROOT)/include/

CXXFLAGS += $(OPTS) -Wno-reorder -g -std=gnu++14 -pthread

CPPFLAGS += -I $(GTEST_INC)
CPPFLAGS += -D _USE_MATH_DEFINES
CPPFLAGS += -MMD -MT $@

LDFLAGS  += -L /usr/local/lib
LDLIBS += -lm
LDLIBS += -lrt

all: test
test: unit_tests
	rm *.log
	./unit_tests
.PHONY: all test

TEST_OBJS = position_test.o
#TEST_OBJS += boundingbox_test.o
TEST_OBJS += point_test.o
TEST_OBJS += twovector_test.o
GTEST_OBJS = gtest.o gtest_main.o
ALL_OBJS+= $(TEST_OBJS) $(GTEST_OBJS)

unit_tests: $(TEST_OBJS) $(GTEST_OBJS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# gtest rules:
# gtest isn't installed as a library by the debian package because it
# depends on compiler flags. Instead, the recommended way to use it is
# to compile its pieces into object files ourselves and then link
# those in.
gtest.o: $(GTEST_DIR) $(GTEST_INC)
	$(CXX) $(CXXFLAGS) -I $(GTEST_INC) -I $(GTEST_DIR) -c -o $@ $(GTEST_DIR)/src/gtest-all.cc
#gtest_main.o: $(GTEST_DIR)/src/gtest_main.cc
gtest_main.o: gtest_main.cc
	$(CXX) $(CXXFLAGS) -I $(GTEST_INC) -c -o $@ $^
