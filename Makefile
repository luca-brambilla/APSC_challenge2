# Luca Brambilla
# 10510718 - 919812

PACS_ROOT		?= /home/luca/Documents/pacs-examples/Examples

CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations

LDFLAGS ?=
LDLIBS  ?= 

# get all files *.cpp
SRCS=$(wildcard *.cpp)
# get the corresponding object file
OBJS = $(SRCS:.cpp=.o)
# get all headers in the working directory
HEADERS=$(wildcard *.hpp)
#
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)

# CPPFLAGS += -I./include -I$(mkBoostInc)
# LDFLAGS  += -L$(mkBoostLib)
# LDLIBS	 += -l boost_iostreams -l boost_system -l boost_filesystem

all: clean $(EXEC)

# COMPILER
%.o: %.cpp %.hpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# LINKER and execute
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
	./$(EXEC)

clean:
	$(RM) *.o
	
distclean: clean
	$(RM) $(EXEC)
	$(RM) *~