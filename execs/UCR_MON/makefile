SOURCES := $(wildcard *.cpp)
TARGETS := $(patsubst %.cpp, %, $(SOURCES))

all:$(TARGETS)

% :: %.cpp
	g++ $< -O3 -o $@
