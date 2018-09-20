TARGET	= mps
OBJS	= main.o

OUT_DIR	= out

CC	= gcc
CXX	= g++

CFLAGS	=
CXXFLAGS= -std=c++17
LDFLAGS	= -lstdc++fs

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

default: $(TARGET)

run: $(TARGET)
	./$< input/dambreak.prof $(OUT_DIR)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)
