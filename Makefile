TARGET	= mps
OBJS	= main.o

CC	= gcc
CXX	= g++

%.o: %.c
	$(CC) -c $< -o $@

%.o: %.cc
	$(CXX) -c $< -o $@

default: $(TARGET)

run: $(TARGET)
	./$< input/dambreak.prof

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS)
