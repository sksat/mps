TARGET	= mps
OBJS	= main.o

OUT_DIR	= out

PLOT_SIZE	= 2.5
D_PLOT		= 1

CC	= gcc
CXX	= g++

CFLAGS	=
CXXFLAGS= -std=c++17 -O3 -fopenmp
LDFLAGS	= -lstdc++fs -fopenmp

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

default: $(TARGET)

run: $(TARGET)
	./$< input/dambreak.prof $(OUT_DIR)

plot:
	gnuplot -c anim.gp $(OUT_DIR) $(D_PLOT) $(PLOT_SIZE)

vtk:
	ls $(OUT_DIR)/*.prof | xargs -i util/prof2vtk.rb "{}"

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)
