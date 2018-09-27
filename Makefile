TARGET	= mps
OBJS	= main.o

OUT_DIR	= out

PLOT_SIZE	= 2.5
D_PLOT		= 1

CC	= gcc

OPT = -O3 -march=native

CFLAGS	= $(OPT)
CXXFLAGS= -std=c++17 $(OPT)
LDFLAGS	= -lm -lstdc++ -lstdc++fs

BENCH= off
BENCH_NUM=10
BENCH_STATUS=$(CC) $(OPT)
BENCH_LIST= -O0,\
			-O1,\
			-O2,\
			-O2 -march=native,\
			-O3,\
			-O3 -march=native,\
			-Ofast,\
			-Ofast -march=native
ifeq ($(BENCH),on)
	CFLAGS += -DBENCH
	CXXFLAGS += -DBENCH
endif

OMP= on
ifeq ($(OMP),on)
	CFLAGS	+= -fopenmp -DOPENMP
	CXXFLAGS+= -fopenmp -DOPENMP
	LDFLAGS	+= -fopenmp
	BENCH_STATUS+=:openmp
endif

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(CXXFLAGS) -c $< -o $@

default: $(TARGET)

clean:
	rm -rf $(TARGET)
	rm -rf *.o

src_only:
	make clean
	rm -rf $(OUT_DIR)

full:
	make clean
	make
full_run:
	make full
	make run

bench:
	make full BENCH=on
	@t=0.0; \
	  for i in {1..$(BENCH_NUM)}; \
	  { ti=`make run | grep "simulation time" | grep -o "[.0-9]*"`; t=`echo "scale=10; $$t+$$ti" | bc`; }; \
	  echo "[`echo $(BENCH_STATUS)`] time:`echo "scale=10; $$t/$(BENCH_NUM)" | bc`"

bench_opt:
	@IFS=',';for opt in `echo $(BENCH_LIST)`; { make -s bench OPT="$$opt"; }

bench_all:
	@make bench_opt -s CC=gcc	OMP=off
	@make bench_opt -s CC=clang	OMP=off
	@make bench_opt -s CC=icc	OMP=off
	@make bench_opt -s CC=gcc	OMP=on
	@make bench_opt -s CC=clang	OMP=on
	@make bench_opt -s CC=icc	OMP=on

run: $(TARGET)
	./$< input/dambreak.prof $(OUT_DIR)

plot:
	gnuplot -c anim.gp $(OUT_DIR) $(D_PLOT) $(PLOT_SIZE)

vtk:
	ls $(OUT_DIR)/*.prof | xargs -i util/prof2vtk.rb "{}"

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)
