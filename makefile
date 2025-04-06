# Object and source structure
SRCDIR = src
INCLUDEDIR = src
OBJS = $(SRCDIR)/main.o $(SRCDIR)/tsp.o $(SRCDIR)/tsp_heuristics.o $(SRCDIR)/tsp_greedy.o $(SRCDIR)/tsp_cplex.o $(SRCDIR)/utils.o
HEADERS = $(INCLUDEDIR)/main.o $(INCLUDEDIR)/tsp.o $(INCLUDEDIR)/tsp_heuristics.o $(INCLUDEDIR)/tsp_greedy.o $(INCLUDEDIR)/tsp_cplex.o $(INCLUDEDIR)/utils.o

EXE = tsp_solver
setting = -1
OS := $(shell uname)

# ---------------------------------------------------------------------
# CPLEX dynamic setup
# ---------------------------------------------------------------------
CPLEX_HOME ?= $(shell echo $$CPLEX_HOME)

ifndef CPLEX_HOME
$(error CPLEX_HOME is not set. Please export it: export CPLEX_HOME=/path/to/cplex)
endif

ifeq ($(OS),Linux)
	setting = 1
	CC = gcc
	AR = ar rc
	LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -ldl
	INC = -I. -Iinclude -I${CPLEX_HOME}/include/ilcplex
endif

# ---------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------
CFLAGS = -Iinclude -Wall -O3
RM = rm -f

.SUFFIXES: .o .c .cpp

$(SRCDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

clean:
	$(RM) $(OBJS) $(EXE)

again:
	make clean
	make

who:
	@echo "User: $(USER)"
	@echo "OS: $(OS) (setting = $(setting))"
	@echo "CPLEX_HOME: $(CPLEX_HOME)"
