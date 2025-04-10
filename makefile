# Object and source structure
SRCDIR = src
INCLUDEDIR = include
OBJS = $(SRCDIR)/main.o $(SRCDIR)/tsp.o $(SRCDIR)/tsp_heuristics.o $(SRCDIR)/tsp_greedy.o $(SRCDIR)/tsp_cplex.o $(SRCDIR)/utils.o
HEADERS = $(INCLUDEDIR)/main.h $(INCLUDEDIR)/tsp.h $(INCLUDEDIR)/tsp_heuristics.h $(INCLUDEDIR)/tsp_greedy.h $(INCLUDEDIR)/tsp_cplex.h $(INCLUDEDIR)/utils.h

EXE = tsp_solver
setting = -1
OS := $(shell uname)

ifndef CPLEX_HOME
    $(warning CPLEX_HOME is not set. Please enter the path to your CPLEX installation:)
    $(shell read CPLEX_HOME_INPUT; export CPLEX_HOME=$$CPLEX_HOME_INPUT)
    CPLEX_HOME := $(shell echo $$CPLEX_HOME) # Update the makefile variable
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
CFLAGS = -Iinclude -Wall -O3 $(INC)
RM = rm -f

.SUFFIXES: .o .c .cpp

$(SRCDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

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
