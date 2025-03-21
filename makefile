# Compiler
CC = gcc

# Compiler flags
#CFLAGS = -Wall -Wextra -O2 
CFLAGS = -Iinclude

# Executable name
TARGET = my_program

# Source files
SRCS = src/main.c src/tsp.c src/chrono.c src/vns.c src/tsp_greedy.c src/tabu.c

# Object files
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Link the object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

# Compile the source files into object files
src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run the program with command-line arguments
run: $(TARGET)
	./$(TARGET) $(ARGS)

# Clean up the build files
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean run
