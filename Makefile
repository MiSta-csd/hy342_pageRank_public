#    CS342 - Programming Assignment1, 2023
#    Author: Minos Stavrakakis - csd4120

# the compiler to use
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -pthread
  
#files to link:
LFLAGS = -lm -latomic

# Source files
SRCS = main.c 

# target file
TARGET = cs342_ass1.exec
  
all: $(TARGET)
  
$(TARGET): $(SRCS)
	@$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LFLAGS)

clean:
	@rm $(TARGET)
