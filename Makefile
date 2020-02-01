CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = --std="c++17" -L/home/srossd/muser/src/api -I/home/srossd/muser/src/api

# the build target executable:
TARGET = main

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp -lmuser2_api

clean: 
	$(RM) $(TARGET)