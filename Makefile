# CLOVERS Makefile
# Ab initio prediction of overprinted genes using the Z-curve method

# Compiler settings
CXX := g++

# Allow environment to override flags
CXXFLAGS ?= -DZLIB -fopenmp -mavx -mfma -static -O3
LDFLAGS ?= -lz

# Append extra flags if provided
CXXFLAGS += $(EXTRA_CXXFLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)

# Directories
SRC_DIR := src
INC_DIR := include
BIN_DIR := bin
BUILD_DIR := build

# Source files
SOURCES := $(SRC_DIR)/Main.cpp \
           $(SRC_DIR)/BioIO.cpp \
           $(SRC_DIR)/BioUtil.cpp \
           $(SRC_DIR)/Encoding.cpp \
           $(SRC_DIR)/Model.cpp \
           $(SRC_DIR)/svm.cpp

# Object files
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# Binary data file
META_BIN := $(BIN_DIR)/meta.bin
META_OBJ := $(BUILD_DIR)/meta.bin.o

# Target executable
TARGET := clovers

# Default target
.PHONY: all clean

all: $(BUILD_DIR) $(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir $(BUILD_DIR)

# Convert binary data to object file
# Note: objcopy generates symbols based on the input file path.
# For bin/meta.bin, it generates _binary_bin_meta_bin_start/end
# We need to rename them to match the expected _binary_meta_bin_start/end
$(META_OBJ): $(META_BIN)
	objcopy --input binary --output elf64-x86-64 --binary-architecture i386:x86-64 $< $@
	objcopy --redefine-sym _binary_bin_meta_bin_start=_binary_meta_bin_start $@
	objcopy --redefine-sym _binary_bin_meta_bin_end=_binary_meta_bin_end $@

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Link executable
$(TARGET): $(META_OBJ) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Clean build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Rebuild
rebuild: clean all
