# Directories
SRC_DIR = ./
BIN_DIR = bin
TEST_DIR = tests
DEST_DIR = ~/bin

# Compiler settings
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
CC = gcc
LDFLAGS = -lm -lz -lpthread

# Source path
VPATH = $(SRC_DIR)

# Output binaries
ALL = $(BIN_DIR)/FasTAN $(BIN_DIR)/aln2bed

# Source files
FASTAN_SRCS = $(SRC_DIR)/FasTAN.c $(SRC_DIR)/alncode.c $(SRC_DIR)/align.c \
              $(SRC_DIR)/GDB.c $(SRC_DIR)/gene_core.c $(SRC_DIR)/ONElib.c
FASTAN_HDRS = $(SRC_DIR)/alncode.h $(SRC_DIR)/GDB.h $(SRC_DIR)/ONElib.h \
              $(SRC_DIR)/align.h $(SRC_DIR)/gene_core.h

ALN2BED_SRCS = $(SRC_DIR)/aln2bed.c $(SRC_DIR)/alncode.c $(SRC_DIR)/GDB.c \
               $(SRC_DIR)/gene_core.c $(SRC_DIR)/ONElib.c
ALN2BED_HDRS = $(SRC_DIR)/alncode.h $(SRC_DIR)/GDB.h $(SRC_DIR)/ONElib.h \
               $(SRC_DIR)/align.h $(SRC_DIR)/gene_core.h

# Default target
all: $(BIN_DIR) $(ALL)

# Create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Dependencies
$(SRC_DIR)/GDB.c: $(SRC_DIR)/gene_core.c $(SRC_DIR)/gene_core.h
$(SRC_DIR)/GDB.h: $(SRC_DIR)/gene_core.h

# FasTAN binary
$(BIN_DIR)/FasTAN: $(FASTAN_SRCS) $(FASTAN_HDRS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $(FASTAN_SRCS) $(LDFLAGS)

# aln2bed binary
$(BIN_DIR)/aln2bed: $(ALN2BED_SRCS) $(ALN2BED_HDRS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $(ALN2BED_SRCS) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(BIN_DIR)/FasTAN $(BIN_DIR)/aln2bed
	rm -fr $(BIN_DIR)/*.dSYM
	rm -f FasTAN.tar.gz

# Install binaries to destination
install: all
	cp $(BIN_DIR)/FasTAN $(BIN_DIR)/aln2bed $(DEST_DIR)

# Create source package
package:
	make clean
	tar -zcf FasTAN.tar.gz LICENSE README.md Makefile \
		$(SRC_DIR)/*.h $(SRC_DIR)/*.c

# Run tests
test: all
	@echo "Running tests..."
	@if [ -f $(TEST_DIR)/run_tests.sh ]; then \
		bash $(TEST_DIR)/run_tests.sh; \
	else \
		echo "No tests found in $(TEST_DIR)/"; \
	fi

# Help target
help:
	@echo "FasTAN Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build all binaries (default)"
	@echo "  clean     - Remove built binaries"
	@echo "  install   - Install binaries to $(DEST_DIR)"
	@echo "  test      - Run tests"
	@echo "  package   - Create source package"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Directories:"
	@echo "  Source:   $(SRC_DIR)/"
	@echo "  Binaries: $(BIN_DIR)/"
	@echo "  Tests:    $(TEST_DIR)/"

.PHONY: all clean install package test help
