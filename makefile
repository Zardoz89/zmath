# Makefile to build zmath
#
# make debug => makes debug build of the library
#
# make release => makes release build of the library
#
# make unittest => builds all unittests (for both debug and release)
# and runs them
#
# make doc => makes documentation

DDOC = dmd
DMD = dmd
DFLAGS =

SRC_DIR = ./src
DOC_DIR = ./doc
BIN_DIR = ./build
LIB_DIR = ./lib
IMPORT_DIR = ./import

MODULE = zmath
LIB = lib$(MODULE).a
TARGET = $(LIB_DIR)/$(LIB)
TEST_TARGET = $(BIN_DIR)/test

SRC_FILES = $(SRC_DIR)/aux.d $(SRC_DIR)/vector.d $(SRC_DIR)/matrix.d $(SRC_DIR)/quaternion.d $(SRC_DIR)/math3d.d
SRC_TEST_FILES = $(SRC_DIR)/unittest.d

# Build can be debug or release
BUILD =
# 32 or 64 bits
MODEL =

# Gets this makefile name
MAKEFILE:=$(lastword $(MAKEFILE_LIST))

# DDoc flags
DDOCFLAGS =
DDOCFLAGS :=-c -o-

# Set DFLAGS
DFLAGS := -w
ifeq ($(BUILD),debug)
	DFLAGS += -gc -debug
else
	DFLAGS += -O -release
endif

# Set flags if MODEL it's set
ifeq ($(MODEL),)
# Default 32 or 64 bits target of system
else
	DFLGS += -m$(MODEL)
	DDOCFLAGS += -m$(MODEL)
endif

ifeq ($(BUILD),)
# If BUILD it's empty, recall-self setting BUILD to apropiated value

release :
	$(MAKE) --no-print-directory -f$(MAKEFILE) MODEL=$(MODEL) BUILD=release
debug :
	$(MAKE) --no-print-directory -f$(MAKEFILE) MODEL=$(MODEL) BUILD=debug
unittest :
	$(MAKE) --no-print-directory -f$(MAKEFILE) MODEL=$(MODEL) BUILD=debug unittest
#	$(MAKE) --no-print-directory -f$(MAKEFILE) MODEL=$(MODEL) BUILD=release unittest

else

########################################################
# Real compiling stuff

$(BUILD) : $(TARGET)

$(TARGET) : $(SRC_FILES)
	$(DMD) $(DFLAGS) -lib -of$(TARGET) $(SRC_FILES) -H -Hd$(IMPORT_DIR)/$(MODULE)
	$(DMD) $(DFLAGS) $(SRC_FILES) -o- -H -Hd$(IMPORT_DIR)/$(MODULE)


# Make unittest and run it
unittest : $(TARGET) $(SRC_TEST_FILES) $(SRC_FILES)
	$(DMD) $(DFLAGS) -unittest -of$(TEST_TARGET) $(SRC_TEST_FILES) -I$(IMPORT_DIR) -L-L$(LIB_DIR)/ -L-l$(MODULE)
	$(TEST_TARGET)
endif

########################################################
# Documentation generator

DOCS = $(DOC_DIR)/aux.html \
	$(DOC_DIR)/vector.html	\
	$(DOC_DIR)/matrix.html	\
	$(DOC_DIR)/quaternion.html	\
	$(DOC_DIR)/math3d.html	\

doc : $(DOCS)

$(DOC_DIR)/aux.html : aux.d
	$(DDOC) $(DDOCFLAGS) -Df$@ $< -I$(SRC_DIR)

$(DOC_DIR)/vector.html : vector.d
	$(DDOC) $(DDOCFLAGS) -Df$@ $< -I$(SRC_DIR)

$(DOC_DIR)/matrix.html : vector.d matrix.d
	$(DDOC) $(DDOCFLAGS) -Df$@ $^ -I$(SRC_DIR)

$(DOC_DIR)/quaternion.html : vector.d matrix.d quaternion.d
	$(DDOC) $(DDOCFLAGS) -Df$@ $^ -I$(SRC_DIR)

$(DOC_DIR)/math3d.html : matrix.d quaternion.d math3d.d
	$(DDOC) $(DDOCFLAGS) -Df$@ $^ -I$(SRC_DIR)

########################################################
#keep clean your home!
.PHONY: clean
clean :
	rm -rf $(TARGET) $(TEST_TARGET) $(DOC_DIR)/* $(IMPORT_DIR)/$(MODULE)/*

