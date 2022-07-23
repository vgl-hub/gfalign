CXX = g++
INCLUDE_DIR = -I./include -I../gfastats/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalign
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

GA_SUBDIR := $(CURDIR)/GraphAligner
GFASTATS_SUBDIR := $(CURDIR)/../gfastats

LIBS = -lz
LDFLAGS :=

ifeq (,$(shell which conda))
    HAS_CONDA=False
else
    HAS_CONDA=True
    ENV_DIR=$(shell conda info --base)
    MY_ENV_DIR=$(ENV_DIR)/envs/GraphAligner
    CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
endif

OSF :=
ifeq ($(OS),Windows_NT)
	OSF += WIN32
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSF = LINUX
	else ifeq ($(UNAME_S),Darwin)
		OSF = OSX
	endif
endif

GA_LIBSFILES := $(GA_SUBDIR)/$(SOURCE)/* $(GA_SUBDIR)/$(INCLUDE)/*
GFALIGN_LIBSFILES := main

OBJS := stream-obj bed struct log functions
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(GFASTATS_SUBDIR)/$(INCLUDE)/threadpool.h $(BINS) $(INCLUDE)/$(GFALIGN_LIBSFILES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(LIBS)

$(INCLUDE)/%: $(INCLUDE)/%.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $(BINDIR)/$(notdir $@)

$(GA_LIBSFILES): GraphAligner
	@# Do nothing
	
$(BINDIR)%: $(GFASTATS_SUBDIR)/$(SOURCE)/%.cpp $(GFASTATS_SUBDIR)/$(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(GFASTATS_SUBDIR)/$(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: GraphAligner
GraphAligner:
ifeq (True,$(HAS_CONDA))
ifneq ("$(wildcard $(MY_ENV_DIR))","")
	@echo ">>> Found GraphAligner environment in $(MY_ENV_DIR). Skipping installation..."
else
	@echo ">>> Detected conda, but $(CONDA_ENV_NAME) is missing in $(ENV_DIR). Installing ..."
ifeq ($(OSF),LINUX)
	conda env create -f $(GA_SUBDIR)/CondaEnvironment_linux.yml
else ifeq ($(OSF),OSX)
	conda env create -f $(GA_SUBDIR)/CondaEnvironment_osx.yml
else
	@echo ">>> OS $(OSF) not supported by GraphAligner ..."
	exit
endif
endif
else
	@echo ">>> Install conda first."
    exit
endif
	. activate GraphAligner && $(MAKE) -j -C $(GA_SUBDIR)
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
	
clean:
	$(MAKE) -C $(GA_SUBDIR) clean
	$(RM) -r build
	$(RM) -r $(MY_ENV_DIR)
