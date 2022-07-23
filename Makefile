CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalign
BUILD = build/bin
SOURCE = src
INCLUDE = include

GA_SUBDIR := $(CURDIR)/GraphAligner
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

main: $(SOURCE)/main.cpp $(GA_LIBSFILES) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOURCE)/main.cpp -o $(BUILD)/$(TARGET)

$(GA_LIBSFILES): GraphAligner
	@# Do nothing

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
	
clean:
	$(MAKE) -C $(GA_SUBDIR) clean
	$(RM) -r build
	$(RM) -r $(MY_ENV_DIR)