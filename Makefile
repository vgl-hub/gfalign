CXX = g++
INCLUDE_DIR = -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g --std=c++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalign
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

GA_SUBDIR := $(CURDIR)/GraphAligner
GFALIBS_SUBDIR := $(CURDIR)/gfalibs

LIBS = -lz
LDFLAGS= -pthread

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

#GraphAligner
GA_LIBSFILES := $(GA_SUBDIR)/$(SOURCE)/* $(GA_SUBDIR)/$(INCLUDE)/*

#gfalign
GFALIGN_OBJS := main alignments input eval
GFALIGN_BINS := $(addprefix $(BINDIR)/, $(GFALIGN_OBJS))

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

head: $(GFALIGN_BINS) $(GA_LIBSFILES) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

$(GFALIGN_OBJS): %: $(BINDIR)/%
	@

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h $(CURDIR)/gfalibs/include/*.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

$(GA_LIBSFILES): GraphAligner
	@# Do nothing

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"

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
	mv GraphAligner/bin/* build/bin
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
	
clean:
	$(MAKE) -C $(GA_SUBDIR) clean
	$(RM) -r build
	$(RM) -r $(MY_ENV_DIR)
