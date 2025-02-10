CXX = g++
INCLUDE_DIR = -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g --std=c++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalign
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

GFALIBS_SUBDIR := $(CURDIR)/gfalibs

LIBS = -lz
LDFLAGS = -pthread

ifeq (,$(shell which conda))
    HAS_CONDA = False
else
    HAS_CONDA = True
    CONDA_ENV_NAME = GraphAligner
    ENV_DIR = $(shell conda info --base)
    MY_ENV_DIR = $(ENV_DIR)/envs/$(CONDA_ENV_NAME)
    CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
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

#gfalign
GFALIGN_OBJS := main alignments input-gfalign eval
GFALIGN_BINS := $(addprefix $(BINDIR)/, $(GFALIGN_OBJS))

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

head: $(GFALIGN_BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

all: head validate regenerate

$(GFALIGN_OBJS): %: $(BINDIR)/%
	@

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h $(CURDIR)/gfalibs/include/*.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

$(GA_LIBSFILES): GraphAligner
	@# Do nothing

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"
	
validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp

regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp

.PHONY: GraphAligner
.ONESHELL:

SHELL = /bin/bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

GraphAligner:
ifeq (True,$(HAS_CONDA))
ifneq ("$(wildcard $(MY_ENV_DIR))","")
	@echo ">>> Found GraphAligner environment in $(MY_ENV_DIR). Skipping installation..."
else
	@echo ">>> Detected conda, but $(CONDA_ENV_NAME) is missing in $(ENV_DIR). Installing ..."
	conda install -c bioconda graphaligner
endif
else
	@echo ">>> Install conda first."
	exit
endif
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
	
clean:
	$(RM) -r build
	$(RM) -r $(MY_ENV_DIR)
	$(MAKE) -C $(GFALIBS_DIR) clean
