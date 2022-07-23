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

GA_LIBSFILES := $(GA_SUBDIR)/$(SOURCE)/* $(GA_SUBDIR)/$(INCLUDE)/*

main: $(SOURCE)/main.cpp $(GA_LIBSFILES) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOURCE)/main.cpp -o $(BUILD)/$(TARGET)

$(GA_LIBSFILES): GraphAligner
	@# Do nothing

.PHONY: GraphAligner
GraphAligner:
	cd $(GA_SUBDIR)
	conda env create -f CondaEnvironment_osx.yml
	cd ..
	source activate GraphAligner
	$(MAKE) -C $(GA_SUBDIR)/bin/GraphAligner
	
$(BUILD):
	-mkdir -p $@
	
clean:
	$(MAKE) -C $(GA_SUBDIR) clean
	$(RM) -r build
