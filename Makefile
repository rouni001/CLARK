SHELL := /bin/sh

CXX ?= g++
CXXFLAGS ?= -O3
CPPFLAGS ?=
LDFLAGS ?=

BUILD_DIR ?= build
EXE_DIR ?= exe
OPENMP ?= auto

OPENMP_TEST := $(shell tmp="$${TMPDIR:-/tmp}/clark-openmp-test-$$$$"; printf 'int main(){return 0;}\n' | $(CXX) -x c++ -fopenmp - -o "$$tmp" >/dev/null 2>&1 && { rm -f "$$tmp"; echo yes; } || { rm -f "$$tmp"; echo no; })

ifeq ($(OPENMP),1)
  OPENMP_FLAGS := -fopenmp
else ifeq ($(OPENMP),true)
  OPENMP_FLAGS := -fopenmp
else ifeq ($(OPENMP),auto)
  ifeq ($(OPENMP_TEST),yes)
    OPENMP_FLAGS := -fopenmp
  else
    OPENMP_FLAGS :=
  endif
else
  OPENMP_FLAGS :=
endif

HELPERS := \
	$(EXE_DIR)/getTargetsDef \
	$(EXE_DIR)/getAccssnTaxID \
	$(EXE_DIR)/getfilesToTaxNodes \
	$(EXE_DIR)/getAbundance \
	$(EXE_DIR)/getConfidenceDensity \
	$(EXE_DIR)/getGammaDensity \
	$(EXE_DIR)/makeSummaryTables \
	$(EXE_DIR)/converter \
	$(EXE_DIR)/exeSeq \
	$(EXE_DIR)/dscriptMaker \
	$(EXE_DIR)/getTargetSpecificKmersStat \
	$(EXE_DIR)/extractSeqs

VARIANTS := $(EXE_DIR)/CLARK $(EXE_DIR)/CLARK-l $(EXE_DIR)/CLARK-S

.PHONY: all helpers variants clean test openmp-status

all: openmp-status helpers variants

helpers: $(HELPERS)

variants: $(VARIANTS)

openmp-status:
	@if [ "$(OPENMP_FLAGS)" = "-fopenmp" ]; then \
		echo "OpenMP: enabled"; \
	else \
		echo "OpenMP: disabled (compiler did not accept -fopenmp; build will be single-threaded)"; \
	fi

$(EXE_DIR):
	mkdir -p "$(EXE_DIR)"

$(BUILD_DIR):
	mkdir -p "$(BUILD_DIR)"

$(EXE_DIR)/getTargetsDef: src/getTargetsDef.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getAccssnTaxID: src/getAccssnTaxID.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getfilesToTaxNodes: src/getfilesToTaxNodes.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getAbundance: src/getAbundance.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getConfidenceDensity: src/getConfidencedensity.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getGammaDensity: src/getGammadensity.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/makeSummaryTables: src/file.cc src/makeSamplesSummaryTables.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/converter: src/main_spaced.cc src/kmersConversion.cc src/contiguousToSpaced_hh.hh src/hashTable_hh.hh | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" src/main_spaced.cc src/kmersConversion.cc $(LDFLAGS)

$(EXE_DIR)/exeSeq: src/getSeqFiles.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/dscriptMaker: src/dscriptMaker.cc src/file.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/getTargetSpecificKmersStat: src/file.cc src/getTargetSpecificKmersStat.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(EXE_DIR)/extractSeqs: src/file.cc src/extractSequences.cc | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o "$@" $^ $(LDFLAGS)

$(BUILD_DIR)/default/.prepared: src/*.cc src/*.hh src/parameters.hh | $(BUILD_DIR)
	rm -rf "$(BUILD_DIR)/default"
	mkdir -p "$(BUILD_DIR)/default"
	cp src/*.hh "$(BUILD_DIR)/default/"
	cp src/main.cc src/analyser.cc src/file.cc src/kmersConversion.cc src/FileHandler*.cc "$(BUILD_DIR)/default/"
	cp src/parameters.hh "$(BUILD_DIR)/default/parameters.hh"
	touch "$@"

$(BUILD_DIR)/light/.prepared: src/*.cc src/*.hh src/parameters_hh | $(BUILD_DIR)
	rm -rf "$(BUILD_DIR)/light"
	mkdir -p "$(BUILD_DIR)/light"
	cp src/*.hh "$(BUILD_DIR)/light/"
	cp src/main.cc src/analyser.cc src/file.cc src/kmersConversion.cc src/FileHandler*.cc "$(BUILD_DIR)/light/"
	cp src/parameters_hh "$(BUILD_DIR)/light/parameters.hh"
	touch "$@"

$(BUILD_DIR)/spaced/.prepared: src/*.cc src/*.hh src/parameters_shh | $(BUILD_DIR)
	rm -rf "$(BUILD_DIR)/spaced"
	mkdir -p "$(BUILD_DIR)/spaced"
	cp src/*.hh "$(BUILD_DIR)/spaced/"
	cp src/main.cc src/analyser.cc src/file.cc src/kmersConversion.cc src/FileHandler*.cc "$(BUILD_DIR)/spaced/"
	cp src/parameters_shh "$(BUILD_DIR)/spaced/parameters.hh"
	touch "$@"

$(EXE_DIR)/CLARK: $(BUILD_DIR)/default/.prepared | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -o "$@" $(BUILD_DIR)/default/*.cc $(LDFLAGS) $(OPENMP_FLAGS)

$(EXE_DIR)/CLARK-l: $(BUILD_DIR)/light/.prepared | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -o "$@" $(BUILD_DIR)/light/*.cc $(LDFLAGS) $(OPENMP_FLAGS)

$(EXE_DIR)/CLARK-S: $(BUILD_DIR)/spaced/.prepared | $(EXE_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) -o "$@" $(BUILD_DIR)/spaced/*.cc $(LDFLAGS) $(OPENMP_FLAGS)

test: all
	tests/run_tests.sh

clean:
	rm -rf "$(BUILD_DIR)" "$(EXE_DIR)"
