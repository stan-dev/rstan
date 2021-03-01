# Set OS if there is no default:
uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')
OS = win
ifeq ($(uname_S), Darwin)
  OS = mac
endif
ifeq ($(uname_O), Cygwin)
  OS = win
endif
ifeq ($(uname_S), MINGW32_NT-6.1)
  OS = win
endif
ifeq ($(uname_S), Linux)
  OS = linux
endif

ifdef R_HOME
  R = $(R_HOME)/bin/R 
  RSCRIPT = $(R_HOME)/bin/Rscript 
else
  # R = R --arch x86_64
  R = R 
  RSCRIPT = Rscript 
endif

R_MAKEVARS_USER = $(CURDIR)/R_Makevars
export R_MAKEVARS_USER

BUILD_OPTS ?= 

RSTANVER := $(shell $(RSCRIPT) -e "cat(read.dcf(file = './rstan/DESCRIPTION')[1, deparse(quote(Version))])")
STANPKG := rstan_$(RSTANVER).tar.gz 

GIT = $(strip $(shell $(WH) git))

# The default target of this makefile 
help:

.PHONY: help
help: 
	@echo '--------------------------------------------------------------------------------'
	@echo 'main targets:'
	@echo ' build: to build rstan source package using R CMD build'
	@echo ' install: to install rstan using R CMD INSTALL' 
	@echo ' uninstall: to install rstan using R CMD REMOVE rstan'
	@echo ' clean: to remove the gzip package file'
	@echo ' check: use R CMD check on rstan package'
	@echo ' example_csv: create CSV files included in the package for examples'
	@echo ' install_pre_rpkg: to install R packages that rstan depends on'
	@echo ' test-R: run tests not packaged in rstan for package rstan'
	@echo ' test-cpp: to compile c++ code directly (for dev)'
	@echo '--------------------------------------------------------------------------------'

.PHONY: build check install clean clean-all install_pre_rpkg test-R example_csv

install_pre_rpkg:
	@R -q -e "options(repos=structure(c(CRAN = 'https://cran.rstudio.com'))); for (pkg in c('inline', 'Rcpp', 'RcppEigen', 'RUnit', 'BH', 'StanHeaders', 'RInside')) if (!require(pkg, character.only = TRUE))  install.packages(pkg, dep = TRUE); sessionInfo()"

build $(STANPKG): ./rstan/DESCRIPTION
	$(R) CMD build rstan --md5  $(BUILD_OPTS)   --no-build-vignettes --no-manual

example_csv: 
ifneq ($(OS), win)
	@if test -f "example/rstan_doc_ex_1.csv"; then cp example/rstan_doc_ex_*.csv ./rstan/inst/misc && echo "Copy done.";  \
	else echo "Not found example/rstan_doc_ex_1.csv, run example/examplemodel.R first."; fi
endif

check: $(STANPKG) 
	$(R) CMD check --as-cran $(STANPKG) 

install: $(STANPKG) 
	$(R) CMD INSTALL $(STANPKG)

uninstall:
	$(R) CMD REMOVE rstan

test-R:
	cd tests; R -q -f runRunitTests.R 

clean: 
	rm -f $(STANPKG) 
	rm -rf $(shell find tests/cpp -type f -name '*.d')
	rm -rf $(shell find tests/cpp -type f -name '*.d.*')

clean-all: clean
	rm -f $(shell find . -name 'rstan_*.tar.gz')
	rm -rf rstan.Rcheck
	rm -rf tests/cpp/lib $(GTEST_LIB)
	rm -rf $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) 

# buildbin:   # build a binary version  
# R CMD INSTALL -l ./tmp --build rstan


##
# C++ tests of rstan

STAN_SUBMODULE := ../StanHeaders/inst/include/upstream/
STAN_MATH_SUBMODULE := ../StanHeaders/inst/include/mathlib/

STANLIB := $(STAN_MATH_SUBMODULE)lib
STANSRC := $(STAN_SUBMODULE)src
GTEST_LIB := tests/cpp/gtest.o

LIBPTHREAD = 
ifeq ($(OS), linux)
  LIBPTHREAD += -lpthread
endif

CPPTESTS := $(shell find tests/cpp -type f -name '*_test.cpp')
STAN_INSTANTIATION_FILES := $(patsubst rstan/src/%.cpp,tests/cpp/lib/src/%.o,$(wildcard rstan/src/*.cpp))

## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R  --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R  --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell $(R_HOME)/bin/Rscript  -e 'RInside:::CxxFlags()')
RINSIDELIBS := 		$(shell $(R_HOME)/bin/Rscript  -e 'RInside:::LdFlags()')


$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): CPPFLAGS += -Wall $(RCPPFLAGS)

$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): CXXFLAGS += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS) -I $(STAN_MATH_SUBMODULE) -I rstan/inst/include -I $(STANSRC) $(addprefix -isystem ,$(wildcard $(STANLIB)/*)) -isystem $(GTESTPATH)/include -isystem $(GTESTPATH) -ftemplate-depth=256 -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS

$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): LDFLAGS += $(RLDFLAGS) $(RCPPLIBS) $(RINSIDELIBS) $(LIBPTHREAD)

GTESTPATH = $(wildcard $(STAN_MATH_SUBMODULE)lib/gtest*)

$(STAN_INSTANTIATION_FILES): tests/cpp/lib/src/%.o: rstan/src/%.cpp $(STANLIB) $(STANSRC)
	@mkdir -p $(dir $@)
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<  -I $(STANSRC)

$(GTEST_LIB): $(STANLIB) $(STANSRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(OUTPUT_OPTION) $(GTESTPATH)/src/gtest-all.cc

$(patsubst %.cpp,%.d,$(CPPTESTS)) : %.d : %.cpp $(STANLIB) $(STANSRC)
	@if test -d $(dir $@); \
	then \
	(set -e; \
	rm -f $@; \
	$(COMPILE.cpp) -MM $< > $@.$$$$; \
	sed -e 's,\($(notdir $*)\)\.o[ :]*,$(dir $@)\1$(EXE) $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$);\
	fi

$(patsubst %.cpp,%$(EXE),$(CPPTESTS)) : %$(EXE) : %.cpp %.d $(GTEST_LIB) $(STANLIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OUTPUT_OPTION) $(GTEST_LIB) $< $(GTESTPATH)/src/gtest_main.cc $(LDFLAGS)
ifeq (,$(filter test-cpp,$(MAKECMDGOALS)))
	$@
endif

.PHONY: test-cpp
test-cpp: $(patsubst %.cpp,%$(EXE),$(CPPTESTS))
	$(foreach test,$^,$(test); echo;)

ifneq (,$(filter test-cpp,$(MAKECMDGOALS)))
  -include $(patsubst %.cpp,%.d,$(CPPTESTS))
endif
ifneq (,$(filter $(patsubst %.cpp,%$(EXE),$(CPPTESTS)),$(MAKECMDGOALS)))
  -include $(patsubst %$(EXE),%.d,$(filter $(patsubst %.cpp,%$(EXE),$(CPPTESTS)),$(MAKECMDGOALS)))
endif

