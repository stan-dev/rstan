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
  R = $(R_HOME)/bin/R --vanilla
  RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
else
  # R = R --arch x86_64
  R = R --vanilla
  RSCRIPT = Rscript --vanilla
endif

R_MAKEVARS_USER = $(CURDIR)/R_Makevars
export R_MAKEVARS_USER

BUILD_OPTS ?= 
NOTBUILD ?= FALSE

ifeq ($(OS), win)
  LINKCMD1 = @cmd /c mklink /D ".^\rstan^\inst^\include^\stanlib"  "..^\..^\..^\..^\stan^\lib"
  LINKCMD2 = @cmd /c mklink /D ".^\rstan^\inst^\include^\stansrc"  "..^\..^\..^\..^\stan^\src"
  WH = where
else
  LINKCMD1 = ln -s ../../../../stan/lib ./rstan/inst/include/stanlib
  LINKCMD1b = ln -s ../../../$(TMP_STANLIB)/ ./rstan/inst/include/stanlib
  LINKCMD2 = ln -s ../../../../stan/src ./rstan/inst/include/stansrc
  WH = which
endif

RSTANVER := $(shell $(RSCRIPT) -e "cat(read.dcf(file = './rstan/DESCRIPTION')[1, deparse(quote(Version))])")
STANPKG := rstan_$(RSTANVER).tar.gz 

BOOSTPATH := $(shell find ../stan/lib -path '*lib/boost_*' -regex '.*lib\/boost_[^/]*')
EIGENPATH := $(shell find ../stan/lib -path '*lib/eigen_*' -regex '.*lib\/eigen_[^/]*')
TMP_STANLIB := tmpstanlib#
BOOST_FOR_RSTAN := $(TMP_STANLIB)/$(subst ../stan/lib/,,$(BOOSTPATH))
EIGEN_FOR_RSTAN := $(TMP_STANLIB)/$(subst ../stan/lib/,,$(EIGENPATH))

# a tool in boost 
BCP = $(strip $(shell $(WH) bcp))
GIT = $(strip $(shell $(WH) git))

# The default target of this makefile 
help:

.PHONY: help
help: 
	@echo '--------------------------------------------------------------------------------'
	@echo 'main targets:'
	@echo ' build: to build rstan source package in a gzip form using R CMD build'
	@echo ' install: to install rstan using R CMD INSTALL' 
	@echo ' uninstall: to install rstan using R CMD REMOVE rstan'
	@echo ' clean: to remove the gzip package file'
	@echo ' check: use R CMD check on rstan package'
	@echo ' tmpstanlib: use Boost bcp to strip boost'
	@echo ' example_csv: create CSV files included in the package for examples'
	@echo ' install_pre_rpkg: to install R packages of Rcpp, inline, RcppEigen, BH, and RUnit'
	@echo ' test-cpp: to compile c++ code directly (for dev)'
	@echo '--------------------------------------------------------------------------------'

install_pre_rpkg:
	@R -q -e "options(repos=structure(c(CRAN = 'http://cran.rstudio.com'))); for (pkg in c('inline', 'Rcpp', 'RcppEigen', 'RUnit', 'BH')) if (!require(pkg, character.only = TRUE))  install.packages(pkg, dep = TRUE); sessionInfo()" 


STANLIB := rstan/inst/include/stanlib
STANSRC := rstan/inst/include/stansrc

$(STANLIB):
ifeq ($(OS), win)
	$(LINKCMD1)
else 
	if test -d $(TMP_STANLIB); then $(LINKCMD1b); else $(LINKCMD1); fi
endif

$(STANSRC):
	$(LINKCMD2)


build $(STANPKG): ./rstan/DESCRIPTION $(STANSRC)
	@test -f "../stan/src/stan/version.hpp" || { echo "Stan version header file is not found; maybe the submodule of Stan is not initialized."; exit 1;}
ifeq ($(GIT),)
	@echo "git_head <- function() \"unknown\"" > rstan/R/git_rstan_head.R
	@echo "git_stan_head <- function() \"unknown\"" > rstan/R/git_stan_head.R
else
	@echo "git_head <- function() \"`git rev-parse HEAD`\"" > rstan/R/git_rstan_head.R
	@echo "git_stan_head <- function() \"`cd ..; git submodule status`\"" > rstan/R/git_stan_head.R
endif
ifeq ($(NOTBUILD),FALSE)
	$(R) CMD build rstan --md5  $(BUILD_OPTS)  # --no-build-vignettes --no-manual
endif


example_csv : 
ifneq ($(OS), win)
	@if test -f "example/rstan_doc_ex_1.csv"; then cp example/rstan_doc_ex_*.csv ./rstan/inst/misc && echo "Copy done.";  \
	else echo "Not found example/rstan_doc_ex_1.csv, run example/examplemodel.R first."; fi
endif

check: $(STANPKG) 
	$(R) CMD check --as-cran $(STANPKG) 

install: $(STANPKG) 
ifeq ($(OS), win)
	$(R) CMD INSTALL --merge-multiarch $(STANPKG) 
else
#       $(R) CMD INSTALL --no-multiarch $(STANPKG) 
	$(R) CMD INSTALL $(STANPKG)
endif

uninstall:
	$(R) CMD REMOVE rstan

$(TMP_STANLIB): 
ifneq "$(BCP)" ""
	mkdir -p ${BOOST_FOR_RSTAN} 
	mkdir -p ${EIGEN_FOR_RSTAN} 
	cp -r $(EIGENPATH)/* $(EIGEN_FOR_RSTAN)
	find ../stan/src/ -name \*\.\[ch]pp -exec bcp --scan --boost=${BOOSTPATH} '{}' ${BOOST_FOR_RSTAN} \; &> bcp.log  
	find ./rstan/inst/include -name \*\.\[ch]pp -exec bcp --scan --boost=${BOOSTPATH} '{}' ${BOOST_FOR_RSTAN} \; >> bcp.log  2>&1
else 
	@echo "cannot find boost bcp"
endif 


rmtmpstanlib: 
	rm -rf $(TMP_STANLIB) 
	rm -f bcp.log 

clean: 
	rm -f $(STANPKG) 
	rm -f rstan/R/git_*head.R
	rm -rf $(shell find tests/cpp -type f -name '*.d')
	rm -rf $(shell find tests/cpp -type f -name '*.d.*')

clean-all: clean
	rm -f $(STANPKG); rm -rf $(TMP_STANLIB); rm -f bcp.log 
	rm -rf rstan/R/git_*stan_head.R 
	rm -rf "./$(STANLIB)"
	rm -rf "./$(STANSRC)"
	rm -rf rstan.Rcheck
	rm -rf tests/cpp/lib $(GTEST_LIB)
	rm -rf $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) 

# buildbin:   # build a binary version  
# R CMD INSTALL -l ./tmp --build rstan

.PHONY: build check install clean clean-all $(TMP_STANLIB) rmtmpstanlib install_pre_rpkg

##
# C++ tests of rstan
GTEST_LIB := tests/cpp/gtest.o

CPPTESTS := $(shell find tests/cpp -type f -name '*_test.cpp')
STAN_LIB := tests/cpp/lib/libstan.a
STAN_INSTANTIATION_FILES := $(patsubst rstan/src/%.cpp,tests/cpp/lib/src/%.o,$(wildcard rstan/src/*.cpp))

## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): CPPFLAGS += -Wall $(RCPPFLAGS)

$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): CXXFLAGS += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS) -I rstan/inst/include -I $(STANSRC) $(addprefix -isystem ,$(wildcard $(STANLIB)/*)) -isystem $(GTESTPATH)/include -ftemplate-depth=256 -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS

$(STAN_INSTANTIATION_FILES) $(GTEST_LIB) $(patsubst %.cpp,%$(EXE),$(CPPTESTS)) $(patsubst %.cpp,%.d,$(CPPTESTS)): LDFLAGS += $(RLDFLAGS) $(RCPPLIBS) $(RINSIDELIBS) -L $(dir $(STAN_LIB)) -l stan

GTESTPATH = $(wildcard $(STANLIB)/gtest*)

$(STAN_INSTANTIATION_FILES): tests/cpp/lib/src/%.o: rstan/src/%.cpp $(STANLIB) $(STANSRC)
	@mkdir -p $(dir $@)
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<  -I $(STANSRC)

$(STAN_LIB): $(STAN_INSTANTIATION_FILES)
	@mkdir -p $(dir $@)
	$(AR) -rs $@ $(STAN_INSTANTIATION_FILES)

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

$(patsubst %.cpp,%$(EXE),$(CPPTESTS)) : %$(EXE) : %.cpp %.d $(GTEST_LIB) $(STANLIB) $(STAN_LIB)
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


.PHONY: test-R
test-R:
	cd tests; R -q -f runRunitTests.R 
