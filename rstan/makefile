
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
	@echo '--------------------------------------------------------------------------------'

build $(STANPKG): ./rstan/DESCRIPTION 
	@rm -rf "./rstan/inst/include/stanlib"
	@rm -rf "./rstan/inst/include/stansrc"
ifeq ($(OS), win)
	$(LINKCMD1)
else 
	if test -d $(TMP_STANLIB); then $(LINKCMD1b); else $(LINKCMD1); fi
endif
	$(LINKCMD2)
	@rm -f $(STANPKG) 
ifeq ($(GIT),)
	@echo "git_head <- function() \"unknown\"" > rstan/R/git_rstan_head.R
	@echo "git_stan_head <- function() \"unknown\"" > rstan/R/git_stan_head.R
else
	@echo "git_head <- function() \"`git rev-parse HEAD`\"" > rstan/R/git_rstan_head.R
	@echo "git_stan_head <- function() \"`cd ..; git submodule status`\"" > rstan/R/git_stan_head.R
endif
	$(R) CMD build rstan --md5  # --no-vignettes --no-manual

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
	# $(R) CMD INSTALL --no-multiarch $(STANPKG) 
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

clean-all: 
	rm -f $(STANPKG); rm -rf $(TMP_STANLIB); rm -f bcp.log 
	rm -rf "./rstan/inst/include/stanlib"
	rm -rf "./rstan/inst/include/stansrc"
	rm -rf rstan.Rcheck
  

# buildbin:   # build a binary version  
# R CMD INSTALL -l ./tmp --build rstan

.PHONY: build check install clean clean-all $(TMP_STANLIB) rmtmpstanlib
