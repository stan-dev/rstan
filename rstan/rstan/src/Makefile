all: symbols.rds $(SHLIB) libstanLibrary 
LIBSTAN = ../inst/libstan
SHLIB = rstan.so

USERLIB = libstan.a

OBJECTS = stanc.o chains.o misc.o init.o
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard gm__*.cpp))

LIBSTANOBJECTS = agrad__rev__var_stack.o

PKG_LIBS = $(shell  $(R_HOME)/bin/Rscript --vanilla -e "Rcpp:::LdFlags()") 
PKG_CPPFLAGS += -I"../inst/include/stansrc"  

PKG_CPPFLAGS += -I"../inst/include" 
PKG_CPPFLAGS += -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS
PKG_CPPFLAGS += -ftemplate-depth-256

SHLIB_LDFLAGS = $(SHLIB_CXXLDFLAGS) 
SHLIB_LD = $(SHLIB_CXXLD)

$(USERLIB): $(LIBSTANOBJECTS)
	$(AR) qc $@ $^ 
	@if test -n "$(RANLIB)"; then $(RANLIB) $@; fi


libstanLibrary: $(USERLIB) 
	# -@if test ! -e $(LIBSTAN)$(R_ARCH); then mkdir -p $(LIBSTAN)$(R_ARCH); fi
	# cp $(USERLIB) $(LIBSTAN)$(R_ARCH)
	# rm $(USERLIB) 

# $(USERLIB): $(LIBSTANOBJECTS)
#   $(SHLIB_CXXLD) -s -o $(USERLIB) $(LIBSTANOBJECTS) $(SHLIB_CXXLDFLAGS) $(ALL_LIBS)
#   @if test -e "/usr/bin/install_name_tool"; then /usr/bin/install_name_tool -id $(R_PACKAGE_DIR)/libstan/$(R_ARCH)/$(USERLIB) $(USERLIB); fi

include $(R_SHARE_DIR)/make/shlib.mk

.PHONY: all clean libstanLibrary

clean: 
	rm -f $(OBJECTS) $(SHLIB) $(USERLIB) $(LIBSTANOBJECTS) 

# P.S. While Makefile is used, only the binary of one platform (32bits or
# 64bits) is installed, for example on Mac and Windows. 
