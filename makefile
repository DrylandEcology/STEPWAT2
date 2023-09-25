CC  	=	gcc
CXX   = g++

CFLAGS = \
	-DSQLITE_OMIT_LOAD_EXTENSION \
	-DSQLITE_THREADSAFE=0 \
	-DSQLITE_WITHOUT_ZONEMALLOC \
	-DSTEPWAT \
	-g \
	-O0 \
	-pthread \
	-Wcast-align \
	-Wformat \
	-Wimplicit \
	-Wmissing-prototypes \
	-Wredundant-decls \
	-Wstrict-prototypes \
	-Wunused

sw2 = SOILWAT2
lib_sw2 = lib$(sw2).a
path_sw2 = sw_src
path_sw2lib = $(path_sw2)/bin

std_c++14 = -std=c++14
std_c99 = -std=c99

# Note: `-I$(path_sw2)` is required for `#include "external/pcg/pcg_basic.h"`
# in SOILWAT2 headers that are included by STEPWAT2 code
INC_DIRS = \
	-I. \
	-Isqlite-amalgamation \
	-I$(path_sw2) \
	-I$(path_sw2)/external/googletest/googletest \
	-I$(path_sw2)/external/googletest/googletest/include \
	-Itest

sources_core = \
	sqlite-amalgamation/sqlite3.c \
	ST_environs.c \
	ST_grid.c \
	ST_indivs.c \
	ST_main.c \
	ST_mortality.c \
	ST_output.c \
	ST_params.c \
	ST_resgroups.c \
	ST_species.c \
	ST_sql.c \
	ST_stats.c \
	sxw.c \
	sxw_environs.c \
	sxw_resource.c \
	sxw_soilwat.c \
	sxw_sql.c \
	ST_spinup.c \
	ST_progressBar.c \
	ST_seedDispersal.c \
	ST_colonization.c

sources_test = \
	$(path_sw2)/external/googletest/googletest/src/gtest-all.cc \
	$(path_sw2)/external/googletest/googletest/src/gtest_main.cc \
	test/test_ST_mortality.cc

sw2_sources = \
	src/SW_Output_outarray.c \
	src/SW_Output_outtext.c


objects_core = $(sources_core:%.c=obj/%.o)
objects_core_test = $(sources_core:%.c=obj/%_TEST.o)
objects_test = $(sources_test:%.cc=obj/%.o)

dir_obj = \
	obj \
	obj/sqlite-amalgamation
dir_test = \
	obj/test \
	obj/$(path_sw2)/external/googletest/googletest/src

sw_LDFLAGS = $(LDFLAGS) -L. -L$(path_sw2lib)
sw_LDLIBS = -l$(sw2) $(LDLIBS) -lm



all: $(path_sw2lib)/$(lib_sw2) stepwat

$(path_sw2lib)/$(lib_sw2):
# Note: `-I..` is required for `#include "ST_defines.h"`
# in SOILWAT2 headers that are included by SOILWAT2 code
	@(cd $(path_sw2) && $(MAKE) lib \
		CC="$(CC)" CPPFLAGS="$(CPPFLAGS) -I.." CFLAGS="$(CFLAGS)" set_std="$(std_c99)" AR="$(AR)" \
		sw_sources="$(sw2_sources)")

stepwat: $(path_sw2lib)/$(lib_sw2) $(objects_core) | $(dir_obj)
	$(CC) $(objects_core) $(CFLAGS) $(std_c99) $(CPPFLAGS) $(sw_LDLIBS) $(sw_LDFLAGS) -o stepwat
	-@cp stepwat testing.sagebrush.master
	-@cp stepwat testing.sagebrush.master/Stepwat_Inputs

stepwat_test: $(path_sw2lib)/$(lib_sw2) $(objects_core_test) $(objects_test) | $(dir_obj) $(dir_test)
	$(CXX) $(objects_core_test) $(objects_test) $(std_c++14) $(CFLAGS) $(CPPFLAGS) $(sw_LDLIBS) $(sw_LDFLAGS) -o stepwat_test

obj/%.o: %.c | $(dir_obj)
	$(CC) $(CFLAGS) $(std_c99) $(CPPFLAGS) $(INC_DIRS) -c $< -o $@

obj/%.o: %.cc | $(dir_obj)
	$(CXX) $(CFLAGS) $(std_c++14) $(CPPFLAGS) $(INC_DIRS) -c $< -o $@

obj/%_TEST.o: %.c | $(dir_obj) $(dir_test)
	$(CC) $(CFLAGS) $(std_c99) $(CPPFLAGS) $(INC_DIRS) -DSTDEBUG -c $< -o $@


#--- Create directories
$(dir_obj) $(dir_test):
		-@mkdir -p $@


.PHONY: run_tests
run_tests: stepwat_test
	./stepwat_test

.PHONY: bint_testing_nongridded
bint_testing_nongridded: stepwat
	testing.sagebrush.master/Stepwat_Inputs/stepwat -d testing.sagebrush.master/Stepwat_Inputs -f files.in -o -i

.PHONY: bint_testing_gridded
bint_testing_gridded: stepwat
	testing.sagebrush.master/stepwat -d testing.sagebrush.master -f files.in -g

.PHONY: cleanall
cleanall: clean output_clean

.PHONY: clean
clean: cleanobjs cleanbin documentation_clean

.PHONY: cleanobjs
cleanobjs:
	-@find . -type f -name '*.o' -delete
	-@$(RM) -f $(path_sw2lib)/$(lib_sw2)

.PHONY: cleanbin
cleanbin:
	-@rm -f stepwat
	-@rm -f stepwat_test
	-@rm -f testing.sagebrush.master/stepwat
	-@rm -f testing.sagebrush.master/Stepwat_Inputs/stepwat

.PHONY: output_clean
output_clean:
	-@rm -rf testing.sagebrush.master/Output/*
	-@rm -rf testing.sagebrush.master/Stepwat_Inputs/Output/*

.PHONY : documentation_clean
documentation_clean :
		@rm -rf Documentation/html

.PHONY : documentation
documentation:
		@doxygen doxyfile
		@if open Documentation/html/index.html; \
		  then echo "Success"; \
		  else if echo "Open failed. Attempting fallback method." && xdg-open Documentation/html/index.html; \
		    then echo "Success"; \
		    else echo "Failed to open documentation"; \
		  fi; \
		fi
