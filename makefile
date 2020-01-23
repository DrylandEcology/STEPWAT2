CC  	=	gcc

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

INC_DIRS = \
	-I. \
	-Isqlite-amalgamation \
	-Isw_src \
	-Isw_src/googletest/googletest \
	-Isw_src/googletest/googletest/include \
	-Itest

LIBS = -lm

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
	sw_src/filefuncs.c \
	sw_src/generic.c \
	sw_src/mymemory.c \
	sw_src/pcg/pcg_basic.c \
	sw_src/rands.c \
	sw_src/SW_Carbon.c \
	sw_src/SW_Control.c \
	sw_src/SW_Files.c \
	sw_src/SW_Flow.c \
	sw_src/SW_Flow_lib.c \
	sw_src/SW_Markov.c \
	sw_src/SW_Model.c \
	sw_src/SW_Output.c \
	sw_src/SW_Output_get_functions.c \
	sw_src/SW_Output_outarray.c \
	sw_src/SW_Output_outtext.c \
	sw_src/SW_Site.c \
	sw_src/SW_Sky.c \
	sw_src/SW_SoilWater.c \
	sw_src/SW_VegEstab.c \
	sw_src/SW_VegProd.c \
	sw_src/SW_Weather.c \
	sw_src/Times.c \
	sxw.c \
	sxw_environs.c \
	sxw_resource.c \
	sxw_soilwat.c \
	sxw_sql.c \
	ST_initialization.c \
	ST_progressBar.c \
	ST_seedDispersal.c \
	ST_colonization.c

sources_test = \
	sw_src/googletest/googletest/src/gtest-all.cc \
	sw_src/googletest/googletest/src/gtest_main.cc \
	test/test_ST_mortality.cc

objects_core = $(sources_core:%.c=obj/%.o)
objects_core_test = $(sources_core:%.c=obj/%_TEST.o)
objects_test = $(sources_test:%.cc=obj/%.o)

all: stepwat stepwat_test

stepwat: $(objects_core)
	$(CC) $(objects_core) $(CFLAGS) $(CPPFLAGS) $(LIBS) -o stepwat
	-@cp stepwat testing.sagebrush.master
	-@cp stepwat testing.sagebrush.master/Stepwat_Inputs

stepwat_test: $(objects_core_test) $(objects_test)
	$(CXX) $(objects_core_test) $(objects_test) $(CFLAGS) $(CPPFLAGS) $(LIBS) -o stepwat_test
	-@cp stepwat_test testing.sagebrush.master
	-@cp stepwat_test testing.sagebrush.master/Stepwat_Inputs

obj/%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INC_DIRS) -c $< -o $@

obj/%.o: %.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(INC_DIRS) -std=gnu++11 -c $< -o $@

obj/%_TEST.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INC_DIRS) -DSTDEBUG -c $< -o $@

.PHONY: bint_testing_nongridded
bint_testing_nongridded: stepwat stepwat_test
	testing.sagebrush.master/Stepwat_Inputs/stepwat -d testing.sagebrush.master/Stepwat_Inputs -f files.in -o -i

.PHONY: bint_testing_gridded
bint_testing_gridded: stepwat stepwat_test
	testing.sagebrush.master/stepwat -d testing.sagebrush.master -f files.in -g

.PHONY: cleanall
cleanall: clean output_clean

.PHONY: clean
clean: cleanobjs cleanbin documentation_clean

.PHONY: cleanobjs
cleanobjs:
	-@find . -type f -name '*.o' -delete

.PHONY: cleanbin
cleanbin:
	-@rm -f stepwat
	-@rm -f stepwat_test
	-@rm -f testing.sagebrush.master/stepwat
	-@rm -f testing.sagebrush.master/stepwat_test
	-@rm -f testing.sagebrush.master/Stepwat_Inputs/stepwat
	-@rm -f testing.sagebrush.master/Stepwat_Inputs/stepwat_test

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
