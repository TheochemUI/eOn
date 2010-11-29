OBJECTS := ccl.o potential_base.o spce_ccl.o tip4p.o zhu_philpott_for_eon.o \
   zhu_philpott.o zhu_philpott_parameters.o
CXXFLAGS += -I..
LIB := libplatinum-water.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
.PHONY : clean
