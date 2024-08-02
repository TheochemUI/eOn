OBJECTS := Tip4p_Pt.o zhu_philpott.o zhu_philpott_parameters.o
CXXFLAGS += -I..
LIB := libtip4p_pt.a
$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
.PHONY : clean
