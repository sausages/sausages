# no-missing-field-initializers means we allow struct s={} to initialise
# members of s to zero/null
# Make sure not to include any -O, even -O0, as even this can hide loop counter
CXXFLAGS = -g -Wall -Wextra -Wno-missing-field-initializers

HEADERS = point.h io.h find_sausages.h
OBJS    = main.o io.o find_sausages.o point.o

sausages: $(HEADERS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

.PHONY:
clean:
	rm -f *.o sausages
