CXXFLAGS = -g -O0 -Wall -Wextra

HEADERS = main.h io.h
OBJS    = main.o io.o

sausages: $(HEADERS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

.PHONY:
clean:
	rm -f *.o sausages
