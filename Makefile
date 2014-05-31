CXXFLAGS = -g -O0 -Wall -Wextra

HEADERS = main.h io.h find_sausages.h
OBJS    = main.o io.o find_sausages.o

sausages: $(HEADERS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

.PHONY:
clean:
	rm -f *.o sausages
