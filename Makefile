# no-missing-field-initializers means we allow struct s={} to initialise
# members of s to zero/null
# Make sure not to include any -O, even -O0, as even this can hide loop counter
CXXFLAGS = -g -Wall -Wextra -Wno-missing-field-initializers

SRCDIR   = src
BUILDDIR = build

FILES   = sausages.cpp io.cpp main.cpp point.cpp
SRC     = $(patsubst %, $(SRCDIR)/%, $(FILES))
HEADERS = $(patsubst %.cpp, $(SRCDIR)/%.h, $(FILES))
OBJS    = $(patsubst %.cpp, $(BUILDDIR)/%.o, $(FILES))

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(SRCDIR)/%.h
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c -o $@ $<

sausages: $(HEADERS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

.PHONY:
clean:
	rm -f build/* sausages
