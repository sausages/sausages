# no-missing-field-initializers means we allow struct s={} to initialise
# members of s to zero/null
# Make sure not to include any -O, even -O0, as even this can hide loop counter
CXXFLAGS = -g -Wall
#CXXFLAGS = -Ofast

SRCDIR   = src
BUILDDIR = build

FILES    = eigen-bits io main maths params point sausages
SRC      = $(patsubst %, $(SRCDIR)/%.cpp, $(FILES))
HEADERS  = $(patsubst %, $(SRCDIR)/%.h, $(FILES))
OBJS     = $(patsubst %, $(BUILDDIR)/%.o, $(FILES))

sausages: $(HEADERS) $(OBJS) $(BUILDDIR)/cJSON.o
	$(CXX) $(CXXFLAGS) $(OBJS) $(BUILDDIR)/cJSON.o -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(SRCDIR)/%.h
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILDDIR)/cJSON.o : $(SRCDIR)/cJSON/cJSON.c
	$(CC) $(CXXFLAGS) -c -o $@ $<


DOXYGEN_EXISTS  := $(shell command -v doxygen  2> /dev/null)

docs:
ifndef DOXYGEN_EXISTS
	$(error "Please install doxygen")
endif
	doxygen


.PHONY:
clean:
	rm -rf build sausages docs
