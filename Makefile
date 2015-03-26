# no-missing-field-initializers means we allow struct s={} to initialise
# members of s to zero/null
# Make sure not to include any -O, even -O0, as even this can hide loop counter
CFLAGS =  -g -Wall
#CFLAGS = -Ofast
CXXFLAGS = -std=c++11 $(CFLAGS)

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
	$(CC) $(CFLAGS) -c -o $@ $<


.PHONY:
clean:
	rm -f build/* sausages
