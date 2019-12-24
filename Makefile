ODIR = obj
SDIR = src
BDIR = bin
IDIR = include

OBJ = main.o rrt.o pgm.o
      
BIN = rrt

CFlags=-c -O3 -Wall -fopenmp 
LDFlags= -lm -fopenmp 
CC=g++
RM=rm

all: $(OBJ) $(BIN) 

$(BIN): $(patsubst %,$(ODIR)/%,$(OBJ))
	@mkdir -p $(BDIR)
	$(CC) $^ $(LDFlags) -o $(BDIR)/$@ 

%.o: $(SDIR)/%.c	
	@mkdir -p $(ODIR)
	$(CC) $(CFlags) -I$(IDIR) -o $(ODIR)/$@ $<

%.o: $(SDIR)/%.cpp	
	@mkdir -p $(ODIR)
	$(CC) $(CFlags) -I$(IDIR) -o $(ODIR)/$@ $<
	
clean:
	$(RM) $(patsubst %,$(ODIR)/%,$(OBJ)) $(BDIR)/$(BIN) 
	

