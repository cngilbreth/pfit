include $(HOME)/Projects/make.inc

FC    :=  $(FC)

################################################################################
# Targets
################################################################################

BUILD = build
MOD = mod
BIN = ../bin

$(shell mkdir -p $(MOD) $(BUILD) $(BIN))
vpath %.o $(BUILD)
vpath %.mod $(MOD)
vpath %.f90 ./
vpath pfit $(BIN)

OBJ_PFIT = util.o main.o options.o math.o

all: pfit

pfit: $(OBJ_PFIT)
	$(FC) -o $(BIN)/$@ $+ $(LAPACK95) $(MKL)

clean: 
	rm -f *.o *.mod build/* mod/* pfit 

################################################################################
# File dependencies
################################################################################

$(BUILD)/main.o: math.o util.o options.o
$(BUILD)/options.o: util.o
$(BUILD)/util.o: 
$(BUILD)/math.o: 


################################################################################
# Generic build rules
################################################################################

$(BUILD)/%.o: %.cpp
	$(CXX) -c $< -o $@
$(BUILD)/%.o: %.c
	$(CC) -c $< -o $@
$(BUILD)/%.o: %.f90
	$(FC) -c $< -o $@
$(BUILD)/%.o: %.f
	$(FC) -c $< -o $@
$(BUILD)/%.o: %.F90
	$(FC) -c $< -o $@
$(BUILD)/%.o: %.F
	$(FC) -c $< -o $@
