PYTHONVERSION = 2.7
CXX = g++
CFLAGS = -Wall -g -O2 -fopenmp -I./include -I/usr/include/python$(PYTHONVERSION)
LIBS = -lz -lgsl -lgslcblas -lpython$(PYTHONVERSION)
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
OUT = vlasov

default: all

$(OBJ_DIR)/mymath.o: $(INCLUDE_DIR)/mymath.h $(SRC_DIR)/mymath.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/mymath.cpp -o $(OBJ_DIR)/mymath.o

$(OBJ_DIR)/file_saving.o: $(INCLUDE_DIR)/file_saving.h $(SRC_DIR)/file_saving.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/file_saving.cpp -o $(OBJ_DIR)/file_saving.o

$(OBJ_DIR)/pfunc.o: $(INCLUDE_DIR)/pfunc.h $(SRC_DIR)/pfunc.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pfunc.cpp -o $(OBJ_DIR)/pfunc.o

$(OBJ_DIR)/pyinput.o: $(OBJ_DIR)/pfunc.o $(INCLUDE_DIR)/pyinput.h $(SRC_DIR)/pyinput.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pyinput.cpp -o $(OBJ_DIR)/pyinput.o

$(OBJ_DIR)/mesh.o: $(OBJ_DIR)/pyinput.o $(INCLUDE_DIR)/mesh.h $(SRC_DIR)/mesh.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/mesh.cpp -o $(OBJ_DIR)/mesh.o

$(OBJ_DIR)/fdtd.o: $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/mesh.o $(INCLUDE_DIR)/fdtd.h $(SRC_DIR)/fdtd.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/fdtd.cpp -o $(OBJ_DIR)/fdtd.o

$(OBJ_DIR)/pfc.o: $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/mesh.o $(INCLUDE_DIR)/pfc.h $(SRC_DIR)/pfc.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pfc.cpp -o $(OBJ_DIR)/pfc.o

$(OBJ_DIR)/plasmas.o: $(OBJ_DIR)/mymath.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mesh.o $(INCLUDE_DIR)/plasmas.h $(SRC_DIR)/plasmas.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/plasmas.cpp -o $(OBJ_DIR)/plasmas.o

$(OBJ_DIR)/particle.o: $(INCLUDE_DIR)/particle.h $(SRC_DIR)/particle.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/particle.cpp -o $(OBJ_DIR)/particle.o

$(OBJ_DIR)/testparticles.o: $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mesh.o $(INCLUDE_DIR)/testparticles.h $(SRC_DIR)/testparticles.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/testparticles.cpp -o $(OBJ_DIR)/testparticles.o

$(OBJ_DIR)/solver.o: $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/file_saving.o $(INCLUDE_DIR)/solver.h $(INCLUDE_DIR)/solver.h $(SRC_DIR)/solver.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/solver.cpp -o $(OBJ_DIR)/solver.o

$(OBJ_DIR)/output.o: $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/file_saving.o $(INCLUDE_DIR)/solver.h $(INCLUDE_DIR)/output.h $(SRC_DIR)/output.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/output.cpp -o $(OBJ_DIR)/output.o

$(BIN_DIR)/$(OUT): $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/plasmas.o $(OBJ_DIR)/particle.o $(OBJ_DIR)/testparticles.o $(OBJ_DIR)/solver.o $(OBJ_DIR)/output.o $(SRC_DIR)/main.cpp
	$(CXX) $(SRC_DIR)/main.cpp -o $(BIN_DIR)/$(OUT) $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/plasmas.o $(OBJ_DIR)/particle.o $(OBJ_DIR)/testparticles.o $(OBJ_DIR)/solver.o $(OBJ_DIR)/output.o $(CFLAGS) $(LFLAGS) $(LIBS)

all: $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/plasmas.o $(OBJ_DIR)/particle.o $(OBJ_DIR)/testparticles.o $(OBJ_DIR)/solver.o $(OBJ_DIR)/output.o $(BIN_DIR)/$(OUT)

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(OUT)
