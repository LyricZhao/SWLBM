TARGET = LbmCavity3D
USER = $(shell whoami)
CC = sw5cc
LD = mpicc

CFLAGS = -O3 -I/usr/sw-mpp/mpi2/include/ -lm -msimd -OPT:alias=disjoint
CFLAGS_HOST = -host $(CFLAGS)
CFLAGS_SLAVE = -slave $(CFLAGS)

LDFLAGS = -hybrid

OBJ = LbmCavity3D.o Parallel.o TComputingSlave.o

LIB = lib/liblbm.a

$(TARGET): $(OBJ)
	$(LD) $(OBJ) $(LIB) $(LDFLAGS) -o $(TARGET)
	rm $(OBJ)

LbmCavity3D.o: LbmCavity3D.c
	$(CC) $(CFLAGS_HOST) -c LbmCavity3D.c

Parallel.o: LbmCavity3D.c
	$(CC) $(CFLAGS_HOST) -c Parallel.c

TComputing.o: LbmCavity3D.c
	$(CC) $(CFLAGS_HOST) -c TComputing.c

TComputingSlave.o: LbmCavity3D.c
	$(CC) $(CFLAGS_SLAVE) -c TComputingSlave.c

run: LbmCavity3D
	bsub -I -b -q q_sw_cpc_1 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J test ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(TARGET) $(OBJ)
