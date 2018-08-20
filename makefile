TARGET = LbmCavity3D
USER = $(shell whoami)
CC = sw5cc
LD = mpicc 

CFLAGS =  -O3 -host -I/usr/sw-mpp/mpi2/include/ -lm 
CFLAGS_SLAVE =  -O3 -slave -I/usr/sw-mpp/mpi2/include/ -lm

OBJ = LbmCavity3D.o Parallel.o Work.o WorkSlave.o

LIB = lib/liblbm.a

$(TARGET): $(OBJ)
	$(LD) -hybrid -msimd -O3 $(OBJ) $(LIB) -o $(TARGET) 
	
Parallel.o:Parallel.c
	$(CC) $(CFLAGS) -c $<
LbmCavity3D.o:LbmCavity3D.c
	$(CC) $(CFLAGS) -c $<
Work.o:Work.c
	$(CC) $(CFLAGS) -c $<
WorkSlave.o:WorkSlave.c
	$(CC) -msimd $(CFLAGS_SLAVE) -c $<
run:
	bsub -I -b -q q_sw_cpc_1 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J test ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(TARGET) $(OBJ) 
	
