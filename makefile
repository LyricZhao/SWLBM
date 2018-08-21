TARGET = LbmCavity3D
USER = $(shell whoami)
CC = sw5cc
LD = mpicc

CFLAGS = -O3 -I/usr/sw-mpp/mpi2/include/ -lm -msimd -OPT:alias=disjoint
CFLAGS_HOST = -host $(CFLAGS)
CFLAGS_SLAVE = -slave -ver 5.421-sw-589 $(CFLAGS)

LDFLAGS = -hybrid

OBJ = LbmCavity3D.o Parallel.o TComputing.o TComputingSlave.o
OBJ_N = LbmCavity3D.o Collide.o Parallel.o Stream.o

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

original: $(OBJ_N) makefile
	$(LD) $(OBJ_N) $(LIB) $(LDFLAGS) -o $(TARGET)
	rm $(OBJ_N)

test: TestCode
	$(CC) -O3 -msimd -o TestCode TestCode.c

run_test:
	make test
	bsub -I -b -q q_sw_cpc_1 -cgsp 64 -n 16 -np 4 -share_size 6500 -host_stack 500 -J Lyric_TestCode ./TestCode $(USER)

run_original:
	make original
	bsub -I -b -q q_sw_cpc_1 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J Lyric_OriginalCode ./LbmCavity3D $(USER)

run:
	make $(TARGET)
	bsub -I -b -q q_sw_cpc_1 -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J Lyric_SWLBM_Opt ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(TARGET) $(OBJ)
