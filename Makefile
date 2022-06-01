CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=#-static# -lm
PROGRAM=qTMclust USalign USalignLibTests TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta pdb2ss NWalign HwRMSD cif2pdb
COMMON_DEPS=basic_fun.cpp

all: ${PROGRAM}

qTMclust: ${COMMON_DEPS} qTMclust.cpp HwRMSD.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

USalign: ${COMMON_DEPS} USalign.cpp USalignMain.cpp MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} USalign.cpp USalignMain.cpp -o $@ ${LDFLAGS}

USalignLibTests: ${COMMON_DEPS} USalign.cpp USalignLib.cpp test/USalignLibTests.cpp ${COMMON_DEPS}
	${CC} ${CFLAGS} ${COMMON_DEPS} USalign.cpp USalignLib.cpp test/USalignLibTests.cpp -o $@ ${LDFLAGS}

TMalign: ${COMMON_DEPS} TMalign.cpp param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

TMscore: ${COMMON_DEPS} TMscore.cpp TMscore.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

MMalign: ${COMMON_DEPS} MMalign.cpp MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

se: ${COMMON_DEPS} se.cpp se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

pdb2ss: ${COMMON_DEPS} pdb2ss.cpp se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

pdb2xyz: ${COMMON_DEPS} pdb2xyz.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

xyz_sfetch: xyz_sfetch.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2fasta: ${COMMON_DEPS} pdb2fasta.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

NWalign: ${COMMON_DEPS} NWalign.cpp NWalign.h basic_fun.h pstream.h BLOSUM.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

HwRMSD: ${COMMON_DEPS} HwRMSD.cpp HwRMSD.h NWalign.h BLOSUM.h se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h
	${CC} ${CFLAGS} ${COMMON_DEPS} $@.cpp -o $@ ${LDFLAGS}

cif2pdb: cif2pdb.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM}
