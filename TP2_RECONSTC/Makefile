CC=gcc
CFLAGS=-Wall 
LDFLAGS= $(CFLAGS)
LDLIBS=-lm
RM=rm -f

# ajouter ci-dessous le nom du nouvel objet
OBJECTS=limace.o salade.o fonc.o reconst.o

# ajouter ci-dessous le nom du nouvel executable
EXE=reconst

.PHONY: all
all: $(EXE)

COMMONOBJ=limace.o salade.o fonc.o

reconst: $(COMMONOBJ) reconst.o 
# ajouter ici une ligne pour le nouvel executable

limace.o: limace.h
salade.o: limace.h salade.h
fonc.o: limace.h salade.h fonc.h
reconst.o: limace.h salade.h fonc.h
# ajouter ici une ligne pour le nouvel objet


clean:
	$(RM) $(OBJECTS) $(EXE)
