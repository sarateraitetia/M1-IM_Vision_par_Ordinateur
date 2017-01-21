CC=gcc
CFLAGS=-Wall 
LDFLAGS= $(CFLAGS)
LDLIBS=-lm
RM=rm -f

# ajouter ci-dessous le nom du nouvel objet
OBJECTS=limace.o salade.o fonc.o projet.o

# ajouter ci-dessous le nom du nouvel executable
EXE=projet

.PHONY: all
all: $(EXE)

COMMONOBJ=limace.o salade.o fonc.o

projet: $(COMMONOBJ) projet.o 
# ajouter ici une ligne pour le nouvel executable

limace.o: limace.h
salade.o: limace.h salade.h
fonc.o: limace.h salade.h fonc.h
projet.o: limace.h salade.h fonc.h
# ajouter ici une ligne pour le nouvel objet


clean:
	$(RM) $(OBJECTS) $(EXE)
