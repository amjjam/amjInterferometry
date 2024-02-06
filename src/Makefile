prefix=/opt/amj

ifndef prefix_amjCom
	prefix_amjCom=$(prefix)
endif

LIBS=-lamjCom
LDFLAGS=-L$(prefix_amjCom)/lib

FLAGS=-Wall -g -fPIC -I$(prefix_amjCom)/include
CPPFLAGS=$(FLAGS)

AMJINTERFEROMETRY_SRC=amjInterferometry.C
AMJINTERFEROMETRY_OBJ=$(AMJINTERFEROMETRY_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)

all: libamjInterferometry.so

libamjInterferometry.so: $(AMJINTERFEROMETRY_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^ $(LDFLAGS) $(LIBS)

install: libamjInterferometry.so ../include/amjInterferometry.H
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	cp libamjInterferometry.so $(prefix)/lib
	cp ../include/amjInterferometry.H $(prefix)/include

uninstall:
	- rm $(prefix)/lib/libamjInterferometry.so
	- rm $(prefix)/include/amjInterferometry.H

clean:
	- rm $(AMJINTERFEROMETRY_OBJ) libamjInterferometry.so

