MKOCT_FLAGS = $(shell pkg-config --libs cfitsio --cflags cfitsio)

all: read_fits_image.oct save_fits_image.oct save_fits_image_multi_ext.oct

clean:
	rm -f *.oct *.o

%.oct: %.cc
	mkoctfile --verbose -Wall $(MKOCT_FLAGS) $<

