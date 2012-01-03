#     This file is part of exactcolors.

#     exactcolors is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     exactcolors is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.



# NOTE: CFLAGS are passed here from Makefile in parent directory
# CFLAGS=  -g
# CFLAGS += -std=c99 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings

all: libsewell.a sewell

clean:
	rm -f *.o *.a sewell

sewell: libsewell.a mwss.c
	$(CC) $(CFLAGS) -c mwss.c -o mwss.o
	$(CC) -o sewell mwss.o wstable.o -L. -lsewell

libsewell.a: wstable.c  mwss_ext.c mwss.h mwss_ext.h
	$(CC) $(CFLAGS) -c mwss_ext.c -o mwss_ext.o
	$(CC) $(CFLAGS) -c wstable.c -o wstable.o
	$(AR) rcs libsewell.a wstable.o mwss_ext.o
