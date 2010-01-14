/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2007 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/* This code may be used under the terms of the GNU General Public License  */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

#ifndef  __COLORSAFE_H
#define  __COLORSAFE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
extern int gethostname (char *, int);

#define COLOR_SBUFFER_SIZE (4000)
#define COLOR_SFNAME_SIZE (32)

#define COLOR_BOSS_SEND    'S'
#define COLOR_BOSS_RECEIVE 'R'
#define COLOR_BOSS_YES     'Y'
#define COLOR_BOSS_NO      'N'
#define COLOR_BOSS_EXIT    'X'
#define COLOR_BOSS_PORT ((unsigned short) 24870)

typedef struct COLOR_SFILE {
    int           status;
    int           desc;
    int           type;
    int           chars_in_buffer;
    int           current_buffer_char;     /* only used for reading */
    int           bits_in_last_char;       /* writing: number of empty bits in
                                            * buffer[chars_in_buffer];
                                            * reading: number of full bits in
                                            * buffer[?] */
    int           pos;
    char          fname[COLOR_SFNAME_SIZE];
    char          hname[COLOR_SFNAME_SIZE];
    unsigned char buffer[COLOR_SBUFFER_SIZE];
} COLOR_SFILE;

typedef struct COLOR_SPORT {
    unsigned short port;
    int t;
} COLOR_SPORT;


COLOR_SFILE
    *COLORsafe_sopen (const char *f, const char *s),
    *COLORsafe_snet_open (const char *hname, unsigned short p),
    *COLORsafe_snet_receive (COLOR_SPORT *s);

COLOR_SPORT
    *COLORsafe_snet_listen (unsigned short p);

int
    COLORsafe_swrite (COLOR_SFILE *f, char *buf, int size),
    COLORsafe_swrite_bits (COLOR_SFILE *f, int x, int xbits),
    COLORsafe_swrite_ubits (COLOR_SFILE *f, unsigned int x, int xbits),
    COLORsafe_swrite_char (COLOR_SFILE *f, char x),
    COLORsafe_swrite_string (COLOR_SFILE *f, const char *x),
    COLORsafe_swrite_short (COLOR_SFILE *f, short x),
    COLORsafe_swrite_ushort (COLOR_SFILE *f, unsigned short x),
    COLORsafe_swrite_int (COLOR_SFILE *f, int x),
    COLORsafe_swrite_uint (COLOR_SFILE *f, unsigned int x),
    COLORsafe_swrite_double (COLOR_SFILE *f, double x),
    COLORsafe_sread (COLOR_SFILE *f, char *buf, int size),
    COLORsafe_sread_bits (COLOR_SFILE *f, int *x, int xbits),
    COLORsafe_sread_ubits (COLOR_SFILE *f, unsigned int *x, int xbits),
    COLORsafe_sread_char (COLOR_SFILE *f, char *x),
    COLORsafe_sread_string (COLOR_SFILE *f, char *x, int maxlen),
    COLORsafe_sread_short (COLOR_SFILE *f, short *x),
    COLORsafe_sread_ushort (COLOR_SFILE *f, unsigned short *x),
    COLORsafe_sread_int (COLOR_SFILE *f, int *x),
    COLORsafe_sread_uint (COLOR_SFILE *f, unsigned int *x),
    COLORsafe_sread_double (COLOR_SFILE *f, double *x),
    COLORsafe_sflush (COLOR_SFILE *f),
    COLORsafe_sclose (COLOR_SFILE *f),
    COLORsafe_sbits (unsigned int x);

void
    COLORsafe_snet_unlisten (COLOR_SPORT *s);

#endif /* __COLORSAFE_H */
