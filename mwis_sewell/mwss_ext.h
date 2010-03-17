#ifndef __MWSS_EXT_H
#define __MWSS_EXT_H
/**
    This file is part of exactcolors.

    exactcolors is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exactcolors is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.
*/

typedef int NWT; /* node weight type */

extern
int SEWELL_optimize(int** newset,
                    int*  nnewset,
                    int   ncount, int ecount, const int elist[], NWT nweights[],
                    NWT   goal);

extern
int SEWELL_node_limit(void);

#endif
