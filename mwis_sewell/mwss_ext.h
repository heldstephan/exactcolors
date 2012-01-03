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

/* return value i fa timeout occured: */
#define SEWELL_TIMEOUT 10

/** Compute the maximum-weight stable set for a simple graph.
    @remark: it is required that the input graph is simple, i.e. does not contain parallel edges!
    @param newset: when succeeding newset will point to the array
                   of vertices of an maximum-weight stable set.
    @param nnewset: when succeeding nnewset will point to the
                    number of vertices in the found maximum-weight
                    stable set.
    @param ncount: number of vertices in the input graph.
    @param ecount: number of edges in the input graph.
    @param elist:  list of edges, given as vertex-pairs, of the input graph.
                   This means that elist contains 2*ecount vertex-entries and
                    edge number i connects the vertices elist[2i] and elist[2i+1].
    @param nweights: weights of the vertices.
    @param lower_bound: initial lower bound for the objective. A good initial
                    lower bound can speep-up the computation significantly.
                    If not good bounds are present this parameter should be set to 0.
    @param goal: If a stable set of weight @c goal or greater is found the algorithm stops.
                 This parameter should be set to MWISNW_MAX if the maximum-weight stable
                 set should be found.
    @return      0 iff function completed successfully,
                 SEWELL_TIMEOUT if a timeout occured.
*/

extern
int SEWELL_optimize(int** newset,
                    int*  nnewset,
                    int   ncount, int ecount, const int elist[], NWT nweights[],
                    NWT   lower_bound,
                    NWT   goal);

/** Same as SEWELL_optimize, but with the ability to limit the running time, thus
    serving as a heuristic only.

    @param cpu_limit limit for the number of cpu seconds spend in branch & bound.
                     A negative value meand "no limit".
*/
extern
int SEWELL_heur(int** newset,
                int*  nnewset,
                int   ncount, int ecount, const int elist[], NWT nweights[],
                NWT   lower_bound,
                NWT   goal,
                double cpu_limit);

extern
int SEWELL_node_limit(void);

#endif
