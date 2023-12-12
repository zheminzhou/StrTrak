#!/usr/bin/env python

import sys, re

refLen = {}

if __name__ == '__main__' :
    cutoff = float(sys.argv[1])
    cin = sys.stdin if len(sys.argv) < 3 else open(sys.argv[2])
    
    for line in cin :
        if line.startswith('@') :
            if line.startswith('@SQ') :
                p = line.strip().split()
                refLen[p[1][3:]] = int(p[2][3:])
            sys.stdout.write(line)
            continue

        part = line.strip().split('\t', 10)
        flag = int(part[1])
        if flag & 4 > 0 :
            continue
        try :
            score = int(re.findall('AS:i:(\d+)', line)[0])
        except :
            continue

        flags = re.findall('(\d+)([MDSH])', part[5])
        gap = [int(flags[0][0]) if flags[0][1] in ('S', 'H') else 0, int(flags[-1][0]) if flags[-1][1] in ('S', 'H') else 0]
        if 'H' in part[5] :
            aln = len(part[9])
        else :
            aln = len(part[9]) - sum(gap)

        ref_aln = sum([ int(flag[0]) for flag in flags if flag[1] in ('M', 'D') ])
        start_gap = min( int(part[3])-1, gap[0] )
        end_gap = min( refLen[part[2]]-(int(part[3])+ref_aln-1), gap[1] )
        gap = [start_gap, end_gap]

        if min(gap) >= 10 :
            continue
        if sum(gap) >= 0.25 * aln and aln < 70 :
            continue
        dist = (aln*2. - score)/6.
        if dist >= aln*cutoff :
            continue
        sys.stdout.write(line)
        sys.stdout.flush()
