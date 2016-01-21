import numpy as np
from tailseqext import PolyALocator

def locate_polya_seqbased(query, max_term_mod, weight):
    H = np.zeros([len(query), len(query)], np.int32)

    for i in range(max_term_mod): # 0-based, left-inclusive
        for j, jbase in enumerate(query): # 0-based, right-inclusive
            if i == j:
                H[i, j] = weight[jbase]
            elif i < j:
                H[i, j] = H[i, j-1] + weight[jbase]

    # find the longest stretch of T
    row, col = np.where(H >= 1)
    if len(row) == 0:
        return 0, 0

    longest = np.argmax(col - row)
    tstart, tend = row[longest], col[longest]

    # trim non-T flankings
    while query[tstart] != 'T' and tstart < tend:
        tstart += 1
    while query[tend] != 'T' and tend > tstart:
        tend -= 1

    return tstart, tend + 1


from tailseq.parsers import parse_sqi


weights = {'T': 1, 'G': -10, 'A': -10, 'C': -10, 'N': -5}
loc = PolyALocator(weights)

for spot in parse_sqi('../sequences/fifth_3T3.lint_sqi.bgz'):
    seq = spot.seq[77:]
    pyresult = locate_polya_seqbased(seq, 30, weights)
    cresult = loc(seq, 30)
    if pyresult != cresult:
        print 'Differs for', seq
        print 'Python:', pyresult
        print 'C:', cresult

print 'Done'
