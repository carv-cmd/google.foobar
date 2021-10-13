from pprint import pprint
from fractions import Fraction


def gauss02(tmtx, amtx):
    _range = range
    _frac = Fraction
    nrows = len(tmtx)
    trange = _range(nrows)
    arange = _range(len(amtx[0]))
    # Record pivot statistics
    audit = {
            'r': list(trange), 
            'c': list(trange), 
            'p': [_frac(0) for _ in trange]
            }
    for mains in trange:
        weight = _frac(0, 1)
        for jx in trange:
            if audit['p'][jx] != 1:
                for ky in trange:
                    sabs = abs(tmtx[jx][ky])
                    if (audit['p'][ky] == 0) and (sabs >= weight):
                        weight, irow, icol = sabs, jx, ky
        audit['p'][icol] += _frac(1, 1)
        if irow != icol:
            for mx in (tmtx, amtx):
                for ex in range(len(mx[0])):
                    mx[irow][ex], mx[icol][ex] = mx[icol][ex], mx[irow][ex]
        audit['r'][mains] = irow
        audit['c'][mains] = icol
        pinvt = _frac(1, 1) / tmtx[icol][icol]
        tmtx[icol][icol] = _frac(1, 1)
        for mx in (tmtx, amtx):
            for ex in range(len(mx[0])):
                mx[icol][ex] *= pinvt
        for rix in trange:
            if rix != icol:
                redux = tmtx[rix][icol]
                tmtx[rix][icol] = _frac(0, 1)
                for mx in (tmtx, amtx):
                    for ex in range(len(mx[0])):
                        mx[rix][ex] -= mx[icol][ex] * redux
    for flip in reversed(trange):
        if audit['r'][flip] != audit['c'][flip]:
            for kips in trange:
                idr, idc = audit['r'][flip], audit['c'][flip]
                tmtx[kips][idr], tmtx[kips][idc] = tmtx[kips][idr], tmtx[kips][idc]

    #print '\n'; pprint(audit, width=120)
    #for foobaz in (tmtx, amtx):
    #    seeinv = [['{}/{}'.format(my.numerator, my.denominator) for my in mx] 
    #            for mx in foobaz]
    #    print '\n::: (a/t) final :::'; pprint(seeinv, width=len(foobaz[0] * 20))

    return amtx


def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix


def matrix_identity(idlen):
    # Generate identity matrices
    fracs = Fraction
    return [[(fracs(0) if x != y else fracs(1)) for x in idlen] for y in idlen]

def gcd(a, b):
    if (b == 0):
        return a
    else:
        print a, b
        return gcd(b, a % b)

def get_gcd(shrink):
    L = len(shrink)
    out = 0
    for i in range(len(shrink)):
        out = gcd(out, shrink[i])
    return out



def matrix_slice(m0slice):
    dmax = get_gcd(m0slice).denominator
    sliced = [(state.numerator, state.denominator) for state in m0slice]
    pvector = [int((num * (dmax // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector


def main_loop(datums):
    # Count transient and absorbing states.
    ntransits = len(datums)
    nstates = len(datums[0])
    nabsorbs = nstates - ntransits

    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    trange = range(ntransits)
    aoffset = range(ntransits, nstates)

    # Isolating transient(Q) & absorbing(R) probability matrices
    Q_matrix = [[datums[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[datums[nyt][nxt] for nxt in aoffset] for nyt in trange]

    qidentity = matrix_identity(trange)
    gausser = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]

    #for foobaz in (gausser, R_matrix):
    #    seeinv = [['{}/{}'.format(my.numerator, my.denominator) for my in mx] 
    #            for mx in foobaz]
    #    print '\n::: (Q/R) matrix :::'; pprint(seeinv, width=len(foobaz[0] * 10))

    # Solve for; ( ( identity(Q) - Q_matrix )^-1 * R )
    # Slice and return relevant probabilities
    iq_inverse = gauss02(gausser, R_matrix)

    return matrix_slice(iq_inverse[0])
    #return matrix_slice(foos)

def solution(m):
    # Extract non-redundant rows
    datarows = filter(any, m)

    # Normalize each rowspan into fractions of that row total.
    datarows = normalize_rowspans(datarows)
    #fuck = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in datarows]
    #print '\n::: INITIAL :::'; pprint(fuck, width=len(datarows[0]) * 10); print

    if len(m) == 1:
        oh_markov = [1, 1]
    else:
        oh_markov = main_loop(datarows)

    return oh_markov


if __name__ == "__main__":
    # ( ::: 5x5 ::: t(2) ::: a(3) ::: [7, 6, 8, 21] ::: )
    sol01 = [ 
            [0, 2, 1, 0, 0], 
            [0, 0, 0, 3, 4], 
	    [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0],
            ]

    # ( ::: 5x5 ::: t(2) ::: a(4) ::: [0, 3, 2, 9, 14] ::: )
    sol02 = [ 
            [0, 1, 0, 0, 0, 1], 
            [4, 0, 0, 3, 2, 0], 
	    [0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0], 
	    [0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0],
            ]

    # ( ::: 10x10 ::: t(9) ::: a(1) ::: [1, 1] ::: )
    sol03 = [ 
            [1, 1, 1, 3, 4, 1, 1, 2], 
            [3, 1, 3, 3, 2, 3, 0, 4], 
            [0, 0, 0, 1, 0, 0, 8, 9], 
            [0, 0, 3, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0], 
            ]
    
    # ( ::: 5x5 ::: t(4) ::: a(1) ::: [] ::: )
    sol04 = [ 
            [9, 2, 1, 5, 0], 
            [1, 7, 2, 3, 0], 
	    [2, 1, 8, 2, 0], 
            [2, 4, 3, 4, 1], 
            [0, 0, 0, 0, 0],
            ]

    # ( ::: 5x5 ::: t(1) ::: a(4) ::: [] ::: )
    sol05 = [ 
            [5, 2, 1, 1, 2], 
            [0, 0, 0, 0, 0], 
	    [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0],
            ]

    # ( ::: 5x5 ::: t(1) ::: a(4) ::: [] ::: )
    sol06 = [ 
            [0, 2, 1, 5, 3], 
            [0, 0, 0, 0, 0], 
	    [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0],
            ]

    # ( ::: 2x2 ::: t(1) ::: a(1) ::: [1, 1] ::: )
    sol07 = [ 
            [0, 2], 
            [0, 0],
            ]

    # ( ::: 2x2 ::: t(1) ::: a(1) ::: [1, 1] ::: )
    sol08 = [[1]]

    # ( ::: 4x4 ::: t(3) ::: a(1) ::: [1, 1] ::: )
    sol09 = [
            [2, 0, 1, 0], 
            [1, 2, 2, 1], 
            [2, 2, 1, 0], 
            [0, 0, 0, 0], 
            ]
    
    # ( ::: 2x2 ::: t(1) ::: a(1) ::: [1, 1] ::: )
    sol10 = [
            [4, 2, 3, 2], 
            [0, 0, 0, 0], 
            [0, 0, 0, 0], 
            [0, 0, 0, 0], 
            ]

    sol33 = [ # [1, 2, 3]
            [1, 2, 3, 0, 0, 0],
            [4, 5, 6, 0, 0, 0], 
            [7, 8, 9, 1, 0, 0], 
            [0, 0, 0, 0, 1, 2], 
            [0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0]
            ]

    sol55 =[ # [1, 2, 3, 4, 5, 15]
            [0, 0, 12, 0, 15, 0, 0, 0, 1, 8], 
            [0, 0, 60, 0, 0, 7, 13, 0, 0, 0], 
            [0, 15, 0, 8, 7, 0, 0, 1, 9, 0], 
            [23, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
            [37, 35, 0, 0, 0, 0, 3, 21, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            ]

    _testing = (sol01, sol02, sol03, sol04, sol05, sol06, sol07, sol08, sol09, sol10)

    for s in (sol33, sol55):  #, sol06, sol10):
        #print '\n>> MARKOV:'; pprint(s, width=len(s) * 8)
        fok = solution(s)
        print '\n>>> RETURNED: {}'.format(fok)
        #print '\n>>>> returned: type({})'.format(type(fok[0]))


