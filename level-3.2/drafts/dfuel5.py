from pprint import pprint
from fractions import Fraction


def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mrows = xrange(len(mtx0))
    mcols = xrange(len(mtx1[0]))
    return [[sum([mtx0[spans][idy] * mtx1[idy][idx] 
        for idy in mrows]) 
        for idx in mcols] 
        for spans in mrows]

def matrix_identity(idlen):
    # Generate identity matrices
    fracs = Fraction
    return [[(fracs(0) if x != y else fracs(1)) for x in idlen] for y in idlen]

def matrix_slice(m2slice):
    sliced = [(state.numerator, state.denominator) for state in m2slice]
    dmax = int(max((maxi[1] for maxi in sliced)))
    #pvector = [(int(num * (dmax // den)) if den != dmax else int(num)) for num, den in sliced]
    pvector = [int((num * (dmax // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector

def gauss_redux0(mtx):
    seeinv = [['{}/{}'.format(my.numerator, my.denominator) for my in mx] for mx in mtx]
    print 'PreGauss'; pprint(seeinv, width=len(mtx[0] * 10))

    _range = range
    _frac = Fraction
    nrows = len(mtx)
    trange = _range(nrows)
    for ix in trange:
        weight = _frac(0, 1)
        if mtx[ix][ix] == 0:
            break
        for jx in trange:
            if jx != ix:
                weight = mtx[jx][ix] / mtx[ix][ix]

            for kx in _range(nrows * 2):
                mtx[jx][kx] = mtx[jx][kx] - weight * mtx[ix][kx]

    invert = []
    invappend = invert.append
    for ivt in trange:
        invappend([mtx[ivt][xi] / mtx[ivt][ivt] for xi in _range(nrows, nrows * 2)])

    seeinv = [['{}/{}'.format(my.numerator, my.denominator) for my in mx] for mx in invert]
    print '\nPostGaussInvert'; pprint(seeinv, width=len(invert[0] * 10))
    return invert

def gauss_redux(qam, rbm):
    #mview = [['{}/{}'.format(my.numerator, my.denominator) for my in mx] for mx in matrix]
    #print '\nPreGauss'; pprint(mview, width=120)
    _range = range
    _frac = Fraction
    nrows = len(qam)
    trange = _range(nrows)
    qrange = _range(len(qam[0]))

    idxrow = list(trange)
    idxcol = list(trange)
    pivots = [0 for _ in trange]
    irow, icol = 0, 0

    for ic in trange:
        big = _frac(0, 1)
        for pchk in trange:
            if pivots[pchk] != 1:
                for ksub in trange:
                    if not pivots[ksub]:
                        chks = abs(qam[pchk][ksub])
                        if chks >= big:
                            big = chks
                            irow = pchk
                            icol = ksub

        pivots[icol] += _frac(1, 1)

        if irow != icol:
            # SWAP(a, b): a, b = b, a
            for runs in qrange:
                qam[irow][runs], qam[icol][runs] =  qam[icol][runs], qam[irow][runs]
            for crun in trange:
                rbm[irow][runs], rbm[icol][runs] =  rbm[icol][runs], rbm[irow][runs]

        idxrow[ic] = irow
        idxcol[ic] = icol

        if qam[icol][irow] == 0:
            print 'singular matrix'
            return [0, 0]

        pivinv = _frac(1, 1) / qam[icol][irow]
        qam[icol][icol] = _frac(1, 1)

        for ipn in trange:
            qam[icol][ipn] *= pivinv
        for ipm in trange:
            rbm[icol][ipm] *= pivinv

        for lx in trange:
            if lx != icol:
                redux = qam[lx][icol]
                qam[lx][icol] = _frac(0, 1)
                for lnx in trange:
                    qam[lnx][lx] -= qam[icol][lx] * redux
                for lmx in trange:
                    rbm[lmx][lx] -= qam[icol][lx] * redux


    for flip in reversed(trange):
        if idxrow[flip] != idxcol[flip]:
            for kips in trange:
                qam[kips][idxrow[flip]], qam[kips][idxcol[flip]] = \
                        qam[kips][idxrow[flip]], qam[kips][idxcol[flip]]

    #return fundamental
    #pprint(qam, width=120)
    #pprint(rbm, width=120)
    return rbm


def main_loop(datums):
    # Count transient and absorbing states.
    ntransits, nstates = len(datums), len(datums[0])
    nabsorbs = nstates - ntransits

    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    trange = range(ntransits)
    aoffset = range(ntransits, nstates)

    # Isolating transient(Q) & absorbing(R) probability matrices
    Q_matrix = [[datums[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[datums[nyt][nxt] for nxt in aoffset] for nyt in trange]

    qidentity = matrix_identity(trange)
    gauss_prep = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]
    #gauss_id = [gauss_prep[tx] + qidentity[tx] for tx in trange]

    # Solve for; ( ( identity(Q) - Q_matrix )^-1 * R )
    #print gauss_redux(Q_matrix, trange)
    #iq_inverse = gauss_redux0(gauss_id)
    iq_inverse = gauss_redux(gauss_prep, Q_matrix)

    print '\nGauss Return: ', iq_inverse

    multi = matrix_multiply(iq_inverse, R_matrix)[0]
    return multi

    # Slice and return relevant probabilities
    #return matrix_slice(multi)

def solution(m):
    # Extract non-redundant rows; count transient and absorbing states.
    datarows = filter(any, m)
    # Normalize each rowspan into fractions of that row total.
    datarows = normalize_rowspans(datarows)
    fuck = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in datarows]
    print '\n::: INITIAL :::'; pprint(fuck, width=len(datarows[0]) * 10); print

    nt = len(datarows)
    na = len(m) - nt

    if nt - na == 1:
        oh_markov = [1, 1]
    #elif (nt == 1) and (datarows[0][0] == 0):
    #    oh_markov = matrix_slice(datarows[0][1:])
    else:
        #oh_markov = main_loop(datarows, trange, roffset)
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
            [3, 1, 2, 4, 1, 1, 2, 0, 0, 0], 
            [0, 4, 6, 3, 0, 0, 0, 1, 4, 2], 
            [0, 1, 0, 3, 1, 1, 5, 1, 0, 7], 
            [3, 0, 0, 3, 0, 2, 0, 4, 1, 0], 
            [1, 2, 0, 3, 7, 1, 3, 1, 0, 0], 
            [0, 1, 0, 3, 0, 1, 0, 2, 4, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
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
    sol08 = [
            [9, 2], 
            [0, 0], 
            ]

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

    _testing = (sol01, sol02, sol03, sol04, sol05, sol06, sol07, sol08, sol09, sol10)

    for s in (sol01, sol02):  #, sol06, sol10):
        #print '\n>> MARKOV:'; pprint(s, width=len(s) * 8)
        fok = solution(s)
        print '\n>>> RETURNED: {}'.format(fok)
        #print '\n>>>> returned: type({})'.format(type(fok[0]))


