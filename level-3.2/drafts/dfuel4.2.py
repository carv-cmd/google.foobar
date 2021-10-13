from pprint import pprint
from fractions import Fraction


def gauss_redux(matrix, trange):
    _range = range
    nrows = len(matrix)
    qidentity = matrix_identity(trange)
    diff_matrix = [[qidentity[y][x] - matrix[y][x] for x in trange] for y in trange]
    inv = [diff_matrix[tx] + qidentity[tx] for tx in trange]
    for ix in _range(nrows):
        weight = 0
        if inv[ix][ix] == 0:
            break
        for jx in _range(nrows):
            if jx != ix:
                weight = inv[jx][ix] / inv[ix][ix]
            for kx in _range(nrows * 2):
                inv[jx][kx] = inv[jx][kx] - weight * inv[ix][kx]
    return [[inv[i][x] / inv[i][i] for x in _range(nrows, nrows * 2)] for i in _range(nrows)]

def one_stop(prow):
    print prow
    head = Fraction(1) - prow[0]
    head = Fraction(head.denominator, head.numerator)
    pvector = [head * ix for ix in prow[1:]]
    return matrix_slice(pvector)

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
    pvector = [int((num * (dmax // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def scaling_table(mx, trange):
    vv = [0 for _ in trange]

    # Iterate over rows tablulating row weight
    for ix in trange:
        weight = 0.0

        for jx in trange:
            setweight = abs(mx[ix][jx])
            weight = (setweight if setweight > weight else weight)

            if setweight > weight:
                weight = setweight

        if not weight:
            break

        vv[ix] = 1.0 / weight

    return vv


def crout_decomp(mx, trange):
    # int: i, imax, j, k
    # float: big(weight), dum, sum, temp
    # float *vv
    # returns: mx, ovekt, 
    nrows = len(mx)
    vvs = scaling_table(mx, trange)
    ovekt = [0 for _ in trange]
    tiny = 1.0**20
    parity = 1.0

    for j in trange:

        for i in range(j):
            sumi = mx[i][j]
            for k in range(i):
                sumi -= mx[i][j] * mx[k][j]
            mx[i][j] = sumi

        big = 0.0
        for ix in range(j, nrows):
            sumi = mx[ix][j]
            for kx in range(j):
                sumi -= mx[ix][kx] * mx[kx][j]
            mx[ix][j] = sumi
            dums = vvs[ix] * abs(sumi)
            if dums >= big:
                big = dums
                imax = ix

        if j != imax:
            for kxt in trange:
                dums = a[imax][kxt]
                mx[imax][kxt] = a[j][kxt]
            parity = not parity
            vvs[imax] = vvs[j]

        ovekt[j] = imax 
        if not mx[j][j]:
            mx[j][j] = tiny

        if j != nrows:
            dums = 1.0 / mx[j][j]
            for irx in range(j + 1, nrows):
                mx[irx][j] *= dums
    
    return mx, ovekt,parity 


def crout_inv(mx, trange):
    mx, indx, dp = crout_decomp(mx, trange)
    itm = 0
    bx = [0 for _ in trange]

    for i in trange:
        ips = indx[i]
        sumi = bx[ips]
        bx[ips] = bx[i]
        if itm:
            for j in range(itm, i - 1):
                sumi -= mx[i][j] * bx[j]
        elif sumi:
            itm = i
        bx[i] = sumi

    for rx in reversed(trange):
        sumi = bx[rx]
        for jx in range(rx + 1):
            sumi -= mx[rx][jx] * bx[jx]
        bx[rx] = sumi / mx[rx][rx]

    print '\n\nbx: ', bx

def main_loop(txr, trange, roffset):
    # Transient(Q) & Absorbing(R) probability matrices must be isolated
    # Canonical form essentially given; split matrix into respective Q/R matrices
    Q_matrix = [[txr[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[txr[nyt][nxt] for nxt in roffset] for nyt in trange]
    _foo = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in Q_matrix]
    print '\nQ-Matrix:'; pprint(_foo, width=len(_foo[0]) * 10)

    # Solve for; ( ( identity(Q) - Q_matrix )^-1 * R )
    qidentity = matrix_identity(trange)
    imq_matrix = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]
    diffs = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in imq_matrix]
    print '\nInverted:'; pprint(diffs, width=120)

    iq_inverse = crout_inv(imq_matrix, trange)

    if True: return
    multi = matrix_multiply(iq_inverse, R_matrix)[0]

    # Slice and return relevant probabilities
    return matrix_slice(multi)


def solution(m):
    # Extract non-redundant rows; count transient and absorbing states.
    datarows = filter(any, m)
    ntransits = len(datarows)
    nabsorbs = len(m) - ntransits

    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    trange = xrange(ntransits)
    roffset = xrange(ntransits, len(m))

    # Normalize each rowspan into fractions of that row total.
    datarows = normalize_rowspans(datarows)

    fuk = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in datarows]

    print '\nDataRows'; pprint(fuk, width=120)

    if nabsorbs == 1:
        oh_markov = [1, 1]
    elif ntransits == 1:
        oh_markov = one_stop(datarows[0])
    else:
        oh_markov = main_loop(datarows, trange, roffset)

    return
    #return oh_markov


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
            [0, 1, 0, 3, 1, 1, 5, 1, 0, 4], 
            [3, 0, 0, 3, 0, 2, 0, 4, 1, 0], 
            [1, 2, 0, 3, 5, 1, 3, 1, 0, 0], 
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
        print '\n>> MARKOV:'; pprint(s, width=len(s) * 8)
        fok = solution(s)
        print '\n>>> RETURNED: {}'.format(fok)


