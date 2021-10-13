from pprint import pprint
from fractions import Fraction, gcd

def print_view(views):
    foo = [['{}/{}'.format(i.numerator, i.denominator) for i in ixi] for ixi in views]
    pprint(foo, width=len(views[0]) * 15)

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def matrix_slice(pvector, dmax=0):
    # Find common denominator in P_Matrix[0].
    for gx in pvector:
        dmax = (gcd(dmax, gx) if gx else dmax)
    dmax = dmax.denominator
    # Scale each numerator by common denominator and insert then append gcd.
    sliced = ((state.numerator, state.denominator) for state in pvector)
    pvector = [int((num * (int(dmax) // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector

def lu_decomp(mat):
    fract = Fraction
    _range = range
    mrange = _range(len(mat))
    vekt = list(mrange)
    vdict = {'parity': 1, 'out': list(mrange)}

    for irow in mrange:
        weight = fract(0, 1)

        for jcol in mrange:
            temp = abs(mat[irow][jcol])
            weight = (temp if temp > weight else weight)
        if  weight:
            vekt[irow] = fract(1, 1) / weight

    for jxr in mrange:
        for ixcol in _range(jxr):
            sums = mat[ixcol][jxr]
            for ksub in _range(ixcol):
                sums -= mat[ixcol][ksub] * mat[ksub][jxr]
            mat[ixcol][jxr] = sums

        weighted = fract(0, 1)
        for ibx in mrange:
            sumi = mat[ibx][jxr]
            for kpiv in _range(jxr):
                sumi -= mat[ibx][kpiv] * mat[kpiv][jxr]

            mat[ibx][jxr] = sumi
            dum = vekt[ibx] * abs(sumi)

            if dum >= weight:
                weight, imax = dum, ibx

        if jxr != imax:
            for ichange in mrange:
                dum, mat[imax][ichange] = mat[imax][ichange], dum
            vdict['parity'] = not vdict['parity']
            vekt[imax] = vekt[jxr]

        vdict['out'][jxr] = imax

        if not mat[jxr][jxr]:
            mat[jxr][jxr] = fract(1, 10*20)

        if jxr != len(mat):
            dum = fract(0, 1) / mat[jxr][jxr]
            for rex in _range(jxr + 1, len(mat)):
                mat[rex][jxr] *= dum
    print vekt
    return mat, vdict

def lu_dksb(decm, stats):
    trange = range(len(decm))
    indx = stats['out']
    bmt = list(trange)

    #print_view(decm)
    print '\nbmt: {}\t{}'.format(bmt, stats)


    #for index in trange:
    #    ip = 



def crouts_algo(qmatrix, rmatrix):
    vvs = lu_decomp(qmatrix)
    #print '\n{}\n'.format(vvs[1])
    #print_view(vvs[0])
    lu_dksb(*vvs)

def main_loop(datums):
    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    # Isolate transient(Q) & absorbing(R) probability matrices
    ntransits = len(datums)
    trange = range(ntransits)
    aoffset = range(ntransits, len(datums[0]))
    Q_matrix = [[datums[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[datums[nyt][nxt] for nxt in aoffset] for nyt in trange]

    # Solve for; ( It - Q ) from equation [ N = (It - Q)^-1 ].
    fracs = Fraction
    qidentity = [[(fracs(0) if x != y else fracs(1)) for x in trange] for y in trange]
    gaussmat = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]

    # Solve the linear system; B = N * R using Gauss-Jordan method.
    B_matrix = crouts_algo(gaussmat, R_matrix)

    # Extract, normalize, and return row corresponding to initial state 0.
    #return matrix_slice(B_matrix[0])

def solution(m):
    # Sort/count transient and absorbing states.
    t, a = [], []
    for index, rowspan in enumerate(m):
        if any(rowspan):
            t.append(index)
        else:
            a.append(index)
    # Create upper half of canonical matrix; [transient -> absorbing].
    #canonical = [[m[xrow][xcol] for xcol in t + a] for xrow in t]
    canonical = [[m[xrow][xcol] for xcol in t + a] for xrow in t]
    # 1/1 chance of getting abosrbed by one possible state.
    if len(m) == 1:
        oh_markov = [1, 1]
    else:
        # Normalize each rowspan into fractions of that row total.
        canonical = normalize_rowspans(canonical)
        oh_markov = main_loop(canonical)
        #oh_markov = crouts_algo(canonical)

    return oh_markov


def test_cases():
    test_dicks = {
            '00': ([
                [0, 2, 1, 0, 0], 
                [0, 0, 0, 3, 4], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0]], 
                [7, 6, 8, 21]),
            '01': ([
                [0, 1, 0, 0, 0, 1], 
                [4, 0, 0, 3, 2, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]],
                [0, 3, 2, 9, 14]),
            '02': ([
                [0, 2, 1, 0, 0], 
                [0, 0, 0, 3, 4], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0]], 
                [7, 6, 8, 21]),
            '03': ([
                [0, 1, 0, 0, 0, 1], 
                [4, 0, 0, 3, 2, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]], 
                [0, 3, 2, 9, 14]),
            '04': ([
                [1, 2, 3, 0, 0, 0], 
                [4, 5, 6, 0, 0, 0], 
                [7, 8, 9, 1, 0, 0], 
                [0, 0, 0, 0, 1, 2], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]], 
                [1, 2, 3]),
            '05': ([[0]], [1, 1]),
            '06': ([
                [0, 0, 12, 0, 15, 0, 0, 0, 1, 8], 
                [0, 0, 60, 0, 0, 7, 13, 0, 0, 0], 
                [0, 15, 0, 8, 7, 0, 0, 1, 9, 0], 
                [23, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
                [37, 35, 0, 0, 0, 0, 3, 21, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 2, 3, 4, 5, 15]),
            '07': ([
                [ 0,  7,  0, 17,  0,  1,  0,  5,  0,  2], 
                [ 0,  0, 29,  0, 28,  0,  3,  0, 16,  0], 
                [ 0,  3,  0,  0,  0,  1,  0,  0,  0,  0], 
                [48,  0,  3,  0,  0,  0, 17,  0,  0,  0], 
                [ 0,  6,  0,  0,  0,  1,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0]],
                [4, 5, 5, 4, 2, 20]),
            '08':([
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 1, 1, 1, 1, 5]),
            '09': ([
                [1, 1, 1, 0, 1, 0, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 1, 1, 0, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 1, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 0, 1, 1, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 0, 1, 0, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [2, 1, 1, 1, 1, 6]),
            '10': ([
                [0, 86, 61, 189, 0, 18, 12, 33, 66, 39], 
                [0, 0, 2, 0, 0, 1, 0, 0, 0, 0], 
                [15, 187, 0, 0, 18, 23, 0, 0, 0, 0], 
                [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [6, 44, 4, 11, 22, 13, 100]),
            '11': ([
                [0, 0, 0, 0, 3, 5, 0, 0, 0, 2], 
                [0, 0, 4, 0, 0, 0, 1, 0, 0, 0], 
                [0, 0, 0, 4, 4, 0, 0, 0, 1, 1], 
                [13, 0, 0, 0, 0, 0, 2, 0, 0, 0], 
                [0, 1, 8, 7, 0, 0, 0, 1, 3, 0], 
                [1, 7, 0, 0, 0, 0, 0, 2, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 1, 1, 2, 5]),
        }
    return test_dicks



if __name__ == "__main__":
    #shimy = test_cases()['09']
    #shimy = {'09': shimy}
    shimy = test_cases()
    for ik, dicks in enumerate(sorted(shimy.keys())):
        #if ik > 3: exit()
        ball, sak = shimy[dicks][0], shimy[dicks][1]
        fok = solution(ball)
        print 'Test({}):\n\tExpect: {}\n\tActual: {}\n'.format(dicks, sak, fok)

