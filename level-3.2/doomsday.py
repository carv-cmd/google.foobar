from pprint import pprint
from fractions import Fraction, gcd


def print_view(views):
    foo = [['{}/{}'.format(i.numerator, i.denominator) for i in ixi] for ixi in views]
    print; pprint(foo, width=len(views[0]) * 20 + 10)


def gauss_redux(q_mx, r_mx):
    # Implemented from "Numerical Recipes in C" Gauss-Jordan elimination example.
    # Modifes matrices in place, returning the funadmental matrix for slicing.
    # Back substitution eliminated through prior canonical sorting.

    _range = xrange
    _frac = Fraction
    trange = _range(len(q_mx))
    addone, zed = _frac(1, 1), _frac(0, 1)

    # Keep records of pivot operations
    # Since inverse is redundant, isolate just the B matrix we want.
    audit = {
            'r': list(trange), 
            'c': list(trange),
            'p': [_frac(0) for _ in trange]
            }

    # Main iterator to perform reductions
    for mains in trange:
        weight = zed
        # Find r/c dominant pivot
        for jx in trange:
            if audit['p'][jx] != 1:
                for ky in trange:
                    sabs = abs(q_mx[jx][ky])
                    if (audit['p'][ky] == 0) and (sabs >= weight):
                        weight, irow, icol = sabs, jx, ky

        # Record new pivot plan
        audit['p'][icol] += addone
        audit['r'][mains] = irow
        audit['c'][mains] = icol

        # Perform columnwise reductions on both Q and R matrices.
        pinvt, q_mx[icol][icol] = addone/q_mx[icol][icol], addone
        for mx in (q_mx, r_mx):
            for cx in xrange(len(mx[0])):
                mx[icol][cx] *= pinvt

        # Perform rowise reductions on both Q and R matrices, skipping diagonal.
        for rox in trange:
            if rox != icol: 
                redux, q_mx[rox][icol] = q_mx[rox][icol], zed
                for mx in (q_mx, r_mx):
                    for cx in xrange(len(mx[0])):
                        mx[rox][cx] -= mx[icol][cx] * redux
    return r_mx


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


def main_loop(datums):
    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    ntransits = len(datums)
    trange = range(ntransits)
    aoffset = range(ntransits, len(datums[0]))

    # Isolate transient(Q) & absorbing(R) probability matrices.
    Q_matrix = [[datums[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[datums[nyt][nxt] for nxt in aoffset] for nyt in trange]
    
    # Solve for; ( It - Q ) from equation ( N = (It - Q)^-1 ).
    fracs = Fraction
    qidentity = [[(fracs(0) if x != y else fracs(1)) for x in trange] for y in trange]
    iQ_matrix = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]
    
    # Solves linear system; B = N^-1 * R.
    B_matrix = gauss_redux(iQ_matrix, R_matrix)
    
    # Extract, normalize, and return row corresponding to initial state 0.
    return matrix_slice(B_matrix[0])


def solution(m):
    # 1/1 chance of getting abosrbed by one possible state.
    if len(m) == 1:
        oh_markov = [1, 1]
    else:
        t, a = [], []
        # Sort/count transient and absorbing states.
        for index, rowspan in enumerate(m):
            if any(rowspan):
                t.append(index)
            else:
                a.append(index)

        # Create upper half of canonical matrix; [transient -> absorbing].
        canonical = [[m[xrow][xcol] for xcol in t + a] for xrow in t]

        # Normalize each rowspan into fractions of that row total.
        canonical = normalize_rowspans(canonical)

        # Solve B = N * R linear system
        oh_markov = main_loop(canonical)

    return oh_markov


def test_cases():
    return {
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
            '05': ([
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
            '06': ([
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
            '07':([
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
            '08': ([
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
            '09': ([
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
            '10': ([
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


def main():
    shimy = test_cases()
    for vit in sorted(shimy.keys()):
        vitlen, expected, actual = (len(shimy[vit][0]), shimy[vit][1], solution(shimy[vit][0]))
        print """\r->>> Test( \033[38;5;196m{}\033[00m |  {}x{} ) <<<-
        \r~> Expecting: \033[38;5;87m{}\033[00m
        \r~~> Returned: \033[38;5;201m{}\033[00m\n""".format(
                vit, vitlen, vitlen, expected, actual)

if __name__ == "__main__":
    main()
    
