from argparse import ArgumentParser
from pprint import pprint
from random import seed, randint
from math import log, log1p

def args_parse():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-R')
    _parser.add_argument('-M')
    _parser.add_argument('--func')
    _parsed = _parser.parse_args()
    return int(_parsed.R), int(_parsed.M), int(_parsed.func)

def printer43(sync, pivot):
    _cul, _clr = '\033[38;5;202m', '\033[00m'
    _msg = '> {}iniSort: [ Head: {} | pivot: {} ]{}\n'
    print _msg.format(_cul, sync, pivot, _clr)

def printer44(mode, tmph, ngen, sync, pivot):
    col, clr = '\033[38;5;{};48;5;{}m', '\033[00m'
    c0, c1 = col.format(mode[0], 16), col.format(mode[1], 16)
    p01, p02 = (c0, tmph, pivot, clr), (c1, ngen, sync, clr)
    _msg1 = '{} lastHead( {} ), pivot( {} ) {}'.format(*p01)
    _msg2 = '{} igen=div( {} ) | sync=mod( {} ) {}'.format(*p02)        
    _out = '>> divmod({})\n=> [{}]'.format(_msg1, _msg2)
    print _out

def solution00(x, y):

    _divmod = divmod
    mx, fy = int(x), int(y)
    sync, pivot = ((fy, mx) if mx > fy else (mx, fy))
    printer43(sync, pivot)
    kthgen = -1

    while sync != 1:
        try:
            sync, pivot = pivot, sync
            tmpsync = sync
            igen, sync = _divmod(sync, pivot)
            printer44((201, 46, 1), tmpsync, igen, sync, pivot)
            kthgen += igen
        except ZeroDivisionError:
            kthgen = None
            break
    try:
        kthgen = str(kthgen + pivot)
    except TypeError:
        kthgen = 'impossible'
    return kthgen


def solution01(m, f):
    mx, fy = int(m), int(f)
    kthgen = -1
    while 1 not in [mx, fy]:
        imax = max(mx, fy)
        imin = min(mx, fy)
        try:
            kthgen += imax // imin
            mx = imax % imin
            fy = imin
        except ZeroDivisionError:
            return 'impossible'
    dominant = max(mx, fy)
    return str(kthgen + dominant)


# Subtracts least common divisor from x and y until sync/pivot irreducible.
# Pivot will always be the max(mx, fy) at iter start in while block.
# Swapping pivot and sync always assigns where; sync=max(x, y) and pivot=min(x, y).
# Previous min(x, y) is stored in pivot, sync becomes new min(x, y).
# Quotient is equal to number of generations possible w/ current sync units.
# Remainder is sync(new-min) to be swapped with the pivot(new-max).
def solution02(x, y):

    _divmod = divmod
    mx, fy = int(x), int(y)
    sync, pivot = ((mx, fy) if mx < fy else (fy, mx))
    kthgen = -1
    
    # Base case, sync/pivot irreducible therefore solution exists.
    while sync != 1:
        try:
            # Always pivots element were reducing from into sync.
            sync, pivot = pivot, sync
            igen, sync = _divmod(sync, pivot)
            kthgen += igen

        # ZeroDivisionError implies solution does not exist.
        except ZeroDivisionError:
            kthgen = None
            break

    # Raised ZeroDivisionError, solution doesnt exist.
    # Else, combine kth-generations w/ final remainder.
    if kthgen is None:
        kthgen = 'impossible'
    else:
        kthgen = str(kthgen + pivot)
    return kthgen


def main(n=4, multiplier=4, _func=0):
    if not _func:
        func = solution00
    elif _func == 1:
        func = solution01
    elif _func == 2:
        func = solution02

    tests = [
            (('2', '4'), 'impossible'), 
            (('376', '40'), 'impossible'), 
            (('4', '7'), '4'), 
            (('2', '1'), '1'), 
            (('3', '2'), '2'), 
            (('5', '7'), '4'), 
            (('2', '5'), '3'), 
            (('7', '9'), '5'), 
            (('148', '197'), '52'), 
            (('24387', '2242'), '33'), 
            (('489', '433'), '17'), 
            (('772', '243'), '18'), 
            (('512', '289'), '14'), 
            (('319', '323'), '83'),
            (('53543', '916308'), '54'),
            #(('54368795463', '94321679320'), '131'),
            ]

    if n != 0:
        seed(512)
        for i in range(n):
            mrand, frand = randint(1, 10**multiplier), \
                    randint(1, 10**multiplier)
            tests.append(((mrand, frand), 'rand'))
    
    _col1 = '\033[38;5;227;48;5;16m'
    _col2 = '\033[38;5;51;48;5;16m'
    _clr = '\033[00m'
    test_set = []
    testpend = test_set.append
    for item, ex in tests:
        mxm, fxf = item[0], item[1]
        #_msg = '\n> {}Head: {}{}\n> {}Pivot: {}{}'
        #print _msg.format(_col1, mxm, _clr, _col1, fxf, _clr)
        #foobar = func(mxm, fxf)
        #_msg2 = '\n>>> Actual/Expected:{} [ {} / {} ]{}\n'
        #print _msg2.format(_col2, foobar, ex, _clr)
        #testpend(((str(mxm), str(fxf)), foobar))
        func(mxm, fxf)

    #pprint(test_set, width=60)
    print '\n->>> Jeff Zuckermusk ->>>'


def coloop():
    for i in range(0, 256):
        _clr = '\033[00m'
        _cols = '>>> \033[38;5;{};48;5;16m'.format(i)
        print '{}* Color({}) *{}'.format(_cols, i, _clr)


if __name__ == "__main__":
    _range, _multiplier, _func = args_parse()
    main(_range, _multiplier, _func)
    #coloop()








