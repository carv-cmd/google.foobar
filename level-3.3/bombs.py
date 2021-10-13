from argparse import ArgumentParser
from pprint import pprint
from random import seed, randint
from math import log, log1p

def args_parse():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-m')
    _parser.add_argument('-f')
    _parsed = _parser.parse_args()
    return int(_parsed.m), int(_parsed.f)

def fib(f):
    if f == 0:
        return 0
    elif f == 1:
        return 1
    else:
        try:
            return MEMO[f]
        except KeyError:
            newfib = fib(f-1) + fib(f-2)
            MEMO[f] = newfib
            return newfib


def notfib(m, f):
    mx, fy = int(m), int(f)
    target = m + f
    igen = -1
    while 1 not in [mx, fy]:
        if mx + fy > m + f:
            print '!!! Caught Exception  !!!'
            return 'impossible'
        else:
            imax, imin = ((mx, fy) if mx > fy else (fy, mx))
            mx, fy = imax % imin, imin
            igen += imax // imin
            print '|| mx: {} || fy: {} || igen: {} ||'.format(mx, fy, igen)
    print '\nMX: {}\nFY: {}'.format(mx, fy)
    return str(igen + max(mx, fy))

    
def solution01(m, f):
    col = '\033[38;5;201m'
    print '\n>>> Solving: {}( {} / {} )\033[00m'.format(col, m, f)
    mx, fy = int(m), int(f)
    kthgen = -1
    while 1 not in [mx, fy]:
        imax, imin = ((mx, fy) if mx > fy else (fy, mx))
        try:
            kthgen += imax // imin
            mx = imax % imin
            fy = imin
        except ZeroDivisionError:
            return 'impossible'
    dominant = max(mx, fy)
    print '\nDominant Clanker: {}'.format(dominant)
    return str(kthgen + dominant)

def printer(col, state, head, tail, igen, kgen):
    _color = '\033[38;5;{};48;5;{}m'.format(*col)
    _clear = '\033[00m'
    print '>> {}{}: [ head( {} ) : tail( {} ) : igen( {} ) : kthgen( {} ){}'.format(\
            _color, state, head, tail, igen, kgen, _clear)

def solution02(m, f):
    _divmod = divmod
    head, tail = int(m), int(f)
    kthgen = -1
    if head + tail == 3:
        return str(1)
    while head != 1:
        try:
            #printer((15, 202), 'init', head, tail, head//tail, kthgen)
            head, tail = ((head, tail) if head > tail else (tail, head))
            #printer((16, 202), 'Mids', head, tail, head//tail, kthgen)
            igen, head = _divmod(head, tail)
            kthgen += igen
            #printer((202, 0), 'POST', head, tail, igen, kthgen)
        except ZeroDivisionError:
            return 'impossible'
    print '\033[38;5;201m# Exit_State: head( {} ) / tail( {} )\033[00m'.format(
            head, tail)
    #return str(kthgen + (head if head > tail else tail))
    return str(kthgen + tail)


def solution03(m, f):
    _divmod = divmod
    head, tail = int(m), int(f)
    kthgen = -1
    
    # Orientation irrelevant when sum is 3, output always equals one.
    if head + tail == 3:
        kthgen, tail = 1, None
    else:
        # If head equals 1, found viable combination, break loop
        while head != 1:
            try:
                # Shuffle head and tail always assigning the largest to head
                # Divmod[0]: how many generations spawned
                
                printer((11, 16), 'init', head, tail, head//tail, kthgen)
                head, tail = ((head, tail) if head > tail else (tail, head))
                
                printer((201, 16), 'Mids', head, tail, head//tail, kthgen)
                igen, head = _divmod(head, tail)
                kthgen += igen

                printer((202, 16), 'POST', head, tail, igen, kthgen)
            
            except ZeroDivisionError:
                kthgen = 'impossible'
                break

    try:
        # Possible, kthgen + tail is the total generations spawned
        tgens = str(kthgen + tail)

    # Impossible and sum(3) exceptions.
    except TypeError:
        tgens = kthgen

    return tgens

def main(n=4, multiplier=4):
    func = solution03
    _col = '\033[38;5;39;48;5;16m'

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
            ]

    if n != 0:
        seed(512)
        for i in range(n):
            mrand, frand = randint(1, 10**multiplier), randint(1, 10**multiplier)
            tests.append(((mrand, frand), 'long'))
    
    test_set = []
    testpend = test_set.append
    for item, ex in tests:
        mxm, fxf = item[0], item[1]
        print '> Solving:{} ( {} / {} ) \033[00m'.format(_col, mxm, fxf)
        foobar = func(mxm, fxf)
        testpend(((str(mxm), str(fxf)), foobar))
        print '>>> Actual/Expected:{} [ {} / {} ] \033[00m\n'.format(\
                _col, foobar, ex)

    #pprint(test_set, width=60)
    #print(test_set)


xm, yf = args_parse()
main(xm, yf)
#main()





























