from argparse import ArgumentParser
from random import seed, sample
from time import sleep
from collections import defaultdict

def args_parse():
    _parser = ArgumentParser()
    _parser.add_argument('-m')
    _parser.add_argument('-f')
    _parsed = _parser.parse_args()
    return int(_parsed.m), int(_parsed.f)


def rand_sample(ir=1000000,  ik=10000, _seeder=1024):
    seed(_seeder)
    return sample(range(ir), ik)


def bsearch(array, elem):
    lo = 0
    mid = len(array) // 2
    hi = len(array) - 1
    while True:
        if lo >= mid >= hi:
            if array[mid] != elem:
                return None
            else:
                return array[mid]
        if array[mid] == elem:
            return array[mid]
        elif array[mid] < elem:
            lo = mid + 1
        else:
            hi = mid - 1
        mid = (hi + lo) // 2


def bfinder(element):
    _test_set = rand_sample()
    print _test_set[:50]
    print 'Sorting'
    _test_set = sorted(_test_set)
    _fin = 4265
    _foo = bsearch(_test_set, _fin)
    print '\nsearch({}): -found-> {}'.format(_fin, _foo)


def newgen(seedm, seedf):
    # parbombF = seedm + seedf
    # parbombM = seedf
    return parbombF, parbombM


def constructor(pNodeM, pNodeF):
    #############################################
    # Gen_Classes:
    #############################################
    # M_cycle:
    #   cNodeM = pNodeM + pNodeF
    #   cNodeF = pNodeF 
    #   mdominant_gen = cNodeM + cNodeF
    #   return mdominant_gen
    #   *mdominant_gen = ( ( pNodeM + pNodeF ) + ( pNodeF ) )
    #############################################
    # F_cycle:
    #   childM = pNodeM
    #   childF = pNodeM + pNodeF
    #   fdominant_gen = cNodeM + cNodeF
    #   return fdominant_gen
    #   *fdominant_gen = ( ( pNodeM ) + ( pNodeM + pNodeF ) )
    #############################################
    # F-cycle && M-cycle : symetric
    # gen_total = SIGMA[mdominant, fdominant] 
    # gen_total = (
    #   ((pNodeM_i + pNodeF_i) + (pNodeF_i)) +
    #   ((pNodeM_i) + (pNodeM_i + pNodeF_i))
    # )
    #############################################
    pass

def replicator(m, f):
    modM, modF = 1, 1
    count = 0
    print '>>> Iter({}): [ M: {} / F: {} ]'.format(
            count, modM, modF)
    while True:
        if (m == 0) or (f == 0):
                return 1
        else:
            replicator(m-1, f) 
            replicator(m, f-1)
        sleep(1)
    print '\nFinal-M: {}\nFinal-F: {}'.format(m, f)


def fib(x):
    #print '>>> Fibbing: {}'.format(x)
    if x in [0, 1]: 
        return 1
    else:
        return fib(x-1) + fib(x-2)

def solution(x, y):
    bmin = min(x, y)
    target = x + y
    xt, yt = 1, 1
    count = 0
    while True:
        print 'iter({}): [ xt: {} ] [ yt: {} ] = {}'.format(
                count, xt, yt, xt + yt)
        count += 1
        if xt + yt == target:
            break
        elif xt + yt > target:
            count = 'impossible'
            break
        elif xt + yt < target:
            if xt + yt == bmin:
                xt = bmin
            yt += xt
            target -= xt
        sleep(1)
    return count

def memofib(x, y):

    def fib(mf):
        if not mf:
            return 'impossible'

        if mf in (x, y):
            print '# f({}) in [x({}), y({})] #'.format(mf, x, y)
            return 1

        try:
            return memos[mf]

        except KeyError:
            print '>>> iter_fib >>>'
            memos['gens'] += 1
            sleep(0.5)
            tmp = fib(mf-1) + fib(mf-2)
            memos[mf] = tmp
            return tmp

    memos = {'gens': 0}
    fib(x + y)
    return memos['gens']

def main():
    # Input: -> ('4', '7') ==> 4 
    # Input: -> ('2', '1') ==> 1
    # Input: -> ('3', '2') ==> 2
    # Input: -> ('2', '4') ==> None
    # Input: -> ('2', '3') ==> ? 3 ?
    # Input: -> ('5', '7') ==> ? None ?
    repl = replicator(3, 2)
    print '\nReturn: {}'.format(repl)


marg, farg = args_parse()
#print '\n>>> Naive Fib {}'.format(fib(_initargs))
solute = memofib(marg, farg)
#solute = solution(marg, farg)
print "\n>>> Fibonacci'd: {}".format(solute)















