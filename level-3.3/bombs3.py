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

