from collections import defaultdict


test = ([2, -3, 1, 0, -5], [2, 0 ,2, 2, 0], [-2, -3, 4, -5], [0, 0, 0], [-2])

dick = defaultdict(int)

for b in test:
    for elem in b:
        if elem not in dick:
            dick[elem] = elem
        else:
            dick[elem] *= elem

    dkeys = dick.keys()

    try:
        assert any(dkeys), 'no deek'
        print any(dkeys)

    except AssertionError as e:
        print e
    
    dick.clear()

