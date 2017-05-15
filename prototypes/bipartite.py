#!/usr/bin/python

from collections import defaultdict


def recursive(indices, weights, numEven, numOdd, memo):
    if indices in memo:
        return memo[indices]

    # print indices

    mx = 0
    for p, i in enumerate(indices):
        if i > 0:
            newIndices = list(indices)
            newIndices[p] -= 1
            newIndices = tuple(newIndices)

            val = recursive(newIndices, weights, numEven, numOdd, memo)
            mx = max(val, mx)


    # now for pairwise
    for pair_id in range(numEven*numOdd):
        o = pair_id % numOdd
        e = pair_id / numOdd

        e_index = e
        o_index = numEven + o

        curr_even_pos = indices[e_index]
        curr_odd_pos =  indices[o_index]

        if curr_even_pos > 0 and curr_odd_pos > 0:
            newIndices = list(indices)
            newIndices[e_index] -= 1
            newIndices[o_index] -= 1
            newIndices = tuple(newIndices)

            if True or ((e,o), curr_even_pos, curr_odd_pos) in weights:

                # print 'looking at edge =', ((e,o), curr_even_pos, curr_odd_pos)

                val = recursive(newIndices, weights, numEven, numOdd, memo) +\
                     weights[((e,o), curr_even_pos, curr_odd_pos)]
                mx = max(val, mx)


    memo[indices] = mx
    return mx


def main():

    numEven = 2
    numOdd = 2

    totalLevels = numEven + numOdd
    
    '''
    |even| = p, |odd| = q
    pairs are always (even, odd) are number as follows:
    0: (e_0, o_0)
    1: (e_0, e_1)
    ...
    q-1: (e_0, o_{q-1})
    q:   (e_1, o_0)
    q+1: (e_1, o_1), 
    ...
    2q:  (e_2, o_0)
    ...

    pair_id(even_id, odd_id) = even_id*q + odd_id
    '''


    num = 0
    p, q = numEven, numOdd
    for even_id in range(numEven):
        for odd_id in range(numOdd):
            print num, (even_id, odd_id), 
            id = even_id*q + odd_id
            print id,

            o = id % q
            e = id / q

            print (e,o)
            num += 1





    weights = defaultdict(float)
    '''
    weights[((0,0), 1, 1)] = 1
    # weights[((0,0), 1, 3)] = 10
    weights[((0,0), 3, 3)] = 1
    weights[((1,0), 2, 2)] = 1
    weights[((2,0), 4, 1)] = 7
    '''

    weights[((0,0), 4, 4)] = 1
    weights[((1,0), 5, 5)] = 1
    weights[((1,1), 6, 6)] = 1 
    weights[((0,1), 7, 7)] = 1
    weights[((0,1), 1, 9)] = 10
    

    # matrix should be H([pos(e_0), pos(e_1), ... , pos(e_{q-1}), pos(o_0), pos(e_1), ..., pos(e_{p-1})])

    H = defaultdict(float)
    n = 10


    indices = tuple([n]*numEven + [n]*numOdd)
    v = recursive(indices, weights, numEven, numOdd, {})
    print v




if __name__ == "__main__":
    main()
