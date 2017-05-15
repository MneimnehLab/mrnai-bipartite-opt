#!/usr/bin/python

import sys

def main():
    inp = sys.stdin.readlines()
    sequences = inp[0].strip().split("&")
    names = inp[1].strip().split("&")

    # convert "a,b,c|d,e,f" to even = [a,b,c], odd = [d,e,f]
    evens, odds = map(lambda x: map(int, x.split(",")), inp[2].strip().split("|"))

    # print evens
    # print odds

    print '%d\t%d' % (len(evens), len(odds))
    
    count = 0
    for e in evens:
        for o in odds:
            fname = "output/default-%d_%d_itemized.out" % (e, o)
            with open(fname) as f:
                for line in f:
                    print '%d\t%s' % (count, line.strip())

            count += 1



if __name__ == '__main__':
    main()
