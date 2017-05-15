#!/usr/bin/python

import random
'''
def getNeighbors(elements):
    barrier_index = elements.index(-1)
    
    for i in range(barrier_index): # i before barrier
        for j in range(barrier_index+1,len(elements)): # j after barrier
            new_list = elements[:i] + elements[i:j+1][::-1] + elements[j+1:]
            # print new_list
            yield new_list
'''
def getNeighbors(elements):
    barrier_index = elements.index(-1)
    
    for i in range(barrier_index+1): # i before barrier
        for j in range(barrier_index,len(elements)): # j after barrier
            if i == j: continue

            new_list = elements[:i] + elements[i:j+1][::-1] + elements[j+1:]
            # print new_list
            if new_list[0] == -1 or new_list[-1] == -1:
                continue
            yield new_list


# ensure that there is at least 1 even and atleast 1 odd
def main():
    N = 4

    
    elements = list(range(N))

    # randomly intialize to some labeling by inserting a barrier to random list
    random.shuffle(elements)
    rand_index = random.randint(1, len(elements)-1)
    elements = elements[:rand_index] + [-1] + elements[rand_index:] 
        
    elements = [0,1,-1,2,3]

    print elements
    print ""
    
    for neighbor in getNeighbors(elements):
        barrier_index = neighbor.index(-1)
        even = neighbor[:barrier_index]
        odd = neighbor[barrier_index+1:]
        print even, odd







if __name__ == '__main__':
    main()