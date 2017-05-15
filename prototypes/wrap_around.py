from collections import defaultdict

def main():
    weights = defaultdict(float)
    weights[(0,2,2)] = -1
    weights[(1,3,3)] = -1
    weights[(2,4,4)] = -1
    weights[(3,5,1)] = -1

    H = defaultdict(float)
    n = 5
    # fixed at four levels for testing
    for i in range(0,n+1):
        for j in range(0,n+1):
            for k in range(0,n+1):
                for l in range(0,n+1):
                    # print (i,j,k,l)
                    m = min( H[(i-1,j,k,l)], H[(i,j-1,k,l)], H[(i,j,k-1,l)], H[(i,j,k,l-1)] )
                    m = min(m, H[(i-1,j-1,k,l)] + weights[(0, i,j)])
                    m = min(m, H[(i,j-1,k-1,l)] + weights[(1, j,k)])
                    m = min(m, H[(i,j,k-1,l-1)] + weights[(2, k,l)])
                    # print  H[(i,j,k-1,l-1)] + weights[(2, k,l)]
                    m = min(m, H[(i-1,j,k,l-1)] + weights[(3, l,i)])    # note: weight is l,i rather than i,l

                    H[(i,j,k,l)] = m

    # print H[(1,1,0,0)]
    print H[(n,n,n,n)]



if __name__ == "__main__":
    main()
