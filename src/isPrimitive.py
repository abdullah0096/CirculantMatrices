def isPrimitive(n, p):
    isPrimitiveRoot = 1
    N = n % p
    n1 = numpy.array(range(1, p + 1), dtype=int)
    for i in range(1, p):
        n1[i] = (N ** i) % p
    
    for i in range(1, p):
        isPresent = 0
        for j in range(1, p):
            if n1[j] == i:
                isPresent = 1
                break
        if isPresent != 1:
            isPrimitiveRoot = 0
            break
    
    if isPrimitiveRoot == 1:
        return 1
    else:
        return 0
    
def check_isPrimitive(n):
    for j in range (40, 500):
        print ("\n p :: " + str(j))
        for i in range (11, 200):           
            if is_prime(i):
                ans = isPrimitive((2 ** j), i)
                if ans == 1:
                    print (str(i) + ",\t")
        print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")