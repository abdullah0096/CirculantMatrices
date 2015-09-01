def PKCS_CM(pStart, pEnd, nStart, nEnd, dStart, dEnd):
    with open("PKCS_CM_ALL-" + str(pStart) + "-" + str(pEnd) + "-" + str(nStart) + "-" + str(nEnd) + "-" + str(dStart) + "-" + str(dEnd) + ".txt", 'w') as fp:
        with open("PKCS_CM-" + str(pStart) + "-" + str(pEnd) + "-" + str(nStart) + "-" + str(nEnd) + "-" + str(dStart) + "-" + str(dEnd) + ".txt", 'w') as fp1:
            fiveConditionCnt = 0
            iterationCnt = 0
            for p in range(pStart, pEnd):
                if is_prime(p):
                    for n in range(nStart, nEnd):
                        for d in range(dStart, dEnd):
                            if is_prime(d):
                                iterationCnt = iterationCnt + 1
                                print("\n" + str(p) + "  " + str(n) + "  " + str(d))
                                print("P is Prime")
                                print("d is Prime")
                                print("q is primitive mod d : TODO.......")
                                fp.write(str(p) + "  " + str(n) + "  " + str(d) + "\n")
                                fp.write("\nP is Prime")
                                fp.write("\nd is Prime")
                                fp.write("\nq is primitive mod d : TODO.......")
                                F = GF(p ** n, 'a')['x']
                                irrd = F.irreducible_element(d-1)

                                print("\n" + str(irrd))                        
                                fp.write("\n\n" + str(irrd))                        
                                                
                                mat = companion_matrix(irrd)
                                det_m = mat.determinant()
                                order_m = multiplicative_order(det_m)
                                phi_x = F.coerce(cyclotomic_polynomial(d))
                                from sage.symbolic.expression_conversions import polynomial

                                x_1 = polynomial(x-1, base_ring=GF(p ** n, 'a'))
                                x_d_1 = polynomial((x ** d)-1, base_ring=GF(p ** n, 'a'))

                                sci_x = crt(1, irrd, x_1, phi_x)
                        
                                sci_x = sci_x.mod(x_d_1)
                        
                                print("\n" + str(sci_x))
                                fp.write("\n\n" + str(sci_x))
                        
                                circ_mat = matrix(GF(p ** n, 'a'), d)
                                circ_mat[0] = vector(list(sci_x))

                                row0 = vector(list(sci_x))

                                col = d-1
                                for i in range(1, d):
                                    col2 = col
                                    for j in range(0, d):
                                        circ_mat[i, j] = row0[col2]
                                        col2 = (col2 + 1) % (d)
                                    col = (col -1) % (d)

                                A = matrix(GF(p ** n, 'a'), d)
                                A = circ_mat ** order_m
                                print("\n" + str(A))
                                fp.write("\n\n" + str(A))
                                row_sum = 0
                                for i in range(d):
                                    for j in range(d):
                                        row_sum = row_sum + A[i, j]

                                if row_sum != 1:
                                    print ("\nRow-Sum NOT 1... exiting")
                                    fp.write ("\n\nRow-Sum NOT 1... exiting")
                                    printStars(fp)
                                    continue
                                else:
                                    print("\nRow sum of  A :: %s" % row_sum)
                                    fp.write("\n\nRow sum of  A :: %s" % row_sum)

                                char_A = A.characteristic_polynomial()
                                char_A0 = char_A // x_1
                        
                                if char_A0.is_irreducible() != True:
                                    print ("\nIs XA/x-1 NOT irreduciable :: %s" % char_A0.is_irreducible())
                                    fp.write ("\nIs XA/x-1 NOT irreduciable :: %s" % char_A0.is_irreducible())
                                    printStars(fp)
                                    continue
                                else:
                                    print ("Is XA/x-1 irreduciable :: %s" % char_A0.is_irreducible())
                                    fp.write ("\n\nIs XA/x-1 irreduciable :: %s" % char_A0.is_irreducible())
                            
                                det_A = A.determinant();
                                if det_A != 1:
                                    print ("\nDeterminant of A is NOT 1. Det(A) :: " + str(det_A))
                                    fp.write ("\nDeterminant of A is NOT 1. Det(A) :: " + str(det_A))
                                    printStars(fp)
                                    continue
                                else:
                                    print ("\nDeterminant of A is :: " + str(det_A))
                                    fp.write ("\nDeterminant of A is :: " + str(det_A))                                    
                            
                                orderOf_A = A.multiplicative_order()
                                print ("\nOrder of A :: %s" % orderOf_A)
                                fp.write ("\nOrder of A :: %s" % orderOf_A)
                                fp1.write(str(p) + "  " + str(n) + "  " + str(d) + "\n")
                                fp1.write (str(orderOf_A)+"\n")
                                fac = ecm.factor(orderOf_A)

                                print ("Factors of order of A :: " + str(fac))
                                fp.write ("\n\nFactors of order of A :: " + str(fac))
                                size = len(fac)-1;
                                print("\nLargest Prime Factor :: " + str(fac[size]) + "\t log2 :: " + str(ceil(float(log(fac[size])))))
                                fp.write("\nLargest Prime Factor :: " + str(fac[size]) + "\t log2 :: " + str(ceil(float(log(fac[size])))))
                                fp1.write(str(fac[size]) + "\n" + str(ceil(float(log(fac[size])))))
                                printStars(fp)                                
                                printStars(fp1)
                                fiveConditionCnt = fiveConditionCnt + 1
            print ("\nTotal Number of Matrices satisfying all 5 Conditions :: " + str(fiveConditionCnt))
            print ("\nTotal Number of instances tested :: " + str(iterationCnt))
            fp.write ("\nTotal Number of Matrices satisfying all 5 Conditions :: " + str(fiveConditionCnt))
            fp.write ("\nTotal Number of instances tested :: " + str(iterationCnt) + "\n")            
            
def printStars(fp):
    fp.write("\n***************************************************************************************************\n")
    fp.write("***************************************************************************************************\n")
                      
                      