###########################################################################################
# Compute the Sylvester Denumerant through a q-Partial Fraction of
# the generating function F(x)=p(x)/((1-x)^n0(1-x^n1)^r1...(1-x^nk)^rk)
# INPUT:
#    t: the value for which we want to compute the denumerant
#    p: numerator, a polynomial in the ring Q[x]
#    N: exponents corresponding to the denominator
#        N=[(n1,r1),...,(nk,rk)]
#    m1: m=m1+r1+...rk, where m is the power of (1-x) in the denominator (default value: 0)
#
# OUTPUT: returns a value or prints based on the set parameters mentioned below
#   part: determines the output, (default value: 'full')
#         'full' returns the sylvester denumerant d(t,A)
#         'polynomial' returns the polynomial part W_1(t,A)
#         'periodic' returns the periodic part W_n1(t,A)+\cdots+Wnk(t,A)
#         'periodic list' returns each periodic part in list form [0,W_n1(t,A),...,W_nk(t,A)]
#         'generalized fds' returns the top-oder coefficients of the periodic part (assumes r1=...=rk)
#   print_qPF: if set to True, prints the q-Partial Fraction of F(x)  (default value: False)
#
# EXAMPLE:
#   #Computing d(t,(3,11,11)) for all values in [0,200]
#   P=PolynomialRing(QQ,'x')
#   p=P(1)
#   m1=0
#   N=[(3,1),(11,2)]
#   L1,L2,L3=[],[],[]
#   for t in range(0,201):
#       L1.append(DENUMERANT(t,p,N,m1,part='full'))
#       L2.append(DENUMERANT(t,p,N,m1,part='polynomial'))
#       L3.append(DENUMERANT(t,p,N,m1,part='periodic'))
#
#   #Compute d(73,(9,17,31))
#   DENUMERANT(73,p,[(9,1),(17,1),(31,1)],0,part='full',print_qPF=True)
#
#   #Generalized FDS terms for (3,3,3,4,4,4)
#   P=PolynomialRing(QQ,'x')
#   L=[]
#   for t in range(1,60):
#       L.append(DENUMERANT(t,P(1),[(3,3),(4,3)],0,part='generalized FDS'))
#
###########################################################################################

def DENUMERANT(t,p,N,m1=0,part='full',print_qPF=False):

    # check for N[0] to be pairwise relatively prime and they are not equal to 1
    for i in range(len(N)):
        if N[i][0] !=1:
            for j in range(i+1,len(N)):
                 if gcd(N[i][0],N[j][0]) !=1:
                        print("Error: The numbers must be pairwise relatively prime.")
                        return
        else:
            print("Error: Cannot inlcude 1 in the list.")
            return

    if part == 'generalized FDS':
        for i in range(len(N)):
            if N[i][1]!=N[0][1]:
                print("Error: All multiplicities must be same.")
                return 
            
    m=sum([z[1] for z in N])+m1

    poly=[]
    poly.append(P((1-x)^(m)))
    for i in range(len(N)):
        poly.append(P(sum([x^j for j in range(0,N[i][0])])^(N[i][1])))

    k=len(poly)
    f=[]
    for i in range(0,k):
        s=1
        for j in range(0,k):
            if i!=j:
                s=s*poly[j]
        f.append(eval(p,s,poly[i]))


    S=str("(")+str(f[0])+str(")")+"/"+str("(")+str(poly[0])+str(")\n\n")
    for j in range(1,k):
        f[j]=f[j]*P((1-x)^(N[j-1][1]))
        poly[j]=poly[j]*P((1-x)^(N[j-1][1]))
        S=S+str("  +")+str("(")+str(f[j])+str(")")+"/"+str("(")+str(poly[j])+str(")\n\n")
    if print_qPF:
       print("The q-Partial Fractions:\n"+S)

    # Expressing Polynomial Part in terms of (1-x)^j
    g=f[0]
    poly_part=[]
    while g.degree()!=0:
        d=g.quo_rem(P(1-x))
        poly_part.append(d[1])
        g=d[0]
    poly_part.append(g(0))

    S_poly=str(poly_part[0])
    for j in range(1,len(poly_part)):
        S_poly=S_poly+"+"+str(poly_part[j])+"(1-x)^"+str(j)
    if print_qPF:
        temp_denominator=P((1-x)^m)
        print("Simplified form of the 'polynomial part': \n"+"("+S_poly+")"+"/"+"("+str(temp_denominator)+")\n\n")



    # compute the denumerant
    poly_sum=0
    for j in range(0,m):
        poly_sum=poly_sum+poly_part[j]*binomial(t+m-1-j,t)

    periodic_sum=0
    periodic_sum_list=[0]
    gen_fds=0
    for j in range(1,k):
        temp_sum=0
        L=gen_taylor(f[j],N[j-1][0])
        for i in range(len(L)):
            temp_sum=temp_sum+L[i][t%N[j-1][0]]*binomial(floor(t/N[j-1][0])+N[j-1][1]-1-i,floor(t/N[j-1][0]))
            if i==0: #assumes that r1=...=rk
                gen_fds = gen_fds + L[0][t%N[j-1][0]]
        periodic_sum_list.append(temp_sum)
    periodic_sum=sum(periodic_sum_list)

    if part=='periodic':
        return periodic_sum
    elif part == 'polynomial':
        return poly_sum
    elif part == 'periodic list':
        return periodic_sum_list  #returns a list of contributions of each value
    elif part == 'generalized FDS': #assumes that r1=...=rk
        return gen_fds
    else:
        return (poly_sum+periodic_sum)


###########################################################################################
# Extended Cover-Up Method
# INPUT:
#        r: polynomial in the numerator
#        s: list of pairwise relatively prime polnomial factors in the denominator
# OUTPUT:
#        prints the partial fractions for r/(s[0]*...*s[k-1]), where k=len(s)
# EXAMPLE:
#   P=PolynomialRing(QQ,'x')
#   extended_cover_up(P(1),[P(x^2+1), P(x-1)])
###########################################################################################
def extended_cover_up(r,s):
    #number of factors in the denominator
    k=len(s)

    for i in range(0,k):
        for j in range(i+1,k):
            if gcd(s[i],s[j])!=1:
                print("Invalid input: the factors should be pairwise relatively prime.")
                return

    #computing the numerators in the partial fraction using the extended cover-up method
    f=[]
    for i in range(0,k):
        t=1
        for j in range(0,k):
            if i!=j:
                t=t*s[j]
        f.append(eval(r,t,s[i]))

    #converting to string
    S=str("(")+str(f[0])+str(")")+"/"+str("(")+str(s[0])+str(")\n")
    for j in range(1,k):
         S=S+str("+")+str("(")+str(f[j])+str(")")+"/"+str("(")+str(s[j])+str(")\n")

    print("The partial fraction is:\n"+S)


# eval function with assumption gcd(s,a)=1
def eval(r,s,a):
    h=xgcd(s,a)
    return (r*h[1])%a

# a generalized taylor series with terms (1-x^b)^j
def gen_taylor(f,b):
    P=PolynomialRing(QQ,x)
    L=[]
    for j in range(0,floor(f.degree()/b)+1):
        L.append(((-1)^(j)/factorial(j))*eval(f,P(1),P(1-x^b)))
        f=Diff_b(P(f),b)
    return L

def Diff_b(f,b):
    P=PolynomialRing(QQ,x)
    ff=list(f)
    g=P(0)
    for j in range(0,len(ff)):
        g=g+P(ff[j])*diff_b(j,b)
    return g

def diff_b(k,b):
    P=PolynomialRing(QQ,x)
    if k>=b:
        return floor(k/b)*P(x^(k-b))
    else:
        return P(0)











