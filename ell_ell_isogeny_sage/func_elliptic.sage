

#construct big finite field-------


#return finite field F_{p^4}.
def GFp4pow(p:int):
    Fp2=GF(p**2)
    z=Fp2(0)
    while z.is_square():
        z=Fp2.random_element()
    t=ZZ(z + z**p)
    n=ZZ(z**(p+1))
    R = PolynomialRing(Fp2,'x')
    x=R.gen()
    f=x**4-t*x**2+n
    Fp4=GF(p**4,modulus=f,name='a')
    #Fp4=Fp2.extension(f)
    assert(order(Fp4)==p**4)
    return Fp4,x,Fp2



def In_subfield(a):
    K=parent(a)
    p=K.characteristic() 
    assert(order(K)==p**4)
    assert(a in K)
    return a**(p**2)==a




def In_subfield_tc(tc:Coord):
    t1=tc.numer[1]/tc.numer[0]
    t2=tc.numer[2]/tc.numer[0]
    t3=tc.numer[3]/tc.numer[0]
    return In_subfield(t1) and In_subfield(t2) and In_subfield(t3)


#Cyclic isogeny--------------------


def Decomp_degree(N:int):
    fac=list(factor(N))
    compo_list=[]
    for i in range(0,len(fac)):
        for j in range(0,fac[i][1]):
            compo_list.append(fac[i][0])
    return compo_list





def Elliptic_Cyclic(E,ker,P,Q):
    deg=order(ker)
    fac_list=Decomp_degree(deg)
    #print("Elliptic isogeny",fac_list)
    s=1
    for i in range(0,len(fac_list)):
        l=fac_list[i]
        #print(l)
        k=deg//(s*l)
        decomp_ker=k*ker
        isogeny=EllipticCurveIsogeny(E,decomp_ker)
        E=isogeny.codomain()
        ker=isogeny(ker)
        P=isogeny(P)
        Q=isogeny(Q)
        s*=l
    return E,P,Q










#theta function on elliptic cueve--------------------------



def Elliptic_to_Legendre(E):
    coeff=E.a_invariants()
    #y^2=x^3+ax^2+bx+c
    a=coeff[2]
    b=coeff[3]
    c=coeff[4]
    fld=E.base_field()
    X=gen(fld['X'])
    f=X**3+a*X**2+b*X+c
    roots=f.roots()
    assert(len(roots)==3)
    x1=roots[0][0]
    x2=roots[1][0]
    x3=roots[2][0]
    lmd=(x3-x1)/(x2-x1)
    E_lmd=EllipticCurve(fld,[0,(-lmd-1),0,lmd,0])
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd



'''
def Elliptic_to_Legendre_with_basis(E,N:int):
    fld=E.base_field()
    #assert(K==base_field)
    p=fld.characteristic()
    assert((p+1)%N==0)
    lmd,E_lmd,_=Elliptic_to_Legendre(E)
    S1,S2=E_lmd.torsion_basis(N)
    lmd=fld(lmd)
    E_lmd=EllipticCurve(fld,[0,(-lmd-1),0,lmd,0])
    S1=E_lmd(S1)
    S2=E_lmd(S2)
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd,S1,S2
'''


def Elliptic_to_Legendre_with_basis(E,N):
    coff=E.a_invariants()
    #y^2=x^3+ax^2+bx+c
    K=E.base_ring()
    Fp2=K.subfield(2)
    a=Fp2(coff[2])
    b=Fp2(coff[3])
    c=Fp2(coff[4])
    Y=gen(Fp2['Y'])
    f=Y**3+a*Y**2+b*Y+c
    roots=f.roots()
    x1=roots[0][0]
    x2=roots[1][0]
    x3=roots[2][0]
    lmd=(x3-x1)/(x2-x1)
    E_lmd=EllipticCurve(Fp2,[0,(-lmd-1),0,lmd,0])
    S1,S2=E_lmd.torsion_basis(N)
    E_lmd=EllipticCurve(K,[0,(-lmd-1),0,lmd,0])
    S1=E_lmd(S1)
    S2=E_lmd(S2)
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd,S1,S2





def Legendre_to_Elliptic(lm):
    E=EllipticCurve(parent(lm),[0,(-lm-1),0,lm,0])
    return E




def Legendre_to_lv2tnp(lm):
    K=parent(lm)
    sq_rt_lm  =sqrt(lm)
    sq_rt_lmm1=sqrt(lm-1)
    assert(sq_rt_lm**2==lm)
    assert(sq_rt_lmm1**2==lm-1)
    thnp0_sq=sq_rt_lm
    thnp1_sq=sq_rt_lmm1
    thnp2_sq=K(1)
    lv2tnp=[(thnp0_sq+thnp2_sq),thnp1_sq]
    return lv2tnp,sq_rt_lm,sq_rt_lmm1




def Lv2tnp_to_Legendre(lv2tnp:list):
    assert(len(lv2tnp)==2)
    a=lv2tnp[0]
    b=lv2tnp[1]
    sq_rt_lm=(a**2+b**2)/(a**2-b**2)
    sq_rt_lmm1=(1+sq_rt_lm)*(lv2tnp[1]/lv2tnp[0])
    lm=sq_rt_lm**2
    assert(sq_rt_lmm1**2+1==lm)
    return lm,sq_rt_lm,sq_rt_lmm1
    


def Is_isomorphic_Legendre(lmd_1,lmd_2):
    if lmd_1 in {lmd_2,1/lmd_2,1-lmd_2,1/(1-lmd_2),1-(1/lmd_2),lmd_2/(lmd_2-1)}:
        return True
    else:
        return False





#-----------------------




#Legendre representaion to lv2 theta coordinate.
def Leg_lv2(lm,sq_rt_lm,sq_rt_lmm1,lv2tnp:list,pt):
    assert(sq_rt_lm**2==lm)
    assert(sq_rt_lmm1**2==lm-1)
    u=pt[0]
    v=pt[1]
    w=pt[2]
    if w==0:
        return lv2tnp
    assert(w==1)
    assert((v**2)==u*(u-1)*(u-lm))
    thc0_sq=sq_rt_lm*(u-1)
    thc1_sq=sq_rt_lmm1*u
    thc2_sq=sq_rt_lm*thc0_sq-sq_rt_lmm1*thc1_sq
    thc3_sq=sq_rt_lm*thc1_sq-sq_rt_lmm1*thc0_sq
    #lv2=[(parent(u))(1),(thc1_sq+thc3_sq)/(thc0_sq+thc2_sq)]
    lv2=[(thc0_sq+thc2_sq),(thc1_sq+thc3_sq)]
    return lv2





def Lv2_to_Legendre(lm,sq_rt_lm,sq_rt_lmm1,lv2tc:list):
    r=lv2tc[1]/lv2tc[0]
    u=(sq_rt_lm*sq_rt_lmm1+r*(lm+sq_rt_lm))/(r*(1+sq_rt_lm)-sq_rt_lmm1)
    return u
    


#--------------------------------------


def Product_theta(lv2_1:list,lv2_2:list):
    assert(len(lv2_1)==2)
    assert(len(lv2_2)==2)
    lv2=[lv2_1[0]*lv2_2[0],
         lv2_1[0]*lv2_2[1],
         lv2_1[1]*lv2_2[0],
         lv2_1[1]*lv2_2[1]]
    return lv2



def Theta_Hadamard(lv2_x:list):
    assert(len(lv2_x)==4)
    x=lv2_x[0]
    y=lv2_x[1]
    z=lv2_x[2]
    w=lv2_x[3]
    dula_lv2_x=[x+y+z+w,
                x-y+z-w,
                x+y-z-w,
                x-y-z+w]
    return dula_lv2_x

#------------------------------------





#----------------------------------------------------------------



#cf.[LR16]section5
#sum_{t}(chi(t)lv2_x{i+t}lv2_x(t))
def Product_term(lv2_x,i,chi):
    K=parent(lv2_x[0])
    sum=K(0)
    for t in range(0,4):
        chi_t=(-1)**((chi//2)*(t//2)+(chi%2)*(t%2))
        ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
        sum+=(chi_t*lv2_x[ipt]*lv2_x[t])
    return sum


#from lv2 theta null point to lv(2,2) theta null squard.
def Lv2tnp_to_lv22tnpsq(lv2tnp):
    lv22tnpsq={}
    for chi in range(0,4):
        for i in range(0,4):
            lv22tnpsq[(chi,i)]=Product_term(lv2tnp,i,chi)
    return lv22tnpsq
            


even_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==0}

odd_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==1}



#check if product of elliptic curves without theta structure.
def Is_Ellipitc_product(lv2tnp):
    lv22tnpsq=Lv2tnp_to_lv22tnpsq(lv2tnp)
    num_zero=len({k for k in lv22tnpsq if lv22tnpsq[k]==0})
    assert(num_zero==6 or num_zero==7)
    for i_chi in odd_theta_lv22: #odd theta=0
        i=i_chi[0]
        chi=i_chi[1]
        assert(lv22tnpsq[(i,chi)]==0)
    if lv22tnpsq[(3,3)]==0:
        return True,True #elliptic product with product theta struture.
    for i_chi in even_theta_lv22: #odd theta=0
        i=i_chi[0]
        chi=i_chi[1]
        if lv22tnpsq[(i,chi)]==0:
            return True,False,i,chi #elliptic product with non-product theta.
    return False,0 #Jacobian.
   


