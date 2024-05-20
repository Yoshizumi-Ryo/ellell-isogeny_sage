

#construct big finite field-------

#return finite field F_{p^4}.
def GFp4pow(p):
    Fp2=GF(p**2)
    z=Fp2(0)
    while z.is_square():
        z=Fp2.random_element()
    t=ZZ(z + z**p)
    n=ZZ(z**(p+1))
    X=gen(ZZ['X'])
    f=X**4-t*X**2+n
    Fp4=GF(p**4,modulus=f,name='a')
    return Fp4




#Cyclic isogeny--------------------


def Decomp_degree(N):
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
    coff=E.a_invariants()
    #y^2=x^3+ax^2+bx+c
    a=coff[2]
    b=coff[3]
    c=coff[4]
    K=E.base_ring()
    X=gen(K['X'])
    f=X**3+a*X**2+b*X+c
    roots=f.roots()
    x1=roots[0][0]
    x2=roots[1][0]
    x3=roots[2][0]
    lmd=(x3-x1)/(x2-x1)
    E_lmd=EllipticCurve(K,[0,(-lmd-1),0,lmd,0])
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd



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
    sqrt_lm=sqrt(lm)
    sqrt_lmm1=sqrt(lm-1)
    assert(sqrt_lm**2==lm)
    assert(sqrt_lmm1**2==lm-1)
    thnp0_sq=sqrt_lm
    thnp1_sq=sqrt_lmm1
    thnp2_sq=K(1)
    thnp3_sq=K(0)
    lv2tnp={"0":(thnp0_sq+thnp2_sq),
           "1":thnp1_sq}
    return lv2tnp,sqrt_lm,sqrt_lmm1


def lv2tnp_to_Legendre(lv2tnp):
    a=lv2tnp[0]
    b=lv2tnp[1]
    sqrt_lm=(a**2+b**2)/(a**2-b**2)
    sqrt_lmm1=(1+sqrt_lm)*(lv2tnp[1]/lv2tnp[0])
    lm=sqrt_lm**2
    assert(sqrt_lmm1**2+1==lm)
    return lm,sqrt_lm,sqrt_lmm1
    


def is_isomorphic_Legendre(lmd_1,lmd_2):
    if lmd_1 in {lmd_2,1/lmd_2,1-lmd_2,1/(1-lmd_2),1-(1/lmd_2),lmd_2/(lmd_2-1)}:
        return True
    else:
        return False


#-----------------------




#Legendre representaion to lv2 theta coordinate.
def Leg_lv2(lm,sqrt_lm,sqrt_lmm1,lv2tnp,pt):
    u=pt[0]
    v=pt[1]
    w=pt[2]
    if w==0:
        return lv2tnp
    assert(w==1)
    assert((v**2)==u*(u-1)*(u-lm))
    thc0_sq=sqrt_lm*(u-1)
    thc1_sq=sqrt_lmm1*u
    thc2_sq=sqrt_lm*thc0_sq-sqrt_lmm1*thc1_sq
    thc3_sq=sqrt_lm*thc1_sq-sqrt_lmm1*thc0_sq
    lv2={"0":1,
         "1":(thc1_sq+thc3_sq)/(thc0_sq+thc2_sq)}
    return lv2



def lv2_to_Legendre(lm,sqrt_lm,sqrt_lmm1,lv2tc):
    r=lv2tc[1]/lv2tc[0]
    u=(sqrt_lm*sqrt_lmm1+r*(lm+sqrt_lm))/(r*(1+sqrt_lm)-sqrt_lmm1)
    return u
    




def Many_Legendre_lv2(lm,sqrt_lm,sqrt_lmm1,lv2tnp,pts_set):
    lv2_coord={}
    for pt in pts_set:
        lv2_coord[pt]=Leg_lv2(lm,sqrt_lm,sqrt_lmm1,lv2tnp,pt)
    return lv2_coord

#--------------------------------------


def product_theta(lv2_1,lv2_2):
    lv2=[lv2_1['0']*lv2_2['0'],
         lv2_1['0']*lv2_2['1'],
         lv2_1['1']*lv2_2['0'],
         lv2_1['1']*lv2_2['1']]
    return lv2

