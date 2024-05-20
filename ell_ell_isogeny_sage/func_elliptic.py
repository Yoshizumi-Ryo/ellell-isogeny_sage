

#construct big finite field-------


#return finite field F_{p^4}.
def GFp4pow(p):
    Fp2=GF(p**2)
    z=Fp2(0)
    while z.is_square():
        z=Fp2.random_element()
    t=ZZ(z + z**p)
    n=ZZ(z**(p+1))
    R = PolynomialRing(GF(p),'x')
    x=R.gen()
    f=x**4-t*x**2+n
    Fp4=GF(p**4,modulus=f,name='a')
    return Fp4,x




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



def Elliptic_to_Legendre(E,base_field):
    coeff=E.a_invariants()
    #y^2=x^3+ax^2+bx+c
    a=coeff[2]
    b=coeff[3]
    c=coeff[4]
    K=E.base_field()
    assert(base_field==K)
    X=gen(K['X'])
    f=X**3+a*X**2+b*X+c
    roots=f.roots()
    assert(len(roots)==3)
    x1=roots[0][0]
    x2=roots[1][0]
    x3=roots[2][0]
    lmd=(x3-x1)/(x2-x1)
    E_lmd=EllipticCurve(K,[0,(-lmd-1),0,lmd,0])
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd




def Elliptic_to_Legendre_with_basis(E,N,base_field):
    K=E.base_field()
    assert(K==base_field)
    p=K.characteristic()
    assert((p+1)%N==0)
    Fp2=K.subfield(2)
    lmd,E_lmd,_=Elliptic_to_Legendre(E,Fp2)
    S1,S2=E_lmd.torsion_basis(N)
    lmd=K(lmd)
    E_lmd=EllipticCurve(K,[0,(-lmd-1),0,lmd,0])
    S1=E_lmd(S1)
    S2=E_lmd(S2)
    iso_E_Elmd=E.isomorphisms(E_lmd)[1]
    return lmd,E_lmd,iso_E_Elmd,S1,S2





def Elliptic_to_Legendre_with_basis_2(E,N,base_field):
    K=E.base_field()
    assert(K==base_field)
    p=K.characteristic()
    assert((p+1)%N==0)
    lmd,E_lmd,_=Elliptic_to_Legendre(E,K)
    S1,S2=E_lmd.torsion_basis(N)
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
    lv2tnp={"0":(thnp0_sq+thnp2_sq),
            "1":thnp1_sq}
    return lv2tnp,sq_rt_lm,sq_rt_lmm1




def Lv2tnp_to_Legendre(lv2tnp):
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
def Leg_lv2(lm,sq_rt_lm,sq_rt_lmm1,lv2tnp,pt):
    assert(sq_rt_lm**2==lm)
    assert(sq_rt_lmm1**2==lm-1)
    u=pt[0]
    v=pt[1]
    w=pt[2]
    if w==0:
        return lv2tnp
    assert(w==1)
    #assert((v**2)==u*(u-1)*(u-lm))
    thc0_sq=sq_rt_lm*(u-1)
    thc1_sq=sq_rt_lmm1*u
    thc2_sq=sq_rt_lm*thc0_sq-sq_rt_lmm1*thc1_sq
    thc3_sq=sq_rt_lm*thc1_sq-sq_rt_lmm1*thc0_sq
    lv2={"0":(parent(u))(1),
         "1":(thc1_sq+thc3_sq)/(thc0_sq+thc2_sq)}
    return lv2



def Lv2_to_Legendre(lm,sq_rt_lm,sq_rt_lmm1,lv2tc):
    r=lv2tc[1]/lv2tc[0]
    u=(sq_rt_lm*sq_rt_lmm1+r*(lm+sq_rt_lm))/(r*(1+sq_rt_lm)-sq_rt_lmm1)
    return u
    




def Many_Legendre_lv2(lm,sq_rt_lm,sq_rt_lmm1,lv2tnp,pts_set):
    lv2_coord={}
    for pt in pts_set:
        lv2_coord[pt]=Leg_lv2(lm,sq_rt_lm,sq_rt_lmm1,lv2tnp,pt)
    return lv2_coord

#--------------------------------------


def Product_theta(lv2_1,lv2_2):
    lv2=[lv2_1['0']*lv2_2['0'],
         lv2_1['0']*lv2_2['1'],
         lv2_1['1']*lv2_2['0'],
         lv2_1['1']*lv2_2['1']]
    return lv2

#------------------------------------



def Construct_pt(ell_data:list,points:list):
    assert(len(ell_data)==8)
    assert(len(points)==2)
    lm_1        =ell_data[0]
    sq_rt_lm_1  =ell_data[1]
    sq_rt_lmm1_1=ell_data[2]
    lv2tnp_E_1  =ell_data[3]
    
    lm_2        =ell_data[4]
    sq_rt_lm_2  =ell_data[5]
    sq_rt_lmm1_2=ell_data[6]
    lv2tnp_E_2  =ell_data[7]
    assert(sq_rt_lm_1**2  ==lm_1)
    assert(sq_rt_lmm1_1**2==lm_1-1)
    assert(sq_rt_lm_2**2  ==lm_2)
    assert(sq_rt_lmm1_2**2==lm_2-1)
    if (points[0]==0):
        lv2_x1=lv2tnp_E_1
    else:
        lv2_x1=Leg_lv2(lm_1,sq_rt_lm_1,sq_rt_lmm1_1,lv2tnp_E_1,points[0])
    if (points[1]==0):
        lv2_x2=lv2tnp_E_2
    else:
        lv2_x2=Leg_lv2(lm_2,sq_rt_lm_2,sq_rt_lmm1_2,lv2tnp_E_2,points[1])
    lv2_x12=Product_theta(lv2_x1,lv2_x2)
    tc_x12=Coord(lv2_x12,1)
    return tc_x12




#----------------------------------------------------------------



#cf.[LR16]section5
#sum_{t}(chi(t)lv2_x{i+t}lv2_x(t))
def product_term(lv2_x,i,chi):
    K=parent(lv2_x[0])
    sum=K(0)
    for t in range(0,4):
        chi_t=(-1)**((chi//2)*(t//2)+(chi%2)*(t%2))
        ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
        sum+=(chi_t*lv2_x[ipt]*lv2_x[t])
    return sum


#from lv2 theta null point to lv(2,2) theta null squard.
def lv2tnp_to_lv22tnpsq(lv2tnp):
    lv22tnpsq={}
    for chi in range(0,4):
        for i in range(0,4):
            lv22tnpsq[(chi,i)]=product_term(lv2tnp,i,chi)
    return lv22tnpsq
            


even_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==0}

odd_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==1}



#check if product of elliptic curves without theta structure.
def Is_Ellipitc_product(lv2tnp):
    lv22tnpsq=lv2tnp_to_lv22tnpsq(lv2tnp)
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
   


