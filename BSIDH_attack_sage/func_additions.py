#functions of operations (e.g. Diff_Add) on Kummer surface


#correspondece of indecies. (0,0)=0, (0,1)=1, (1,0)=2, (1,1)=3.
#In this code, we use {0,1,2,3}=range(0,4)


even_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==0}

odd_theta_lv22={(i,j) for i,j in itertools.product(range(0,4),range(0,4)) if ((i//2)*(j//2)+(i%2)*(j%2))%2==1}



#useful functions=========================

#is equal as projective point.
def peq_lv2(lv2_tc1,lv2_tc2):
    assert(len(lv2_tc1)==4)
    assert(len(lv2_tc2)==4)
    assert(type(lv2_tc1)==list)
    assert(type(lv2_tc2)==list)
    assert(lv2_tc1[0]!=0)
    assert(lv2_tc2[0]!=0)
    ratio=lv2_tc1[0]/lv2_tc2[0]
    c=0
    for i in range(0,4):
        if (lv2_tc1[i]==(ratio*(lv2_tc2[i]))):
            c+=1
    if c==4:
        return True
    else:
        return False
       

  
      
       
       
#is equal as affine point.
def aeq_lv2(lv2_tc1,lv2_tc2):
    return lv2_tc1==lv2_tc2


#for const in base field.
def Mult_Const(const,lv2_x):
    new_lv2_x=[const*lv2_x[i] for i in range(0,4)]
    assert(peq_lv2(lv2_x,new_lv2_x))
    return new_lv2_x




#[DMPR23]=========================

def theta_square(lv2_x):
    sq_lv2_x=[lv2_x[0]**2,
              lv2_x[1]**2,
              lv2_x[2]**2,
              lv2_x[3]**2]
    return sq_lv2_x


def theta_Hadamard(lv2_x):
    x=lv2_x[0]
    y=lv2_x[1]
    z=lv2_x[2]
    w=lv2_x[3]
    dula_lv2_x=[x+y+z+w,
                x-y+z-w,
                x+y-z-w,
                x-y-z+w]
    return dula_lv2_x




#cf[LR16]section5.
def DiffAdd(lv2_x,lv2_y,lv2_xmy,HSN):
    HS_lv2_x=theta_Hadamard(theta_square(lv2_x))
    HS_lv2_y=theta_Hadamard(theta_square(lv2_y))
    z00_chi=[HS_lv2_x[i]*HS_lv2_y[i]/HSN[i] for i in range(0,4)]
    kappa_ii=theta_Hadamard(z00_chi)
    lv2_xpy=[kappa_ii[i]/(4*lv2_xmy[i]) for i in range(0,4)]
    return lv2_xpy



def Double(lv2tnp,lv2_x,HSN):
    HS_lv2_x=theta_Hadamard(theta_square(lv2_x))
    z00_chi=[(HS_lv2_x[i]**2)/HSN[i] for i in range(0,4)]
    kappa_ii=theta_Hadamard(z00_chi)
    lv2_2x=[kappa_ii[i]/(4*lv2tnp[i]) for i in range(0,4)]
    return lv2_2x
    
    
    

#calculate k*x.
def Mult(lv2tnp,k,lv2_x,HSN):
    if k==0:
        return lv2tnp
    if k==1:
        return lv2_x
    if k==2:
        return Double(lv2tnp,lv2_x,HSN)
    if k==3:
        return DiffAdd(Double(lv2tnp,lv2_x,HSN),lv2_x,lv2_x,HSN)
    bit_k=(ZZ(k)).digits(2)
    x=lv2_x
    y=Double(lv2tnp,lv2_x,HSN)
    xmy=lv2_x
    for i in range(1,len(bit_k)):
        bit_i=bit_k[len(bit_k)-i-1]
        if bit_i==1:
            x_0=DiffAdd(x,y,xmy,HSN)
            y_0=Double(lv2tnp,y,HSN)
            x=x_0
            y=y_0
        else:
            x_0=Double(lv2tnp,x,HSN)
            y_0=DiffAdd(x,y,xmy,HSN)
            x=x_0
            y=y_0
    return x




def Is_given_order(lv2tnp,order,lv2_x,HSN):
    assert(order>=1)
    if order==1:
        return peq_lv2(lv2tnp,lv2_x)
    for i in factor(order):
        s=order//i[0]
        if peq_lv2(lv2tnp,Mult(lv2tnp,s,lv2_x,HSN)):
            return False
    if peq_lv2(lv2tnp,Mult(lv2tnp,order,lv2_x,HSN)):
        return True
    return False







#Slow!Don't use.
def Mult_slow(lv2tnp,k,lv2_x,HSN):
    if k==0:
        return lv2tnp
    if k==1:
        return lv2_x
    if k>=2:
        t1=lv2_x
        t2=lv2tnp
        for i in range(2,k+1):
            t3=DiffAdd(t1,lv2_x,t2,HSN)
            t2=t1
            t1=t3
        return t3
    


#calculate kx+y from x,y,x-y.
def kxpy_xmy(lv2tnp,k,lv2_x,lv2_y,lv2_xmy,HSN):
    if k==0:
        return lv2_y
    bit_k=(ZZ(k)).digits(2)
    X=lv2_x
    Y=lv2_y
    Z=lv2_xmy
    for i in range(0,len(bit_k)):
        if bit_k[i]==1:
            X0=Double(lv2tnp,X,HSN)
            Y=DiffAdd(X,Y,Z,HSN)
            X=X0
        else:
            X0=Double(lv2tnp,X,HSN)
            Z=DiffAdd(X,Z,Y,HSN)
            X=X0
    return Y
    



#calculate kx+y from x,y,x+y.
def kxpy_xpy(lv2tnp,k,lv2_x,lv2_y,lv2_xpy,HSN):
    if k==0:
        return lv2_y
    if k==1:
        return lv2_xpy
    if k>=2:
        #kx+y=(k-1)x+(x+y).
        return kxpy_xmy(lv2tnp,(k-1),lv2_x,lv2_xpy,lv2_y,HSN)




#Riemann-relatioin==============


#construct all Riemann relation indecies.
RR_set={(i,j,k,l) for i in range(0,4) for j in range(0,4) for k in range(0,4) for l in range(0,4) if (i//2+j//2+k//2+l//2)%2==0 and (i%2+j%2+k%2+l%2)%2==0}



#cf.[LR16]section5
#sum_{t}(chi(t)lv2_x{i+t}lv2_x(t))
def product_term(lv2_x,i,chi):
    sum=0
    for t in range(0,4):
        chi_t=(-1)**((chi//2)*(t//2)+(chi%2)*(t%2))
        ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
        sum+=(chi_t*lv2_x[ipt]*lv2_x[t])
    return sum


#from lv2 theta null point to lv(2,2) theta null squard.
def lv2tnp_to_lv22tnpsq(lv2tnp):
    lv22tnpsq={}
    for i in range(0,4):
        for j in range(0,4):
            lv22tnpsq[(i,j)]=product_term(lv2tnp,j,i)
    return lv22tnpsq
            




#term appeard in Riemann relation.
def RR_product_term(lv2_x,lv2_y,chi,i,j):
    sum=0
    for t in range(0,4):
        chi_t=(-1)**((chi//2)*(t//2)+(chi%2)*(t%2))
        ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
        jpt=2*((j//2+t//2)%2)+(j%2+t%2)%2
        sum+=chi_t*lv2_x[ipt]*lv2_y[jpt]
    return sum



#check if given orderd 8 theta coordinates satisfy all Riemann relations.
def Is_RR(lv2_1,lv2_2,lv2_3,lv2_4,lv2_5,lv2_6,lv2_7,lv2_8):
    for chi in range(0,4):
        for ijkl in RR_set:
            LHS1=RR_product_term(lv2_1,lv2_2,chi,ijkl[0],ijkl[1])
            LHS2=RR_product_term(lv2_3,lv2_4,chi,ijkl[2],ijkl[3])
            RHS1=RR_product_term(lv2_5,lv2_6,chi,ijkl[0],ijkl[1])
            RHS2=RR_product_term(lv2_7,lv2_8,chi,ijkl[2],ijkl[3])
            if LHS1*LHS2!=RHS1*RHS2:
                return False
    return True

            

#we won't use.
def Normalize_for_RR_DiffAdd(lv2tnp,lv2_x,lv2_y,lv2_xpy,lv2_xmy):
    chi=0
    LHS_1=RR_product_term(lv2_xpy,lv2_xmy,chi,0,0)
    LHS_2=RR_product_term(lv2tnp ,lv2tnp ,chi,0,0)
    RHS_1=RR_product_term(lv2_y  ,lv2_y  ,chi,0,0)
    RHS_2=RR_product_term(lv2_x  ,lv2_x  ,chi,0,0)
    normalize_coff=(RHS_1*RHS_2)/(LHS_1*LHS_2)
    if normalize_coff==1:
        print("already excellent.")
        return lv2_xpy
    return Mult_Const(normalize_coff,lv2_xpy)





#check if product of elliptic curves without theta structure.
def Is_Ellipitc_product(lv2tnp):
    lv22tnpsq=lv2tnp_to_lv22tnpsq(lv2tnp)
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
   
   
    

#Three way addition.
def Extended_Addition(lv2tnp,lv2_x,lv2_y,lv2_z,lv2_xpy,lv2_ypz,lv2_zpx):
    LHS2=[RR_product_term(lv2_y,lv2_z,0,0,0),
          RR_product_term(lv2_y,lv2_z,1,0,0),
          RR_product_term(lv2_y,lv2_z,2,0,0),
          RR_product_term(lv2_y,lv2_z,3,0,0)]
    RHS2=[RR_product_term(lv2_zpx,lv2_xpy,0,0,0),
          RR_product_term(lv2_zpx,lv2_xpy,1,0,0),
          RR_product_term(lv2_zpx,lv2_xpy,2,0,0),
          RR_product_term(lv2_zpx,lv2_xpy,3,0,0)]   
    lv2_xpypz=[]
    for i in range(0,4):
        sum=0
        for chi in range(0,4):
            RHS1=RR_product_term(lv2tnp,lv2_ypz,chi,i,i)
            LHS1=(RHS1*RHS2[chi])/LHS2[chi]
            sum+=LHS1
        lv2_xpypz.append(sum/(4*lv2_x[i]))
    assert(Is_RR(lv2_xpypz,lv2_x,lv2_y,lv2_z,lv2tnp,lv2_ypz,lv2_zpx,lv2_xpy))
    return lv2_xpypz
        
            
    


#[cf.LR16]section5.
#from x,y,we compute set {x+y,x-y}. We can't distinguish x+y,x-y.
def Normal_Addition(lv2tnp,lv2_x,lv2_y):
    assert(not Is_Ellipitc_product(lv2tnp)[0])
    z_i_chi={i_chi:(product_term(lv2_x,i_chi[0],i_chi[1])*product_term(lv2_y,i_chi[0],i_chi[1])/product_term(lv2tnp,i_chi[0],i_chi[1])) for i_chi in even_theta_lv22}
    #kappa_{i,j}
    k_ij={(0,0):(z_i_chi[(0,0)]+z_i_chi[(0,1)]+z_i_chi[(0,2)]+z_i_chi[(0,3)]),
          (1,1):(z_i_chi[(0,0)]-z_i_chi[(0,1)]+z_i_chi[(0,2)]-z_i_chi[(0,3)]),
          (2,2):(z_i_chi[(0,0)]+z_i_chi[(0,1)]-z_i_chi[(0,2)]-z_i_chi[(0,3)]),
          (3,3):(z_i_chi[(0,0)]-z_i_chi[(0,1)]-z_i_chi[(0,2)]+z_i_chi[(0,3)]),
          (2,0):(z_i_chi[(2,0)]+z_i_chi[(2,1)]),
          (3,1):(z_i_chi[(2,0)]-z_i_chi[(2,1)]),
          (1,0):(z_i_chi[(1,0)]+z_i_chi[(1,2)]),
          (3,2):(z_i_chi[(1,0)]-z_i_chi[(1,2)]),
          (3,0):(z_i_chi[(3,0)]+z_i_chi[(3,3)]),
          (1,2):(z_i_chi[(3,0)]-z_i_chi[(3,3)])}
    
    for A in {(2,0),(3,1),(1,0),(3,2),(3,0),(1,2)}:
        k_ij[(A[1],A[0])]=k_ij[A]
    a=1 #alpha.
    rootin=k_ij[(a,0)]**2-(k_ij[(a,a)]*k_ij[(0,0)])
    sqrtrootin=sqrt(rootin)
    assert(rootin==sqrtrootin**2)
    Z=(k_ij[(a,0)]+sqrtrootin)/k_ij[(0,0)]
    Zd=2*(k_ij[(a,0)]/k_ij[(0,0)])-Z
    lv2_xpy=[(Z*k_ij[(0,i)]-k_ij[(a,i)])/(Z*k_ij[(0,0)]-k_ij[(a,0)]) for i in range(0,4)]
    lv2_xmy=[(Zd*k_ij[(0,i)]-k_ij[(a,i)])/(Zd*k_ij[(0,0)]-k_ij[(a,0)]) for i in range(0,4)]
    lv2_xpy=Normalize_for_RR_DiffAdd(lv2tnp,lv2_x,lv2_y,lv2_xpy,lv2_xmy)
    assert(Is_RR(lv2_xpy,lv2_xmy,lv2tnp,lv2tnp,lv2_y,lv2_y,lv2_x,lv2_x))
    return lv2_xpy,lv2_xmy


'''
#normal addition by solving 7-variables 17 quadratic simultaneous equation.
def Normal_Addition_2(lv2tnp,lv2_x,lv2_y):
    Poly.<X_1,X_2,X_3,Y_0,Y_1,Y_2,Y_3>=PolynomialRing(K,7)
    lv2_1=[1,X_1,X_2,X_3]
    lv2_2=[Y_0,Y_1,Y_2,Y_3]
    polyset=set()
    for chi in range(0,4):
        for ijkl in RR_set:
            LHS1=RR_product_term(lv2_1,lv2_2  ,chi,ijkl[0],ijkl[1])
            LHS2=RR_product_term(lv2tnp,lv2tnp,chi,ijkl[2],ijkl[3])
            RHS1=RR_product_term(lv2_y,lv2_y  ,chi,ijkl[0],ijkl[1])
            RHS2=RR_product_term(lv2_x,lv2_x  ,chi,ijkl[2],ijkl[3])
            polyset.add(LHS1*LHS2-RHS1*RHS2)
    polylist=list(polyset)
    print(len(polylist))
    I=Poly.ideal(polylist)
    B = I.groebner_basis()
    I=Poly.ideal(B)
    assert(I.dimension()==0)
    V=I.variety()
    sol=V[0]
    lv2_xpy=[K(1),sol[X_1],sol[X_2],sol[X_3]]
    lv2_xmy=[sol[Y_0],sol[Y_1],sol[Y_2],sol[Y_3]]
    return lv2_xpy,lv2_xmy
'''





#[LR16]
#calculate x+z from x,z,x+y,y+z.
def Compatible_Addition(lv2tnp,lv2_x,lv2_z,lv2_xpy,lv2_ypz):
    assert(not Is_Ellipitc_product(lv2tnp)[0])
    lv2_xpz_1,lv2_xpz_2=Normal_Addition(lv2tnp,lv2_x,lv2_z) #(x+z) or (x-z)
    #(x+y)+(x+z)=2x+y+z or y-z.
    #(x+y)+(x-z)=2x+y-z or y+z.
    lv2_tc1,lv2_tc2=Normal_Addition(lv2tnp,lv2_xpy,lv2_xpz_1)
    lv2_tc3,lv2_tc4=Normal_Addition(lv2tnp,lv2_xpy,lv2_xpz_2)
    if peq_lv2(lv2_ypz,lv2_tc1) or peq_lv2(lv2_ypz,lv2_tc2):
        return lv2_xpz_2
    elif peq_lv2(lv2_ypz,lv2_tc3) or peq_lv2(lv2_ypz,lv2_tc4):
        return lv2_xpz_1
    assert(False)





#we wont't use.
def Normalize_for_RR_Threeway(lv2tnp,lv2_x,lv2_y,lv2_z,lv2_xpy,lv2_xpz,lv2_ypz,lv2_xpypz):
    chi=0
    LHS_1=RR_product_term(lv2_xpypz,lv2_x  ,chi,0,0)
    LHS_2=RR_product_term(lv2_y    ,lv2_z  ,chi,0,0)
    RHS_1=RR_product_term(lv2tnp   ,lv2_ypz,chi,0,0)
    RHS_2=RR_product_term(lv2_xpz  ,lv2_xpy,chi,0,0)
    normalize_coff=(RHS_1*RHS_2)/(LHS_1*LHS_2)
    if normalize_coff==1:
        print("already excellent.")
        return lv2_xpypz
    return Mult_Const(normalize_coff,lv2_xpypz)



      

