


#Auxically path======================================



#compute all integers x such that x^2+1=0 mod M. (0<x<M)
def All_4throot_mod(M:int):
    assert(M>=2)
    #set of prime factors of M.
    pr_fac_M={list(factor(M))[i][0] for i in range(0,len(list(factor(M))))}
    if M%4==0:
        return set() #empty set.
    for pr in pr_fac_M:
        if pr%4==3:
            return set()  #empty set.
    fac_M=list(factor(M))
    num_pfac=len(fac_M)
    CRT=[fac_M[i][0]**fac_M[i][1] for i in range(0,num_pfac)]
    X_M=[ZZ(Mod(-1,CRT[i]).sqrt())for i in range(0,num_pfac)]
    set_All_4throot_modM=set()
    for bit in subsets(range(0,num_pfac)):
        Y_M=list()
        for i in range(0,num_pfac):
            if i in bit:
                sign=1
            else:
                sign=-1
            Y_M.append((sign*X_M[i])%CRT[i])
        set_All_4throot_modM.add(CRT_list(Y_M,CRT))
    return set_All_4throot_modM

        

#deteremine if M is B-smooth except for one prime factor.
def Smoothness(M:int,B:int):
    prime_under={pr for pr in range(2,B+1) if is_prime(pr)}
    for pr in prime_under:
        while (M%pr==0):
            M=M//pr
    if M==1 or is_prime(M):
        return True
    else:
        return False



#return r,s such that r^2+s^2=M where gcd(r,s)=1.
def Primitive_Cornacchia(M:int):
    if M==0:
        return True,0,0
    if M==1:
        return True,1,0
    if M==2:
        return True,1,1
    assert(M>=3)
    if not(Smoothness(M,2**(15))):
        return False,0,0
    set_All_4throot_modM=All_4throot_mod(M)
    for r0 in set_All_4throot_modM:
        if (2*r0)<=M:
            r_0=r0
            assert((r_0**2+1)%M==0)
            r_1=M%r_0
            while (r_1**2>=M):
                r_2=r_0%r_1
                r_0=r_1
                r_1=r_2
            assert(r_1**2 < M)
            if is_square(M-r_1**2):
                s=isqrt(M-r_1**2)
                assert(s**2==(M-r_1**2))
                assert(r_1**2+s**2==M)
                assert(GCD(r_1,s)==1)
                return True,r_1,s
    #Is you come here, you can't find.
    return False,0,0






def FullRepresentInteger(C:int,p:int):
    assert(C>p)
    assert(p%4==3)
    for i in range(1,2**(30)):
        cd=floor(sqrt(4*(C/p)))
        zd=randint(-cd,cd)
        cdd=floor(sqrt((4*(C/p))-zd**2))
        td=randint(-cdd,cdd)
        c=4*C-p*(zd**2+td**2)
        TF,xd,yd=Primitive_Cornacchia(c)
        if TF:
            assert(xd**2+yd**2==c)
            if ((xd-td)%2==0) and ((yd-zd)%2==0):
                x=(xd-td)//2
                y=(yd-zd)//2
                z=zd
                t=td
                deg=(xd**2+yd**2+p*(zd**2+td**2))//4
                assert(deg==C)
                return (x,y,z,t)
    #If you come here, you can't find the example.
    assert(False)





def Image_by_RepInt(E,rep_int,P,Q,end_i):
    K=E.base_field()
    E0p=EllipticCurve([K(1),K(0)])
    assert(E==E0p)
    end_j=EllipticCurveHom_frobenius(E)
    end_k=end_i*end_j
    x=rep_int[0]
    y=rep_int[1]
    z=rep_int[2]
    t=rep_int[3]
    hP=P.division_points(2)[0]
    hQ=Q.division_points(2)[0]
    assert(2*hP==P)
    assert(2*hQ==Q)    
    termP_3=(end_i(hP)+end_j(hP))
    termP_4=(hP+end_k(hP))
    termQ_3=(end_i(hQ)+end_j(hQ))
    termQ_4=(hQ+end_k(hQ))
    assert(order(hP)%order(termP_3)==0)
    assert(order(hP)%order(termP_4)==0)
    assert(order(hQ)%order(termQ_3)==0)
    assert(order(hQ)%order(termQ_4)==0)
    img_P=x*P+y*end_i(P)+z*termP_3+t*termP_4
    img_Q=x*Q+y*end_i(Q)+z*termQ_3+t*termQ_4
    assert(order(P)%order(img_P)==0)
    assert(order(Q)%order(img_Q)==0)
    return img_P,img_Q
    





#where P,Q are basis of E[N_A].
def Auxiliary_path(E,N_A:int,N_B:int,P,Q):
    K=E.base_field()
    E0p=EllipticCurve([K(1),K(0)])
    assert(E==E0p)
    assert(order(P)==N_A)
    assert(order(Q)==N_A)
    assert(N_A>N_B)
    assert(GCD(N_A,N_B)==1)
    a=N_A-N_B
    inv_N_B=(N_B).inverse_mod(N_A)
    assert((N_B*inv_N_B)%N_A==1)
    PdivN_B=inv_N_B*P
    QdivN_B=inv_N_B*Q
    assert(N_B*PdivN_B==P)
    assert(N_B*QdivN_B==Q)
    assert(order(PdivN_B)==N_A)
    assert(order(QdivN_B)==N_A)
    assert(a*N_B>p)
    aut_i=E.isomorphisms(E)[3]
    assert(negation_morphism(E)==aut_i**2)
    for counter in range(0,2**(10)):
        gamma_repint=FullRepresentInteger(a*N_B,p)
        gamma_PdivN_B,gamma_QdivN_B=Image_by_RepInt(E,gamma_repint,PdivN_B,QdivN_B,aut_i)
        assert(order(gamma_PdivN_B)==N_A)
        assert(order(gamma_QdivN_B)==N_A)
        S1,S2=E.torsion_basis(N_B)
        assert(order(S1)==N_B)
        assert(order(S2)==N_B)
        bkdd_1,bkdd_2=Image_by_RepInt(E,gamma_repint,S1,S2,aut_i)
        if order(bkdd_1)==N_B:
            bkdd=bkdd_1
            break
        elif order(bkdd_2)==N_B:
            bkdd=bkdd_2
            break
    assert(order(bkdd)==N_B)
    E_cd,alpha_P,alpha_Q=Elliptic_Cyclic(E,bkdd,gamma_PdivN_B,gamma_QdivN_B)
    assert(order(alpha_P)==N_A)
    assert(order(alpha_Q)==N_A)
    return E_cd,alpha_P,alpha_Q

#==========================================



#for given theta null point on product of ellipitc curves, then this function gives theta nill of product of theta.
#cf[DMPR23]p18.
def Theta_Split(lv2tnp:list,lv2_x:list,lv2_y:list,zeta_4):
    assert(zeta_4**2==-1)
    assert(Is_Ellipitc_product(lv2tnp)[0])
    if Is_Ellipitc_product(lv2tnp)[1]:  #already split.
        return lv2tnp,lv2_x,lv2_y
    else:
        i=Is_Ellipitc_product(lv2tnp)[2] #(i,j) is zero even theta.
        j=Is_Ellipitc_product(lv2tnp)[3]
    if (i,j)==(0,0):
        lv2tnp[2]*=zeta_4
        lv2tnp[3]*=zeta_4
        lv2_x[2] *=zeta_4
        lv2_x[3] *=zeta_4
        lv2_y[2] *=zeta_4
        lv2_y[3] *=zeta_4
    if Is_Ellipitc_product(lv2tnp)[1]:
        return lv2tnp,lv2_x,lv2_y
    else:
        i=Is_Ellipitc_product(lv2tnp)[2] #(i,j) is zero even theta.
        j=Is_Ellipitc_product(lv2tnp)[3]
    if i!=0 and j==0:
        assert(i!=0)
        lv2tnp=Theta_Hadamard(lv2tnp)
        lv2_x =Theta_Hadamard(lv2_x)
        lv2_y =Theta_Hadamard(lv2_y)        
    if Is_Ellipitc_product(lv2tnp)[1]:
        return lv2tnp,lv2_x,lv2_y
    else:
        i=Is_Ellipitc_product(lv2tnp)[2] #(i,j) is zero even theta.
        j=Is_Ellipitc_product(lv2tnp)[3]
    assert(j!=0)
    if j==1:
        (lv2tnp[1],lv2tnp[3])=(lv2tnp[3],lv2tnp[1])
        (lv2_x[1] ,lv2_x[3]) =(lv2_x[3] ,lv2_x[1])
        (lv2_y[1] ,lv2_y[3]) =(lv2_y[3] ,lv2_y[1])
    elif j==2:
        (lv2tnp[2],lv2tnp[3])=(lv2tnp[3],lv2tnp[2])
        (lv2_x[2] ,lv2_x[3]) =(lv2_x[3] ,lv2_x[2])
        (lv2_y[2] ,lv2_y[3]) =(lv2_y[3] ,lv2_y[2])
    if Is_Ellipitc_product(lv2tnp)[1]:
        return lv2tnp,lv2_x,lv2_y
    else:
        i=Is_Ellipitc_product(lv2tnp)[2] #(i,j) is zero even theta.
        j=Is_Ellipitc_product(lv2tnp)[3]
    assert(i==0 and j==3)
    lv2tnp[1]*=zeta_4
    lv2tnp[2]*=zeta_4
    lv2_x[1] *=zeta_4
    lv2_x[2] *=zeta_4
    lv2_y[1] *=zeta_4
    lv2_y[2] *=zeta_4
    assert(Is_Ellipitc_product(lv2tnp)[1])
    return lv2tnp,lv2_x,lv2_y
    
    
#---------------------

def Is_correct_isogeny(E_dm,E_cd,deg:int,ker,P_dm,Q_dm,P_cd,Q_cd):
    assert(order(ker)==deg)
    E_cd_1,P_cd_1,Q_cd_1=Elliptic_Cyclic(E_dm,ker,P_dm,Q_dm)
    assert(parent(E_cd.j_invariant())==parent(E_cd_1.j_invariant()))
    if E_cd.j_invariant()!=E_cd_1.j_invariant():
        return False,0
    assert(E_cd.j_invariant()==E_cd_1.j_invariant())
    assert(len(E_cd_1.isomorphisms(E_cd))==2)
    iso_Ecd1_Ecd=E_cd_1.isomorphisms(E_cd)[0]
    assert(iso_Ecd1_Ecd(P_cd_1) in E_cd)
    assert(P_cd in E_cd)
    assert(iso_Ecd1_Ecd(Q_cd_1) in E_cd)
    assert(Q_cd in E_cd)
    if   (iso_Ecd1_Ecd(P_cd_1)==P_cd)  and (iso_Ecd1_Ecd(Q_cd_1)==Q_cd):
        return True,ker
    elif (iso_Ecd1_Ecd(P_cd_1)==-P_cd) and (iso_Ecd1_Ecd(Q_cd_1)==-Q_cd):
        return True,ker
    else:
        return False,0
    
    
    

#============================================



def Construct_pt(ell_data:list,points:list,K):
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
    tc_x12=Coord([K(lv2_x12[i]) for i in range(0,4)],K(1))
    return tc_x12


#----------------------------------------------------------------


#Main attack function.

def Attack_prepare(E_B,E_pr,N_A:int,N_B:int,PA_EB,QA_EB,alpha_PA,alpha_QA,K):
    assert(PA_EB    in E_B)
    assert(QA_EB    in E_B)
    assert(alpha_PA in E_pr)
    assert(alpha_QA in E_pr)
    assert(order(PA_EB)   ==N_A)
    assert(order(QA_EB)   ==N_A)
    assert(order(alpha_PA)==N_A)
    assert(order(alpha_QA)==N_A)
    assert(gcd(N_A,N_B)==1)
    #assert(PA_EB.weil_pairing(QA_EB,N_A)*alpha_PA.weil_pairing(alpha_QA,N_A)==1)
    lm_pr,E_pr_L,iso_Epr_EprL=Elliptic_to_Legendre(E_pr)
    lm_B,E_B_L,iso_EB_EBL,S1,S2=Elliptic_to_Legendre_with_basis(E_B,N_B)
    #print("take basis fin")
    #S1,S2=E_B_L.torsion_basis(N_B)
    assert(order(S1)==N_B)
    assert(order(S2)==N_B)
    alpha_PA=iso_Epr_EprL(alpha_PA)
    alpha_QA=iso_Epr_EprL(alpha_QA)
    PA_EB   =iso_EB_EBL  (PA_EB)
    QA_EB   =iso_EB_EBL  (QA_EB)
    l=Decomp_degree(N_A)[0]
    k=N_A//l
    assert(is_prime(l))
    assert(is_odd(l))
    #theta on E_pr*E_B---------
    #null
    lv2tnp_Epr,sq_rt_lm_pr,sq_rt_lmm1_pr=Legendre_to_lv2tnp(lm_pr)
    lv2tnp_EB ,sq_rt_lm_B ,sq_rt_lmm1_B =Legendre_to_lv2tnp(lm_B)
    lv2tnp=Product_theta(lv2tnp_Epr,lv2tnp_EB)
    tc_0=NullCoord([K(lv2tnp[i]) for i in range(0,4)],K(1))
    ell_data=[lm_pr,sq_rt_lm_pr,sq_rt_lmm1_pr,lv2tnp_Epr,lm_B ,sq_rt_lm_B ,sq_rt_lmm1_B ,lv2tnp_EB]
    f_1=[alpha_PA ,PA_EB]
    f_2=[alpha_QA ,QA_EB]  
    x  =[E_pr_L(0),S1   ]
    y  =[E_pr_L(0),S2   ]
    assert(order(k*PA_EB)==l)
    assert(order(k*QA_EB)==l)
    tc_f1   =Construct_pt(ell_data,f_1                          ,K)
    tc_f2   =Construct_pt(ell_data,f_2                          ,K)
    tc_f12  =Construct_pt(ell_data,[f_1[0]+f_2[0],f_1[1]+f_2[1]],K)
    tc_x    =Construct_pt(ell_data,x                            ,K)
    tc_xpf1 =Construct_pt(ell_data,[x[0]  +f_1[0],x[1]  +f_1[1]],K)
    tc_xpf2 =Construct_pt(ell_data,[x[0]  +f_2[0],x[1]  +f_2[1]],K)
    tc_y    =Construct_pt(ell_data,y                            ,K)
    tc_ypf1 =Construct_pt(ell_data,[y[0]  +f_1[0],y[1]  +f_1[1]],K)
    tc_ypf2 =Construct_pt(ell_data,[y[0]  +f_2[0],y[1]  +f_2[1]],K)
    assert(Is_Ellipitc_product(tc_0.numer)[0])
    #assert(tc_0.Is_order(tc_f1  ,N_A))
    #assert(tc_0.Is_order(tc_f2  ,N_A))
    #assert(tc_0.Is_order(tc_f12 ,N_A))
    #assert(tc_0.Is_order(tc_x   ,N_B))
    #assert(tc_0.Is_order(tc_xpf1,N_A*N_B))
    #assert(tc_0.Is_order(tc_xpf2,N_A*N_B))
    #assert(tc_0.Is_order(tc_y   ,N_B))
    #assert(tc_0.Is_order(tc_ypf1,N_A*N_B))
    #assert(tc_0.Is_order(tc_ypf2,N_A*N_B))
    return E_B_L,PA_EB,QA_EB,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,tc_y,tc_ypf1,tc_ypf2








def Attack_main(N_A:int,N_B:int,tc_0:NullCoord,tc_f1:Coord,tc_f2:Coord,tc_f12:Coord,tc_x:Coord,tc_xpf1:Coord,tc_xpf2:Coord,tc_y:Coord,tc_ypf1:Coord,tc_ypf2:Coord):
    assert(gcd(N_A,N_B)==1)
    assert(tc_0.Is_order(tc_f1  ,N_A))
    assert(tc_0.Is_order(tc_f2  ,N_A))
    assert(tc_0.Is_order(tc_f12 ,N_A))
    assert(tc_0.Is_order(tc_x   ,N_B))
    assert(tc_0.Is_order(tc_xpf1,N_A*N_B))
    assert(tc_0.Is_order(tc_xpf2,N_A*N_B))
    assert(tc_0.Is_order(tc_y   ,N_B))
    assert(tc_0.Is_order(tc_ypf1,N_A*N_B))
    assert(tc_0.Is_order(tc_ypf2,N_A*N_B))
    assert(tc_0.Is_same_proj(tc_0.Mult(tc_0,2)))
    fac=Decomp_degree(N_A)
    print("isogeny chain:",fac)
    s=1
    for i in range(0,len(fac)):
        l=fac[i]
        assert(is_prime(l))
        k=N_A//(s*l)
        assert(s*l*k==N_A)
        print("ell=",l)
        if i!=0: #if not the first step.
            tc_f12=tc_0.Normal_Add(tc_f1,tc_f2,1)
        #construct kernel e_1,e_2,e_1+e_2.
        tc_e1 =tc_0.Mult(tc_f1 ,k)
        tc_e2 =tc_0.Mult(tc_f2 ,k)
        tc_e12=tc_0.Mult(tc_f12,k)
        tc_e1.order=l
        tc_e2.order=l
        tc_e12.order=l
        tc_f1pe1=tc_0.Mult(tc_f1,k+1)               #f_1+e_1
        tc_f1pe2=tc_0.Kxpy_xpy(k,tc_f2,tc_f1,tc_f12)#f_1+e_2
        tc_f2pe1=tc_0.Kxpy_xpy(k,tc_f1,tc_f2,tc_f12)#f_2+e_1
        tc_f2pe2=tc_0.Mult(tc_f2,k+1)               #f_2+e_2
        #assert(tc_0.Is_order(tc_e1  ,l))
        #assert(tc_0.Is_order(tc_e2  ,l))
        #assert(tc_0.Is_order(tc_e12 ,l))
        #assert(tc_0.Is_order(tc_f1  ,k*l))
        #assert(tc_0.Is_order(tc_f2  ,k*l))
        #assert(tc_0.Is_order(tc_f12 ,k*l))
        #assert(tc_0.Is_order(tc_x   ,N_B))
        #assert(tc_0.Is_order(tc_y   ,N_B))
        if i==0:
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            tc_ype1 =tc_0.Kxpy_xpy(k,tc_f1,tc_y,tc_ypf1)#y+e_1
            tc_ype2 =tc_0.Kxpy_xpy(k,tc_f2,tc_y,tc_ypf2)#y+e_2
        if i!=0:
            tc_xpe1 =tc_0.Normal_Add(tc_x,tc_e1,1)#x+e_1
            tc_xpe2 =tc_0.Compatible_Add(tc_x,tc_e2,tc_xpe1,tc_e12)#x+e_2
            tc_ype1 =tc_0.Normal_Add(tc_y,tc_e1,1)#y+e_1
            tc_ype2 =tc_0.Compatible_Add(tc_y,tc_e2,tc_ype1,tc_e12)#y+e_2
        #assert(tc_0.Is_order(tc_xpe1,l*N_B))
        #assert(tc_0.Is_order(tc_xpe2,l*N_B))
        #assert(tc_0.Is_order(tc_ype1,l*N_B))
        #assert(tc_0.Is_order(tc_ype2,l*N_B))
        basis  =[tc_e1,tc_e2   ,tc_e12  ]
        f1_list=[tc_f1,tc_f1pe1,tc_f1pe2]
        f2_list=[tc_f2,tc_f2pe1,tc_f2pe2]
        x_list =[tc_x ,tc_xpe1 ,tc_xpe2 ]
        y_list =[tc_y ,tc_ype1 ,tc_ype2 ]
        #-----------------------------------------------------------------------------------------
        tc_cd0=CodOne(tc_0,basis)
        lmd_data=Product_power_lambda(basis)
        tc_f1=EvalOne(tc_0,basis,f1_list,lmd_data)
        tc_f2=EvalOne(tc_0,basis,f2_list,lmd_data)
        tc_x =EvalOne(tc_0,basis,x_list ,lmd_data)
        tc_y =EvalOne(tc_0,basis,y_list ,lmd_data)
        tc_0 =tc_cd0
        s*=l
        #assert(tc_0.Is_order(tc_x,N_B))
        #assert(tc_0.Is_order(tc_y,N_B))
        #assert(tc_0.Is_order(tc_f1,k))
        #assert(tc_0.Is_order(tc_f2,k))
    #=====================================
    assert(Is_Ellipitc_product(tc_0.numer)[0])
    assert(tc_0.Is_same_proj(tc_f1))
    assert(tc_0.Is_same_proj(tc_f1))
    return tc_0.numer,tc_x.numer,tc_y.numer








def Attack_last(E_0m,N_A:int,N_B:int,P_A,Q_A,PA_EB,QA_EB,E_B_L,lv2tnp:list,lv2_ZS1:list,lv2_ZS2:list,zeta_4):
    lv2tnp,lv2_ZS1,lv2_ZS2=Theta_Split(lv2tnp,lv2_ZS1,lv2_ZS2,zeta_4) 
    assert(Is_Ellipitc_product(lv2tnp)[0])
    assert(Is_Ellipitc_product(lv2tnp)[1])
    assert(lv2tnp[0]*lv2tnp[3]==lv2tnp[1]*lv2tnp[2])
    #split to elliptic curves.
    lv2tnp_Ecd_1=[lv2tnp[0],lv2tnp[1]]
    lv2tnp_Ecd_2=[lv2tnp[0],lv2tnp[2]]
    lm_1=Lv2tnp_to_Legendre(lv2tnp_Ecd_1)[0]
    lm_2=Lv2tnp_to_Legendre(lv2tnp_Ecd_2)[0]
    
    assert(P_A in E_0m)
    assert(Q_A in E_0m)
    assert(order(P_A)==N_A)
    assert(order(Q_A)==N_A)
    lm_0=Elliptic_to_Legendre(E_0m)[0]
    
    if Is_isomorphic_Legendre(lm_0,lm_1):
        lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd=Lv2tnp_to_Legendre(lv2tnp_Ecd_1)
        lv2_P1_cd=[lv2_ZS1[0],lv2_ZS1[1]]
        lv2_P2_cd=[lv2_ZS2[0],lv2_ZS2[1]]
        u_P1_cd=Lv2_to_Legendre(lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd,lv2_P1_cd)
        u_P2_cd=Lv2_to_Legendre(lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd,lv2_P2_cd)
    elif Is_isomorphic_Legendre(lm_0,lm_2):
        lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd=Lv2tnp_to_Legendre(lv2tnp_Ecd_2)
        lv2_P1_cd=[lv2_ZS1[0],lv2_ZS1[2]]
        lv2_P2_cd=[lv2_ZS2[0],lv2_ZS2[2]]
        u_P1_cd=Lv2_to_Legendre(lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd,lv2_P1_cd)
        u_P2_cd=Lv2_to_Legendre(lm_cd,sq_rt_lm_cd,sq_rt_lmm1_cd,lv2_P2_cd)
    else:
        assert(False)
    E_cd=Legendre_to_Elliptic(lm_cd)
    v_P1_cd=sqrt(u_P1_cd*(u_P1_cd-1)*(u_P1_cd-lm_cd))
    P1=E_cd([u_P1_cd,v_P1_cd,1])
    v_P2_cd=sqrt(u_P2_cd*(u_P2_cd-1)*(u_P2_cd-lm_cd))
    P2=E_cd([u_P2_cd,v_P2_cd,1])
    assert(N_B % order(P1)==0)
    assert(N_B % order(P2)==0)
    #since phi_B is cyclic isogeny, we get one element P generating the kernel.
    assert(P1.weil_pairing(P2,N_B)==1)
    if order(P1)==N_B:
        P=P1
    elif order(P2)==N_B:
        P=P2
    else:
        for coff in range(1,N_B):
            R=P1+coff*P2
            if order(R)==N_B:
                P=R
                break
    assert(order(P)==N_B)
    assert(P1.weil_pairing(P,N_B)==1)
    assert(P2.weil_pairing(P,N_B)==1)   
    #we have to consider the potential of automorphism of E_0.
    iso_Ecd_E0m=E_cd.isomorphisms(E_0m)[1] #isomorphism  E_cd->E_0m.
    ker=iso_Ecd_E0m(P)
    assert(ker in E_0m)
    for aut in E_0m.automorphisms():
        iscorrect=Is_correct_isogeny(E_0m,E_B_L,N_B,aut(ker),P_A,Q_A,PA_EB,QA_EB)
        if iscorrect[0]:
            attacker_ker=iscorrect[1]
    assert(attacker_ker in E_0m)
    return attacker_ker
    
#------------------------





def Attack_total(E_0m,E_B,E_pr,N_A:int,N_B:int,P_A,Q_A,PA_EB,QA_EB,alpha_PA,alpha_QA,zeta_4,K):
    assert(order(P_A)==N_A)
    assert(order(Q_A)==N_A)
    assert(order(PA_EB)==N_A)
    assert(order(QA_EB)==N_A)
    assert(order(alpha_PA)==N_A)
    assert(order(alpha_QA)==N_A)
    assert(zeta_4**2==-1)
    assert(gcd(N_A,N_B)==1)
    #prepare
    E_B_L,PA_EB,QA_EB,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,tc_y,tc_ypf1,tc_ypf2=Attack_prepare(E_B,E_pr,N_A,N_B,PA_EB,QA_EB,alpha_PA,alpha_QA,K)
    #isogeny
    lv2tnp,lv2_ZS1,lv2_ZS2=Attack_main(N_A,N_B,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,tc_y,tc_ypf1,tc_ypf2)
    #last
    attacker_ker=Attack_last(E_0m,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,E_B_L,lv2tnp,lv2_ZS1,lv2_ZS2,zeta_4)
    return attacker_ker




def BSIDH_construct_attack(p,N_A,N_B):
    assert((p-1)%N_A==0)
    assert((p+1)%N_B==0)
    #the 4th primitive root of 1=========
    K,_,Fp2=GFp4pow(p)
    X=gen(K['X'])
    f=X**2+1
    zeta_4=f.roots()[0][0]
    assert(zeta_4**2==-1)
    #===================================
    #public setting.====================
    #print("public setting")
    #y^2=x(x-1)(x+1)=x^3-x.
    E_0m=EllipticCurve(K,[K(-1),K(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    P_B,Q_B=E_0m.torsion_basis(N_B)
    #==================================
    #Bob calculates secretly.===========
    #print("Bob calculates secretly.")
    coff_B=randint(0,(N_B-1))
    R_B=P_B+coff_B*Q_B
    assert(order(R_B)==N_B)
    E_B,PA_EB,QA_EB=Elliptic_Cyclic(E_0m,R_B,P_A,Q_A)
    assert(order(PA_EB)==N_A)
    assert(order(QA_EB)==N_A)
    #===================================
    time_start=time.time()
    #Construction auxiliary path.=========
    #print("Construction auxiliary path.")
    E_0p=EllipticCurve(K,[K(1),K(0)])
    iso_mtop=E_0m.isomorphisms(E_0p)[1]
    E_pr,alpha_PA,alpha_QA=Auxiliary_path(E_0p,N_A,N_B,iso_mtop(P_A),iso_mtop(Q_A))
    #=====================================
    #Main attack.=========================
    #print("Main attack.")
    Attacker_ker=Attack_total(E_0m,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_PA,alpha_QA,zeta_4,K)
    #=====================================
    time_end=time.time()
    time_diff=time_end-time_start
    #Check if this attacke succeed.========
    assert(order(Attacker_ker)==N_B)
    result=(Attacker_ker.weil_pairing(R_B,N_B)==1)
    #=====================================
    return time_diff,result
    


