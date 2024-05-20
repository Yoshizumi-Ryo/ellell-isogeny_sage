


#Auxically path======================================



#compute all integers x such that x^2+1=0 mod M. (0<x<M)
def all_4throot_mod(M):
    assert(M>=2)
    divM=list(factor(M))
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
    set_all_4throot_modM=set()
    for bit in subsets(range(0,num_pfac)):
        Y_M=list()
        for i in range(0,num_pfac):
            if i in bit:
                sign=1
            else:
                sign=-1
            Y_M.append((sign*X_M[i])%CRT[i])
        set_all_4throot_modM.add(CRT_list(Y_M,CRT))
    return set_all_4throot_modM

        

#deteremine if M is B-smooth except for one prime factor.
def Smoothness(M,B):
    prime_under={pr for pr in range(2,B+1) if is_prime(pr)}
    for pr in prime_under:
        while (M%pr==0):
            M=M//pr
    if M==1 or is_prime(M):
        return True
    else:
        return False



#return r,s such that r^2+s^2=M where gcd(r,s)=1.
def primitive_Cornacchia(M):
    if M==0:
        return True,0,0
    if M==1:
        return True,1,0
    if M==2:
        return True,1,1
    assert(M>=3)
    if not(Smoothness(M,2**(15))):
        return False,0,0
    set_all_4throot_modM=all_4throot_mod(M)
    for r0 in set_all_4throot_modM:
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






def FullRepresentInteger(C,p):
    assert(C>p)
    assert(p%4==3)
    for i in range(1,2**(30)):
        cd=floor(sqrt(4*(C/p)))
        zd=randint(-cd,cd)
        cdd=floor(sqrt((4*(C/p))-zd**2))
        td=randint(-cdd,cdd)
        c=4*C-p*(zd**2+td**2)
        TF,xd,yd=primitive_Cornacchia(c)
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
    K=E.base_ring()
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
def Auxiliary_path(E,N_A,N_B,P,Q):
    K=E.base_ring()
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
def Theta_Split(lv2tnp,lv2_x,lv2_y,zeta_4):
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
        lv2tnp=theta_Hadamard(lv2tnp)
        lv2_x =theta_Hadamard(lv2_x)
        lv2_y =theta_Hadamard(lv2_y)        
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

def Is_correct_isogeny(E_dm,E_cd,deg,ker,P_dm,Q_dm,P_cd,Q_cd):
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

#Main attack function.



def attack_prepare(E_B,E_pr,N_A,N_B,PA_EB,QA_EB,alpha_PA,alpha_QA):
    assert(PA_EB    in E_B)
    assert(QA_EB    in E_B)
    assert(alpha_PA in E_pr)
    assert(alpha_QA in E_pr)
    assert(order(PA_EB)   ==N_A)
    assert(order(QA_EB)   ==N_A)
    assert(order(alpha_PA)==N_A)
    assert(order(alpha_QA)==N_A)
    lm_pr,E_pr_L,iso_Epr_EprL=Elliptic_to_Legendre(E_pr)
    #lm_B,E_B_L,iso_EB_EBL    =Elliptic_to_Legendre(E_B)
    #print("take bais")
    lm_B,E_B_L,iso_EB_EBL,S1,S2=Elliptic_to_Legendre_with_basis(E_B,N_B)
    #print("take bais fin")
    #S1,S2=E_B_L.torsion_basis(N_B)
    assert(order(S1)==N_B)
    assert(order(S2)==N_B)
    alpha_PA=iso_Epr_EprL(alpha_PA)
    alpha_QA=iso_Epr_EprL(alpha_QA)
    PA_EB=iso_EB_EBL(PA_EB)
    QA_EB=iso_EB_EBL(QA_EB)
    #theta on E_pr*E_B---------
    #null
    lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr=Legendre_to_lv2tnp(lm_pr)
    lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B =Legendre_to_lv2tnp(lm_B)
    lv2tnp=product_theta(lv2tnp_Epr,lv2tnp_EB)
    #f_1=(alpha_PA,PA_EB)
    lv2_alpha_PA=Leg_lv2(lm_pr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_Epr,alpha_PA)
    lv2_f1_EB   =Leg_lv2(lm_B ,sqrt_lm_B ,sqrt_lmm1_B ,lv2tnp_EB ,PA_EB)
    lv2_f1=product_theta(lv2_alpha_PA,lv2_f1_EB)
    #f_2=(alpha_QA,QA_EB)
    lv2_alpha_QA=Leg_lv2(lm_pr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_Epr,alpha_QA)
    lv2_f2_EB   =Leg_lv2(lm_B ,sqrt_lm_B ,sqrt_lmm1_B ,lv2tnp_EB ,QA_EB)
    lv2_f2=product_theta(lv2_alpha_QA,lv2_f2_EB)
    #f_12=f_1+f_2
    lv2_alpha_PQA=Leg_lv2(lm_pr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_Epr,alpha_PA+alpha_QA)
    lv2_f12_EB   =Leg_lv2(lm_B ,sqrt_lm_B ,sqrt_lmm1_B ,lv2tnp_EB ,PA_EB+QA_EB)
    lv2_f12=product_theta(lv2_alpha_PQA,lv2_f12_EB)
    #ZS1=(0,S1)
    lv2_S1=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S1)
    lv2_ZS1=product_theta(lv2tnp_Epr,lv2_S1)
    #ZS1+f1=(alpha_PA,S1+PA_EB)
    lv2_S1_PAEB=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S1+PA_EB)
    lv2_ZS1pf1=product_theta(lv2_alpha_PA,lv2_S1_PAEB)
    #ZS1+f2=(alpha_QA,S1+QA_EB)
    lv2_S1_QAEB=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S1+QA_EB)
    lv2_ZS1pf2=product_theta(lv2_alpha_QA,lv2_S1_QAEB)
    #ZS2=(0,S2)
    lv2_S2=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S2)
    lv2_ZS2=product_theta(lv2tnp_Epr,lv2_S2)
    #ZS2+f1=(alpha_PA,S2+PA_EB)
    lv2_S2_PAEB=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S2+PA_EB)
    lv2_ZS2pf1=product_theta(lv2_alpha_PA,lv2_S2_PAEB)
    #ZS2+f2=(alpha_QA,S2+QA_EB)
    lv2_S2_QAEB=Leg_lv2(lm_B,sqrt_lm_B,sqrt_lmm1_B,lv2tnp_EB,S2+QA_EB)
    lv2_ZS2pf2=product_theta(lv2_alpha_QA,lv2_S2_QAEB)
    #----------------------------
    return E_B_L,PA_EB,QA_EB,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2,lv2_ZS2,lv2_ZS2pf1,lv2_ZS2pf2
  





def attack_isogeny(N_A,N_B,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2,lv2_ZS2,lv2_ZS2pf1,lv2_ZS2pf2,isogeny_type):
    print("begin with calculation attack isogeny")
    fac=Decomp_degree(N_A)
    print(fac)
    s=1
    for i in range(0,len(fac)):
        time_start=time.time()
        l=fac[i]
        k=N_A//(s*l)
        assert(s*l*k==N_A)
        print("l=",l)
        #print("elliptic product",Is_Ellipitc_product(lv2tnp)[0])
        HSN=theta_Hadamard(theta_square(lv2tnp))
        assert(Is_given_order(lv2tnp,N_B,lv2_ZS1,HSN))
        #assert(peq_lv2((Mult(lv2tnp,l*k,lv2_f1,HSN)),lv2tnp))
        #assert(peq_lv2((Mult(lv2tnp,l*k,lv2_f2,HSN)),lv2tnp))
        if i!=0:
            assert(not Is_Ellipitc_product(lv2tnp)[0])
            lv2_f12   =Normal_Addition(lv2tnp,lv2_f1,lv2_f2)[0]
            lv2_ZS1pf1=Normal_Addition(lv2tnp,lv2_ZS1,lv2_f1)[0]
            lv2_ZS1pf2=Compatible_Addition(lv2tnp,lv2_ZS1,lv2_f2,lv2_ZS1pf1,lv2_f12)
            lv2_ZS2pf1=Normal_Addition(lv2tnp,lv2_ZS2,lv2_f1)[0]
            lv2_ZS2pf2=Compatible_Addition(lv2tnp,lv2_ZS2,lv2_f2,lv2_ZS2pf1,lv2_f12)
        #--------------------------------
        #construct kernel.
        lv2_e1 =Mult(lv2tnp,k,lv2_f1 ,HSN)
        lv2_e2 =Mult(lv2tnp,k,lv2_f2 ,HSN)
        lv2_e12=Mult(lv2tnp,k,lv2_f12,HSN)
        #linear combination of f_1,f_2.
        lincom_f1f2={}
        lincom_f1f2[((k+1),0)]=Mult(lv2tnp,k+1,lv2_f1,HSN)
        lincom_f1f2[(0,(k+1))]=Mult(lv2tnp,k+1,lv2_f2,HSN)
        lincom_f1f2[(k,1)]    =kxpy_xpy(lv2tnp,k,lv2_f1,lv2_f2,lv2_f12,HSN)
        lincom_f1f2[(1,k)]    =kxpy_xpy(lv2tnp,k,lv2_f2,lv2_f1,lv2_f12,HSN)
                
        #ZS1+lincom of f_1,f_2.
        tc_ZS1_lincomf1f2={}
        tc_ZS1_lincomf1f2[(0,0)]=lv2_ZS1 #ZS1
        tc_ZS1_lincomf1f2[(k,0)]=kxpy_xpy(lv2tnp,k,lv2_f1,lv2_ZS1,lv2_ZS1pf1,HSN)  #ZS1+e_1
        tc_ZS1_lincomf1f2[(0,k)]=kxpy_xpy(lv2tnp,k,lv2_f2,lv2_ZS1,lv2_ZS1pf2,HSN)  #ZS1+e_2
        #ZS2+lincom of f_1,f_2.
        tc_ZS2_lincomf1f2={}
        tc_ZS2_lincomf1f2[(0,0)]=lv2_ZS2  #ZS2
        tc_ZS2_lincomf1f2[(k,0)]=kxpy_xpy(lv2tnp,k,lv2_f1,lv2_ZS2,lv2_ZS2pf1,HSN) #ZS2+e_1
        tc_ZS2_lincomf1f2[(0,k)]=kxpy_xpy(lv2tnp,k,lv2_f2,lv2_ZS2,lv2_ZS2pf2,HSN)  #ZS2+e_2
        #=================================
        if isogeny_type=="tl_lpow":
            print("codomain")
            #codomian 
            lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12=codomain_lv2tnp_tl(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
            #-----------------------------------
            #evaluation
            #f_1
            print("eveluation 1/4")
            lv2_f1_cd=evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f1,lincom_f1f2[((k+1),0)],lincom_f1f2[(1,k)],HSN)
            #f_2
            print("eveluation 2/4")
            lv2_f2_cd=evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f2,lincom_f1f2[(k,1)],lincom_f1f2[(0,(k+1))],HSN)
            #ZS1
            print("eveluation 3/4")
            lv2_ZS1_cd=evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS1,tc_ZS1_lincomf1f2[(k,0)],tc_ZS1_lincomf1f2[(0,k)],HSN)
            #ZS2
            print("eveluation 4/4")
            lv2_ZS2_cd=evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS2,tc_ZS2_lincomf1f2[(k,0)],tc_ZS2_lincomf1f2[(0,k)],HSN)
        #-----------------------------
        if isogeny_type=="tl_sumsq":
            print("codomain")
            #codomian 
            lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12=codomain_lv2tnp_tl_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
            #-----------------------------------
            #evaluation
            #f_1
            print("eveluation 1/4")
            lv2_f1_cd=evaluation_lv2_tl_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f1,lincom_f1f2[((k+1),0)],lincom_f1f2[(1,k)],HSN)
            #f_2
            print("eveluation 2/4")
            lv2_f2_cd=evaluation_lv2_tl_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f2,lincom_f1f2[(k,1)],lincom_f1f2[(0,(k+1))],HSN)
            #ZS1
            print("eveluation 3/4")
            lv2_ZS1_cd=evaluation_lv2_tl_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS1,tc_ZS1_lincomf1f2[(k,0)],tc_ZS1_lincomf1f2[(0,k)],HSN)
            #ZS2
            print("eveluation 4/4")
            lv2_ZS2_cd=evaluation_lv2_tl_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS2,tc_ZS2_lincomf1f2[(k,0)],tc_ZS2_lincomf1f2[(0,k)],HSN)
        #-----------------------------
        if isogeny_type=="al_lpow":
            print("codomain")
            #codomian 
            lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow=codomain_lv2tnp_al_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
            #-----------------------------------
            #evaluation
            #f_1
            print("eveluation 1/4")
            lv2_f1_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f1,lincom_f1f2[((k+1),0)],lincom_f1f2[(1,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #f_2
            print("eveluation 2/4")
            lv2_f2_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f2,lincom_f1f2[(k,1)],lincom_f1f2[(0,(k+1))],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #ZS1
            print("eveluation 3/4")
            lv2_ZS1_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS1,tc_ZS1_lincomf1f2[(k,0)],tc_ZS1_lincomf1f2[(0,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #ZS2
            print("eveluation 4/4")
            lv2_ZS2_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS2,tc_ZS2_lincomf1f2[(k,0)],tc_ZS2_lincomf1f2[(0,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
        #-----------------------------
        if isogeny_type=="al_sumsq":
            print("codomain")
            #codomian 
            lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow=codomain_lv2tnp_al_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
            #-----------------------------------
            #evaluation
            #f_1
            print("eveluation 1/4")
            lv2_f1_cd=evaluation_lv2_al_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f1,lincom_f1f2[((k+1),0)],lincom_f1f2[(1,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #f_2
            print("eveluation 2/4")
            lv2_f2_cd=evaluation_lv2_al_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f2,lincom_f1f2[(k,1)],lincom_f1f2[(0,(k+1))],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #ZS1
            print("eveluation 3/4")
            lv2_ZS1_cd=evaluation_lv2_al_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS1,tc_ZS1_lincomf1f2[(k,0)],tc_ZS1_lincomf1f2[(0,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #ZS2
            print("eveluation 4/4")
            lv2_ZS2_cd=evaluation_lv2_al_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS2,tc_ZS2_lincomf1f2[(k,0)],tc_ZS2_lincomf1f2[(0,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
   
        #=================================
        #for next step.
        lv2tnp =lv2tnp_cd
        lv2_f1 =lv2_f1_cd
        lv2_f2 =lv2_f2_cd
        lv2_ZS1=lv2_ZS1_cd
        lv2_ZS2=lv2_ZS2_cd
        s*=l
        
        time_end=time.time()
        time_diff=time_end-time_start
        print("Time this step.",time_diff)
        print("")
    #=====================================
    assert(Is_Ellipitc_product(lv2tnp)[0])
    #print("elliptic product",Is_Ellipitc_product(lv2tnp)[0])
    assert(peq_lv2(lv2tnp,lv2_f1))
    assert(peq_lv2(lv2tnp,lv2_f2))
    return lv2tnp,lv2_ZS1,lv2_ZS2



def attack_last(E_0m,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,E_B_L,lv2tnp,lv2_ZS1,lv2_ZS2,zeta_4):
    lv2tnp,lv2_ZS1,lv2_ZS2=Theta_Split(lv2tnp,lv2_ZS1,lv2_ZS2,zeta_4) 
    assert(Is_Ellipitc_product(lv2tnp)[0])
    assert(Is_Ellipitc_product(lv2tnp)[1])
    assert(lv2tnp[0]*lv2tnp[3]==lv2tnp[1]*lv2tnp[2])
    #split to elliptic curves.
    lv2tnp_Ecd_1=[lv2tnp[0],lv2tnp[1]]
    lv2tnp_Ecd_2=[lv2tnp[0],lv2tnp[2]]
    lm_1=lv2tnp_to_Legendre(lv2tnp_Ecd_1)[0]
    lm_2=lv2tnp_to_Legendre(lv2tnp_Ecd_2)[0]
    
    assert(P_A in E_0m)
    assert(Q_A in E_0m)
    assert(order(P_A)==N_A)
    assert(order(Q_A)==N_A)
    lm_0=Elliptic_to_Legendre(E_0m)[0]
    
    if is_isomorphic_Legendre(lm_0,lm_1):
        lm_cd,sqrt_lm_cd,sqrt_lmm1_cd=lv2tnp_to_Legendre(lv2tnp_Ecd_1)
        lv2_P1_cd=[lv2_ZS1[0],lv2_ZS1[1]]
        lv2_P2_cd=[lv2_ZS2[0],lv2_ZS2[1]]
        u_P1_cd=lv2_to_Legendre(lm_cd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_P1_cd)
        u_P2_cd=lv2_to_Legendre(lm_cd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_P2_cd)
    elif is_isomorphic_Legendre(lm_0,lm_2):
        lm_cd,sqrt_lm_cd,sqrt_lmm1_cd=lv2tnp_to_Legendre(lv2tnp_Ecd_2)
        lv2_P1_cd=[lv2_ZS1[0],lv2_ZS1[2]]
        lv2_P2_cd=[lv2_ZS2[0],lv2_ZS2[2]]
        u_P1_cd=lv2_to_Legendre(lm_cd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_P1_cd)
        u_P2_cd=lv2_to_Legendre(lm_cd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_P2_cd)
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


def Main_Attack(E_0m,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_PA,alpha_QA,zeta_4,isogeny_type):
    #prepare
    E_B_L,PA_EB,QA_EB,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2,lv2_ZS2,lv2_ZS2pf1,lv2_ZS2pf2=attack_prepare(E_B,E_pr,N_A,N_B,PA_EB,QA_EB,alpha_PA,alpha_QA)
    #isogeny
    lv2tnp,lv2_ZS1,lv2_ZS2=attack_isogeny(N_A,N_B,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2,lv2_ZS2,lv2_ZS2pf1,lv2_ZS2pf2,isogeny_type)
    #last
    attacker_ker=attack_last(E_0m,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,E_B_L,lv2tnp,lv2_ZS1,lv2_ZS2,zeta_4)
    return attacker_ker





#--------------------------------


def attack_isogeny_for_time(N_A,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2):
    fac=Decomp_degree(N_A)
    #print(fac)
    s=1
    for i in range(0,1):
        time_start=time.time()
        l=fac[i]
        k=N_A//(s*l)
        assert(s*l*k==N_A)
        #print("l=",l)
        #print("elliptic product",Is_Ellipitc_product(lv2tnp)[0])
        HSN=theta_Hadamard(theta_square(lv2tnp))
        #assert(peq_lv2((Mult(lv2tnp,l*k,lv2_f1,HSN)),lv2tnp))
        #assert(peq_lv2((Mult(lv2tnp,l*k,lv2_f2,HSN)),lv2tnp))
        if i!=0:
            assert(not Is_Ellipitc_product(lv2tnp)[0])
            lv2_f12   =Normal_Addition(lv2tnp,lv2_f1,lv2_f2)[0]
            lv2_ZS1pf1=Normal_Addition(lv2tnp,lv2_ZS1,lv2_f1)[0]
            lv2_ZS1pf2=Compatible_Addition(lv2tnp,lv2_ZS1,lv2_f2,lv2_ZS1pf1,lv2_f12)
        #--------------------------------
        #construct kernel.
        lv2_e1 =Mult(lv2tnp,k,lv2_f1 ,HSN)
        lv2_e2 =Mult(lv2tnp,k,lv2_f2 ,HSN)
        lv2_e12=Mult(lv2tnp,k,lv2_f12,HSN)
        #linear combination of f_1,f_2.
        lincom_f1f2={}
        lincom_f1f2[((k+1),0)]=Mult(lv2tnp,k+1,lv2_f1,HSN)
        lincom_f1f2[(0,(k+1))]=Mult(lv2tnp,k+1,lv2_f2,HSN)
        lincom_f1f2[(k,1)]    =kxpy_xpy(lv2tnp,k,lv2_f1,lv2_f2,lv2_f12,HSN)
        lincom_f1f2[(1,k)]    =kxpy_xpy(lv2tnp,k,lv2_f2,lv2_f1,lv2_f12,HSN)
                
        #ZS1+lincom of f_1,f_2.
        tc_ZS1_lincomf1f2={}
        tc_ZS1_lincomf1f2[(0,0)]=lv2_ZS1 #ZS1
        tc_ZS1_lincomf1f2[(k,0)]=kxpy_xpy(lv2tnp,k,lv2_f1,lv2_ZS1,lv2_ZS1pf1,HSN)  #ZS1+e_1
        tc_ZS1_lincomf1f2[(0,k)]=kxpy_xpy(lv2tnp,k,lv2_f2,lv2_ZS1,lv2_ZS1pf2,HSN)  #ZS1+e_2
        #=================================
        isogeny_type="al_lpow"
        if isogeny_type=="al_lpow":
            #codomian 
            lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow=codomain_lv2tnp_al_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
            #-----------------------------------
            #evaluation
            if i==0:
                #f_1
                lv2_f1_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f1,lincom_f1f2[((k+1),0)],lincom_f1f2[(1,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
                #f_2
                lv2_f2_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_f2,lincom_f1f2[(k,1)],lincom_f1f2[(0,(k+1))],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            #ZS1
            lv2_ZS1_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_ZS1,tc_ZS1_lincomf1f2[(k,0)],tc_ZS1_lincomf1f2[(0,k)],mu1_lpow,mu2_lpow,mu12_lpow,HSN)
            time_end=time.time()
        #=================================
        #for next step.
        lv2tnp =lv2tnp_cd
        lv2_f1 =lv2_f1_cd
        lv2_f2 =lv2_f2_cd
        lv2_ZS1=lv2_ZS1_cd
        s*=l
        time_diff=time_end-time_start
    #=====================================
    #print("elliptic product",Is_Ellipitc_product(lv2tnp)[0])
    return lv2tnp,lv2_f1,lv2_f2,lv2_ZS1








def isogeny_time(lv2tnp,lv2_e1,lv2_e2,lv2_x,l,isogeny_type):
    #codomain------------
    time_codomain_start=time.time()
    HSN=theta_Hadamard(theta_square(lv2tnp))
    lv2_e12   =Normal_Addition(lv2tnp,lv2_e1,lv2_e2)[0]
    if isogeny_type=="tl_lpow":
        lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12=codomain_lv2tnp_tl_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    elif isogeny_type=="tl_sumsq":
        lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12=codomain_lv2tnp_tl_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    elif isogeny_type=="al_lpow":
        lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow=codomain_lv2tnp_al_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    elif isogeny_type=="al_sumsq":
        lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow=codomain_lv2tnp_al_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    time_codomain_fin=time.time()
    print("codomain  ",time_codomain_fin-time_codomain_start)
    #---------------------
    
    '''
    #evaluation------------
    time_evaluation_start=time.time()
    lv2_xpe1=Normal_Addition(lv2tnp,lv2_x,lv2_e1)[0]
    lv2_xpe2=Compatible_Addition(lv2tnp,lv2_x,lv2_e2,lv2_xpe1,lv2_e12)
    if isogeny_type=="tl_lpow":
        lv2_ZS1_cd=evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_x,lv2_xpe1,lv2_xpe2,HSN)
    elif isogeny_type=="tl_sumsq":
        print("no code")
    elif isogeny_type=="al_lpow":
        lv2_ZS1_cd=evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_x,lv2_xpe1,lv2_xpe2,mu1_lpow,mu2_lpow,mu12_lpow,HSN)
    elif isogeny_type=="al_sumsq":
        print("not code")
    time_eveluation_fin=time.time()
    print("evaluation",time_eveluation_fin-time_evaluation_start)
    '''
    #-------------------





