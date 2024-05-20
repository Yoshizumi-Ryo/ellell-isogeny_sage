

#----------------------------------------

def sum_of_square(l):
    assert(is_prime(l))
    assert(l!=2)
    if l%4==1:
        for a in range(0,l):
            if is_square(l-a**2):
                b=sqrt(l-a**2)
                return (a,b)
    else:
        for a in range(0,l):
            for b in range(a,l):
                for c in range(b,l):
                    if is_square(l-a**2-b**2-c**2):
                        d=sqrt(l-a**2-b**2-c**2)
                        if a==0:
                            assert(b!=0)
                            return (b,c,d)
                        else:
                            return (a,b,c,d)
                    
                    

#codomain=====================================

#for normalization l-torsion pts.
def mu_ell_power(lv2tnp,lv2_e,l,HSN):
    half_lm1=(l-1)//2
    for i_0 in range(0,4):
       denom=Mult(lv2tnp,half_lm1+1,lv2_e,HSN)[i_0]
       if denom!=0:
          mu_l_pow=Mult(lv2tnp,half_lm1,lv2_e,HSN)[i_0]/denom
          return mu_l_pow
       

def Is_excellent_kernel(lv2tnp,l,lv2_e,HSN):
    return Mult(lv2tnp,l,lv2_e,HSN)==lv2tnp


#construct k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
def linear_combination_ell_torsion(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    #lincom[(k1,k2)]=theta coordinate of (k1*e1+k2*e2).
    lincom={(0,0):lv2tnp,
            (1,0):lv2_e1,
            (0,1):lv2_e2,
            (1,1):lv2_e12}
    for k2 in range(2,l):
        lincom[(0,k2)]=DiffAdd(lincom[(0,k2-1)],lv2_e2,lincom[(0,k2-2)],HSN)
        lincom[(1,k2)]=DiffAdd(lincom[(1,k2-1)],lv2_e2,lincom[(1,k2-2)],HSN)
    for k2 in range(0,l):
        for k1 in range(2,l):
            lincom[(k1,k2)]=DiffAdd(lincom[(k1-1,k2)],lv2_e1,lincom[(k1-2,k2)],HSN)
    return lincom






#construct k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
def linear_combination_ell_torsion_2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    #lincom[(k1,k2)]=theta coordinate of (k1*e1+k2*e2).
    ld=(l-1)//2
    lincom={(0,0):lv2tnp,
            (1,0):lv2_e1,
            (0,1):lv2_e2,
            (1,1):lv2_e12}
    for k2 in range(2,ld+1):
        lincom[(0,k2)]=DiffAdd(lincom[(0,k2-1)],lv2_e2,lincom[(0,k2-2)],HSN)
        lincom[(1,k2)]=DiffAdd(lincom[(1,k2-1)],lv2_e2,lincom[(1,k2-2)],HSN)
    for k2 in range(0,ld+1):
        for k1 in range(2,l):
            lincom[(k1,k2)]=DiffAdd(lincom[(k1-1,k2)],lv2_e1,lincom[(k1-2,k2)],HSN)
    return lincom


#------------------------


'''
#1
#avoid to take l-th square roots.
def codomain_lv2tnp_al_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    assert(is_prime(l))
    assert(l!=2)
    lincom=linear_combination_ell_torsion(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    mu1_lpow =mu_ell_power(lv2tnp,lv2_e1 ,l,HSN)
    mu2_lpow =mu_ell_power(lv2tnp,lv2_e2 ,l,HSN)
    mu12_lpow=mu_ell_power(lv2tnp,lv2_e12,l,HSN)
    mmm=mu12_lpow/(mu1_lpow*mu2_lpow)
    mu1_lpow_pow =[mu1_lpow**i  for i in range(0,l**2)]
    mu2_lpow_pow =[mu2_lpow**i  for i in range(0,l**2)]
    mmm_lpow_pow =[mmm**i       for i in range(0,l**2)]
    lv2tnp_cd=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            kk1=min(k1,l-k1)
            kk2=min(k2,l-k2)
            print(kk1,kk2)
            coff=mu1_lpow_pow[kk1**2]*mu2_lpow_pow[kk2**2]*mmm_lpow_pow[(kk1*kk2)]
            for i in range(0,4):
                lv2tnp_cd[i]+=coff*(lincom[(kk1,kk2)][i]**l)
    return lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow
'''

#C1
#avoid to take l-th square roots.
def codomain_lv2tnp_al_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    assert(is_prime(l))
    assert(l!=2)
    ld=(l-1)//2
    time_1=time.time()
    lincom=linear_combination_ell_torsion_2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    time_2=time.time()
    mu1_lpow =mu_ell_power(lv2tnp,lv2_e1 ,l,HSN)
    mu2_lpow =mu_ell_power(lv2tnp,lv2_e2 ,l,HSN)
    mu12_lpow=mu_ell_power(lv2tnp,lv2_e12,l,HSN)
    mmm=mu12_lpow/(mu1_lpow*mu2_lpow)
    time_3=time.time()
    mu1_lpow_pow=[mu1_lpow**i  for i in range(0,(l-1)**2+1)]
    mu2_lpow_pow=[mu2_lpow**i  for i in range(0,ld**2+1)]
    mmm_pow     =[mmm**i       for i in range(0,((l-1)*ld)+1)]
    lv2tnp_cd=[lv2tnp[i]**l for i in range(0,4)]
    time_4=time.time()
    coff_set={(k1,0) for k1 in range(1,ld+1)}|{(k1,k2) for k1 in range(0,l) for k2 in range(1,ld+1)}
    for kk in coff_set:
        coff=mu1_lpow_pow[kk[0]**2]*mu2_lpow_pow[kk[1]**2]*mmm_pow[kk[0]*kk[1]]
        for i in range(0,4):
            lv2tnp_cd[i]+=2*coff*(lincom[(kk[0],kk[1])][i]**l)
    time_5=time.time()
    #print(100,time_2-time_1)
    #print(101,time_3-time_2)
    #print(102,time_4-time_3)
    #print(103,time_5-time_4)
    return lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow




#-----------------------

#C3
def codomain_lv2tnp_al_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    assert(is_prime(l))
    assert(l!=2)
    a_i=sum_of_square(l)
    r=len(a_i)
    time_1=time.time()
    lincom=linear_combination_ell_torsion(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    time_2=time.time()
    mu1_lpow =mu_ell_power(lv2tnp,lv2_e1 ,l,HSN)
    mu2_lpow =mu_ell_power(lv2tnp,lv2_e2 ,l,HSN)
    mu12_lpow=mu_ell_power(lv2tnp,lv2_e12,l,HSN)
    mmm=mu12_lpow/(mu1_lpow*mu2_lpow)
    time_3=time.time()
    mu1_lpow_pow=[mu1_lpow**i  for i in range(0,(l-1)**2+1)]
    mu2_lpow_pow=[mu2_lpow**i  for i in range(0,(l-1)**2+1)]
    mmm_pow     =[mmm**i       for i in range(0,(l-1)**2+1)]
    time_4=time.time()
    lv2tnp_cd=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            coff=mu1_lpow_pow[k1**2]*mu2_lpow_pow[k2**2]*mmm_pow[k1*k2]
            ai_P=[Mult(lv2tnp,a_i[i],lincom[(k1,k2)],HSN) for i in range(0,r)]
            for j in range(0,4):
                value=1
                for i in range(0,r):
                    if is_odd(a_i[i]):
                        index=j
                    else:
                        index=0
                    value*=ai_P[i][index]
                lv2tnp_cd[j]+=value*coff
    time_5=time.time()
    #print(100,time_2-time_1)
    #print(101,time_3-time_2)
    #print(102,time_4-time_3)
    #print(103,time_5-time_4)
    return lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow




#take l-th square roots.
def codomain_lv2tnp_tl_lpow(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    assert(is_prime(l))
    assert(l!=2)
    mu1_lpow =mu_ell_power(lv2tnp,lv2_e1 ,l,HSN)
    mu2_lpow =mu_ell_power(lv2tnp,lv2_e2 ,l,HSN)
    mu12_lpow=mu_ell_power(lv2tnp,lv2_e12,l,HSN)
    K=parent(mu1_lpow)
    X=gen(K['X'])
    poly1 =X**l-mu1_lpow
    poly2 =X**l-mu2_lpow
    poly12=X**l-mu12_lpow
    mu1 =poly1.roots()[1][0]
    mu2 =poly2.roots()[1][0]
    mu12=poly12.roots()[1][0]
    assert(mu1**l ==mu1_lpow)
    assert(mu2**l ==mu2_lpow)
    assert(mu12**l==mu12_lpow)
    lv2_e1 =Mult_Const(mu1 ,lv2_e1)
    lv2_e2 =Mult_Const(mu2 ,lv2_e2)
    lv2_e12=Mult_Const(mu12,lv2_e12)
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e1 ,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e2 ,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e12,HSN))
    lincom=linear_combination_ell_torsion(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    lv2tnp_cd=[0,0,0,0]
    for i in range(0,4):
        for k1 in range(0,l):
            for k2 in range(0,l):
                lv2tnp_cd[i]+=lincom[(k1,k2)][i]**l
    return lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12







#take l-th square roots.
def codomain_lv2tnp_tl_sumsq(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN):
    assert(is_prime(l))
    assert(l!=2)
    sum_sq_l=sum_of_square(l)
    r=len(sum_sq_l)
    mu1_lpow =mu_ell_power(lv2tnp,lv2_e1 ,l,HSN)
    mu2_lpow =mu_ell_power(lv2tnp,lv2_e2 ,l,HSN)
    mu12_lpow=mu_ell_power(lv2tnp,lv2_e12,l,HSN)
    K=parent(mu1_lpow)
    X=gen(K['X'])
    poly1 =X**l-mu1_lpow
    poly2 =X**l-mu2_lpow
    poly12=X**l-mu12_lpow
    mu1 =poly1.roots()[1][0]
    mu2 =poly2.roots()[1][0]
    mu12=poly12.roots()[1][0]
    assert(mu1**l ==mu1_lpow)
    assert(mu2**l ==mu2_lpow)
    assert(mu12**l==mu12_lpow)
    lv2_e1 =Mult_Const(mu1 ,lv2_e1)
    lv2_e2 =Mult_Const(mu2 ,lv2_e2)
    lv2_e12=Mult_Const(mu12,lv2_e12)
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e1 ,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e2 ,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e12,HSN))
    lincom=linear_combination_ell_torsion(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,HSN)
    lv2tnp_cd=[0,0,0,0]
    for j in range(0,4):
        for k1 in range(0,l):
            for k2 in range(0,l):
                value=1
                for i in range(0,r):
                    a_i=sum_sq_l[i]
                    ai_P=Mult(lv2tnp,a_i,lincom[(k1,k2)],HSN)
                    if is_odd(a_i):
                        index=j
                    else:
                        index=0
                    value*=ai_P[index]
                lv2tnp_cd[j]+=value
    return lv2tnp_cd,lv2_e1,lv2_e2,lv2_e12



#evaluation=====================



def nu_ell_power_al(lv2tnp,lv2_e,lv2_x,lv2_xpe,l,mu_lpow,HSN):
    #x+ell*e.
    #assert(Is_given_order(lv2tnp,l,lv2_e,HSN))
    lv2_xple=kxpy_xpy(lv2tnp,l,lv2_e,lv2_x,lv2_xpe,HSN)
    assert(peq_lv2(lv2_x,lv2_xple))
    nu_l_pow=lv2_x[0]/((mu_lpow**(l-1))*lv2_xple[0])
    return nu_l_pow



def nu_ell_power_tl(lv2tnp,lv2_e,lv2_x,lv2_xpe,l,HSN):
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e,HSN))
    #x+ell*e.
    lv2_xple=kxpy_xpy(lv2tnp,l,lv2_e,lv2_x,lv2_xpe,HSN)
    assert(peq_lv2(lv2_x,lv2_xple))
    nu_l_pow=lv2_x[0]/(lv2_xple[0])
    return nu_l_pow



def x_plus_linear_combination(lv2tnp,l,lv2_e1,lv2_e2,lv2_x,lv2_xpe1,lv2_xpe2,lv2_xpe12,HSN):
    xplincom={(0,0):lv2_x,
              (1,0):lv2_xpe1,
              (0,1):lv2_xpe2,
              (1,1):lv2_xpe12}
    for k2 in range(2,l):
        xplincom[(0,k2)]=DiffAdd(xplincom[(0,k2-1)],lv2_e2,xplincom[(0,k2-2)],HSN)
        xplincom[(1,k2)]=DiffAdd(xplincom[(1,k2-1)],lv2_e2,xplincom[(1,k2-2)],HSN)
    for k2 in range(0,l):
        for k1 in range(2,l):
            xplincom[(k1,k2)]=DiffAdd(xplincom[(k1-1,k2)],lv2_e1,xplincom[(k1-2,k2)],HSN)
    return xplincom


#E1
def evaluation_lv2_al_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_x,lv2_xpe1,lv2_xpe2,mu1_lpow,mu2_lpow,mu12_lpow,HSN):
    lv2_xpe12=Extended_Addition(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_e12,lv2_xpe2)
    nu1_lpow=nu_ell_power_al(lv2tnp,lv2_e1,lv2_x,lv2_xpe1,l,mu1_lpow,HSN)
    nu2_lpow=nu_ell_power_al(lv2tnp,lv2_e2,lv2_x,lv2_xpe2,l,mu2_lpow,HSN)
    xpKer=x_plus_linear_combination(lv2tnp,l,lv2_e1,lv2_e2,lv2_x,lv2_xpe1,lv2_xpe2,lv2_xpe12,HSN)
    lv2_fx=[0,0,0,0]
    for i in range(0,4):
        for k1 in range(0,l):
            for k2 in range(0,l):
                normalize_coff=(mu1_lpow**(k1*(k1-k2-1)))*(mu2_lpow**(k2*(k2-k1-1)))*(mu12_lpow**(k1*k2))*(nu1_lpow**k1)*(nu2_lpow**k2)
                lv2_fx[i]+=normalize_coff*((xpKer[(k1,k2)][i])**l)
    return lv2_fx








def evaluation_lv2_tl_lpow(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_x,lv2_xpe1,lv2_xpe2,HSN):
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e1,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e2,HSN))
    nu1_lpow=nu_ell_power_tl(lv2tnp,lv2_e1,lv2_x,lv2_xpe1,l,HSN)
    nu2_lpow=nu_ell_power_tl(lv2tnp,lv2_e2,lv2_x,lv2_xpe2,l,HSN)
    K=parent(nu1_lpow)
    X=gen(K['X'])
    poly1 =X**l-nu1_lpow
    poly2 =X**l-nu2_lpow
    nu1 =poly1.roots()[1][0]
    nu2 =poly2.roots()[1][0]
    assert(nu1**l==nu1_lpow)
    assert(nu2**l==nu2_lpow)
    lv2_xpe1=Mult_Const(nu1,lv2_xpe1)
    lv2_xpe2=Mult_Const(nu2,lv2_xpe2)
    lv2_xpe12=Extended_Addition(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_e12,lv2_xpe2)
    #lv2_xpe12=Normalize_for_RR_Threeway(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_xpe2,lv2_e12,lv2_xpe12)
    xpKer=x_plus_linear_combination(lv2tnp,l,lv2_e1,lv2_e2,lv2_x,lv2_xpe1,lv2_xpe2,lv2_xpe12,HSN)
    lv2_fx=[0,0,0,0]
    for i in range(0,4):
        for k1 in range(0,l):
            for k2 in range(0,l):
                lv2_fx[i]+=(xpKer[(k1,k2)][i])**l
    return lv2_fx





def evaluation_lv2_tl_sumsq(lv2tnp,l,lv2_e1,lv2_e2,lv2_e12,lv2_x,lv2_xpe1,lv2_xpe2,HSN):
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e1,HSN))
    assert(Is_excellent_kernel(lv2tnp,l,lv2_e2,HSN))
    sum_sq_l=sum_of_square(l)
    r=len(sum_sq_l)
    nu1_lpow=nu_ell_power_tl(lv2tnp,lv2_e1,lv2_x,lv2_xpe1,l,HSN)
    nu2_lpow=nu_ell_power_tl(lv2tnp,lv2_e2,lv2_x,lv2_xpe2,l,HSN)
    K=parent(nu1_lpow)
    X=gen(K['X'])
    poly1=X**l-nu1_lpow
    poly2=X**l-nu2_lpow
    nu1=poly1.roots()[1][0]
    nu2=poly2.roots()[1][0]
    assert(nu1**l==nu1_lpow)
    assert(nu2**l==nu2_lpow)
    lv2_xpe1=Mult_Const(nu1,lv2_xpe1)
    lv2_xpe2=Mult_Const(nu2,lv2_xpe2)
    lv2_xpe12=Extended_Addition(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_e12,lv2_xpe2)
    #lv2_xpe12=Normalize_for_RR_Threeway(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_xpe2,lv2_e12,lv2_xpe12)
    xpKer=x_plus_linear_combination(lv2tnp,l,lv2_e1,lv2_e2,lv2_x,lv2_xpe1,lv2_xpe2,lv2_xpe12,HSN)
    lv2_fx=[0,0,0,0]
    for j in range(0,4):
        for k1 in range(0,l):
            for k2 in range(0,l):
                value=1
                for i in range(0,r):
                    a_i=sum_sq_l[i]
                    ai_P=Mult(lv2tnp,a_i,xpKer[(k1,k2)],HSN)
                    if is_odd(a_i):
                        index=j
                    else:
                        index=0
                    value*=ai_P[index]
                lv2_fx[j]+=value
    return lv2_fx



