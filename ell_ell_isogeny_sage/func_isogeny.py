


def Sum_of_square(l):
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





def Half_coeff_without0(l):
    ld=(l-1)//2
    coeff_set={(1,0),(0,1),(1,1)}
    for k2 in range(2,ld+1):
        coeff_set.add((0,k2))
    for k2 in range(2,l-1):
        coeff_set.add((1,k2))
    for k1 in range(2,ld+1):
        coeff_set.add((k1,0))
    for k1 in range(2,l):
        coeff_set.add((k1,1))
    for k1 in range(2,l-1):
        coeff_set.add((k1,2))
    for k1 in range(2,ld+1):
        for k2 in range(3,l-k1):
            coeff_set.add((k1,k2))
    for k1 in range(ld+1,l):
        for k2 in range(3,l-k1+1):
            coeff_set.add((k1,k2))
    (len(coeff_set)==(l**2-1)//2)
    return coeff_set



#construct k_1*e_1+k_2*e_2.
def Half_LinCom(tc_0: NullCoord,tc_e1:Coord,tc_e2:Coord,tc_e12:Coord):
    l=tc_e1.order
    assert(is_odd(l))
    assert(is_prime(l))
    assert(tc_e2.order==l)
    assert(tc_e12.order==l)
    ld=(l-1)//2
    #lincom[(k1,k2)]=theta coordinate of (k1*e1+k2*e2).
    lincom={(0,0):tc_0,
            (1,0):tc_e1,
            (0,1):tc_e2,
            (1,1):tc_e12}
    for k2 in range(2,ld+1):
        assert(not (0,k2) in lincom)
        lincom[(0,k2)]=tc_0.Diff_Add(lincom[(0,k2-1)],tc_e2,lincom[(0,k2-2)])
    for k2 in range(2,l-1):
        assert(not (1,k2) in lincom)
        lincom[(1,k2)]=tc_0.Diff_Add(lincom[(1,k2-1)],tc_e2,lincom[(1,k2-2)])
    for k1 in range(2,ld+1):
        assert(not (k1,0) in lincom)
        lincom[(k1,0)]=tc_0.Diff_Add(lincom[(k1-1,0)],tc_e1,lincom[(k1-2,0)])
    for k1 in range(2,l):
        assert(not (k1,1) in lincom)
        lincom[(k1,1)]=tc_0.Diff_Add(lincom[(k1-1,1)],tc_e1,lincom[(k1-2,1)])
    for k1 in range(2,l-1):
        assert(not (k1,2) in lincom)
        lincom[(k1,2)]=tc_0.Diff_Add(lincom[(k1-1,2)],tc_e1,lincom[(k1-2,2)])
    for k1 in range(2,ld+1):
        for k2 in range(3,l-k1):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=tc_0.Diff_Add(lincom[(k1,k2-1)],tc_e2,lincom[(k1,k2-2)])
    for k1 in range(ld+1,l):
        for k2 in range(3,l-k1+1):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=tc_0.Diff_Add(lincom[(k1,k2-1)],tc_e2,lincom[(k1,k2-2)])
    del lincom[(0,0)]
    assert(len(lincom.keys())==(l**2-1)//2)
    assert(Half_coeff_without0(l)==lincom.keys())
    return lincom






def Remain_Half_coeff_without0(l):
    remianed_set={(k1,k2) for k1 in range(0,l) for k2 in range(0,l)}
    remianed_set=remianed_set.difference(Half_coeff_without0(l))
    remianed_set=remianed_set.difference({(0,0)})
    (2*len(remianed_set)+1==l**2)
    return remianed_set
    








#construct k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
def LinCom(tc_0:NullCoord,tc_e1:Coord,tc_e2:Coord,tc_e12:Coord):
    l=tc_e1.order
    assert(is_odd(l))
    ld=(l-1)//2
    lincom=Half_LinCom(tc_0,tc_e1,tc_e2,tc_e12)
    lincom[(0,0)]=tc_0
    lmd1_lpow =tc_0.Lmd_lpow(tc_e1)
    lmd2_lpow =tc_0.Lmd_lpow(tc_e2)
    lmd12_lpow=tc_0.Lmd_lpow(tc_e12)
    for k1 in range(ld+1,l):
        lincom[(k1,0)]=lincom[(l-k1,0)].Mult_frac([lmd1_lpow[1]**(2*k1-l),lmd1_lpow[0]**(2*k1-l)])
    for k2 in range(ld+1,l):
        lincom[(0,k2)]=lincom[(0,l-k2)].Mult_frac([lmd2_lpow[1]**(2*k2-l),lmd2_lpow[0]**(2*k2-l)])
    for k1 in range(1,l):
        for k2 in range(ld+1,l):
            if (k1-k2>=0) and (l-k1-k2>=0):
                lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow[1]**(-(k2-k1))*lmd2_lpow[0]**(k1-k2)*lmd12_lpow[0]**(l-k1-k2),lmd1_lpow[0]**(-(k2-k1))*lmd2_lpow[1]**(k1-k2)*lmd12_lpow[1]**(l-k1-k2)])
            elif (k1-k2>=0) and (l-k1-k2<0):
                lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow[1]**(-(k2-k1))*lmd2_lpow[0]**(k1-k2)*lmd12_lpow[1]**(-(l-k1-k2)),lmd1_lpow[0]**(-(k2-k1))*lmd2_lpow[1]**(k1-k2)*lmd12_lpow[0]**(-(l-k1-k2))])
            elif (k1-k2<0) and (l-k1-k2>=0):
                lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow[0]**(k2-k1)*lmd2_lpow[1]**(-(k1-k2))*lmd12_lpow[0]**(l-k1-k2),lmd1_lpow[1]**(k2-k1)*lmd2_lpow[0]**(-(k1-k2))*lmd12_lpow[1]**(l-k1-k2)])
            else:
                lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow[0]**(k2-k1)*lmd2_lpow[1]**(-(k1-k2))*lmd12_lpow[1]**(-(l-k1-k2)),lmd1_lpow[1]**(k2-k1)*lmd2_lpow[0]**(-(k1-k2))*lmd12_lpow[0]**(-(l-k1-k2))])                
    assert(len(lincom)==l**2)
    return lincom







#construct k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
def XpLinCom(tc_0:NullCoord,tc_e1:Coord,tc_e2:Coord,tc_e12:Coord,tc_x:Coord,tc_xpe1:Coord,tc_xpe2:Coord):
    l=tc_e1.order
    assert(is_odd(l))
    tc_xpe12=tc_0.Extended_Addition(tc_x,tc_e1,tc_e2,tc_xpe1,tc_e12,tc_xpe2)
    #tc_xpe12=tc_0.Extended_Addition_2(tc_x,tc_e1,tc_e2,tc_xpe1,tc_e12,tc_xpe2)
    #assert(tc_xpe12_a==tc_xpe12)
    #xplincom[(k1,k2)]=theta coordinate of (x+k1*e1+k2*e2).
    xplincom={}
    xplincom={(0,0):tc_x,
              (1,0):tc_xpe1,
              (0,1):tc_xpe2,
              (1,1):tc_xpe12}
    for k2 in range(2,l):
        xplincom[(0,k2)]=tc_0.Diff_Add(xplincom[(0,k2-1)],tc_e2,xplincom[(0,k2-2)])
        xplincom[(1,k2)]=tc_0.Diff_Add(xplincom[(1,k2-1)],tc_e2,xplincom[(1,k2-2)])
    for k2 in range(0,l):
        for k1 in range(2,l):
            assert(not (k1,k2) in xplincom.keys())
            xplincom[(k1,k2)]=tc_0.Diff_Add(xplincom[(k1-1,k2)],tc_e1,xplincom[(k1-2,k2)])
    assert(len(xplincom.keys())==l**2)
    return xplincom







def Product_power_lambda(tc_0:NullCoord,basis):
    assert(len(basis)==3)
    [tc_e1,tc_e2,tc_e12]=basis
    l=tc_e1.order
    assert(type(tc_e1) ==Coord)
    assert(type(tc_e2) ==Coord)
    assert(type(tc_e12)==Coord)
    assert(tc_e2.order ==l)
    assert(tc_e12.order==l)
    lmd1_lpow =tc_0.Lmd_lpow(tc_e1)
    lmd2_lpow =tc_0.Lmd_lpow(tc_e2)
    lmd12_lpow=tc_0.Lmd_lpow(tc_e12)
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    den_pow         =Power_simultaneously(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,l**2)
    lmd2_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,l**2)
    lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],(l-1)**2)
    lmd_pow_product=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            lmd_pow_product[(k1,k2)]=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow[k1**2+k2**2+k1*k2]]
    lmd_pow_product[(l,0)]=lmd1_lpow_pow[l**2]
    lmd_pow_product[(0,l)]=lmd2_lpow_pow[l**2]
    assert(len(lmd_pow_product)==l**2+2)
    return lmd_pow_product


#============================================================================


def Codomain_common(tc_0:NullCoord,basis:list):
    [tc_e1,tc_e2,tc_e12]=basis
    h_lincom=Half_LinCom(tc_0,tc_e1,tc_e2,tc_e12) 
    l=tc_e1.order
    assert(is_prime(l))
    assert(is_odd(l))
    assert(tc_e2.order==l)
    assert(tc_e12.order==l)
    assert(len(h_lincom)==(l**2-1)//2)
    lmd1_lpow =tc_0.Lmd_lpow(tc_e1)
    lmd2_lpow =tc_0.Lmd_lpow(tc_e2)
    lmd12_lpow=tc_0.Lmd_lpow(tc_e12)
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    return l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den




#take l-th power.
def Codomain_C1(tc_0:NullCoord,basis:list):
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    den_pow         =Power_simultaneously(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
    lmd2_lpow_pow   =Power_simultaneously(lmd2_lpow[0]   ,(l-1)**2)
    lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],(l-1)**2)                                  
    h_lincom_lpow={key:Coord([h_lincom[key].numer[i]**l for i in range(0,4)],h_lincom[key].denom**l) for key in h_lincom.keys()}
    assert(len(h_lincom_lpow)==(l**2-1)//2)
    for (k1,k2) in h_lincom:
        coeff=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow[k1**2+k2**2+k1*k2]]
        h_lincom_lpow[(k1,k2)]=h_lincom_lpow[(k1,k2)].Mult_frac(coeff)
    assert(len(h_lincom_lpow)==(l**2-1)//2)
    h_lincom_lpow[(0,0)]=Coord([tc_0.numer[i]**l for i in range(0,4)],tc_0.denom**l)
    h_lincom_lpow,common_denom_lpow=Dict_common_denom(h_lincom_lpow)
    assert(len(h_lincom_lpow)==(l**2+1)//2)
    tc_f0=NullCoord(h_lincom_lpow[(0,0)].numer,1)
    assert(len(Half_coeff_without0(l))==(l**2-1)//2)
    for key in Half_coeff_without0(l):
        for i in range(0,4):
            tc_f0.numer[i]+=2*h_lincom_lpow[key].numer[i]
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_f0







#l is sum of sq. 
def Codomain_C2(tc_0:NullCoord,basis:list):
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    l,lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    a_u=Sum_of_square(l)
    r=len(a_u)    
    ss1=[(a_u[u]*(l-1))%l         for u in range(0,r)]
    tt1=[(a_u[u]*(l-1)-ss1[u])//l for u in range(0,r)]
    max_exp=(l-1)**2+l*(sum([tt1[u]**2 for u in range(0,r)]))-2*(l-1)*(sum([a_u[u]*tt1[u] for u in range(0,r)]))
    assert(max_exp<=r*(l-1))
    lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,max_exp)
    lmd2_lpow_pow   =Power_simultaneously(lmd2_lpow[0]   ,max_exp)
    lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],max_exp)
    den_pow         =Power_simultaneously(den            ,max_exp*3) 
    #extension of linear combination
    assert(len(lincom)==(l**2-1)//2)
    lincom[(0,0)]=tc_0
    assert(len(lincom)==(l**2+1)//2)
    for (k1,k2) in Remain_Half_coeff_without0(l):
        assert(not ((k1,k2) in lincom))
        a1  =l-2*k1
        a2  =l-2*k2
        adiv=l-k1-k2
        aden=a1+a2+adiv
        assert(aden==3*adiv)
        if k1==0:
            assert(not (0,k2) in lincom)
            assert(a2<=0)
            lincom[(0,k2)]=lincom[(0,l-k2)].Mult_frac([den_pow[-a2],lmd2_lpow_pow[-a2]])
        elif k2==0:
            assert(not (k1,0) in lincom)
            assert(a1<=0)
            lincom[(k1,0)]=lincom[(l-k1,0)].Mult_frac([den_pow[-a1],lmd1_lpow_pow[-a1]])
        elif (a1>=0) and (a2>=0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow_pow[a1]*lmd2_lpow_pow[a2]*den_pow[-aden],lmd_div_lpow_pow[-adiv]])
        elif (a1>=0) and (a2<0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd1_lpow_pow[a1]*den_pow[-aden],lmd2_lpow_pow[-a2]*lmd_div_lpow_pow[-adiv]])
        elif (a1<0) and (a2>=0):
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([lmd2_lpow_pow[a2]*den_pow[-aden],lmd1_lpow_pow[-a1]*lmd_div_lpow_pow[-adiv]])
        else:
            assert((a1<0) and (a2<0))
            assert(not (k1,k2) in lincom)
            lincom[(k1,k2)]=lincom[(l-k1,l-k2)].Mult_frac([den_pow[-aden],lmd1_lpow_pow[-a1]*lmd2_lpow_pow[-a2]*lmd_div_lpow_pow[-adiv]])
    assert(len(lincom)==l**2)
    #assert(LinCom(tc_0,tc_e1,tc_e2,tc_e12)==lincom)
    #-------------------------------------------------
    #calculate theta null point.
    pre_tc_f0=[[prod([tc_0.numer[(a_u[u]%2)*i] for u in range(0,r)]),(tc_0.denom)**r] for i in range(0,4)] 
    for (k1,k2) in Half_coeff_without0(l):
        s1=[(a_u[u]*k1)%l        for u in range(0,r)]
        t1=[(a_u[u]*k1-s1[u])//l for u in range(0,r)]
        s2=[(a_u[u]*k2)%l        for u in range(0,r)]
        t2=[(a_u[u]*k2-s2[u])//l for u in range(0,r)]
        for u in range(0,r):
            assert(a_u[u]*k1==t1[u]*l+s1[u])
            assert(a_u[u]*k2==t2[u]*l+s2[u])
            assert(0<=s1[u]<l)
            assert(0<=s2[u]<l)
        h_1  =k1**2+l*(sum([t1[u]**2    for u in range(0,r)]))-2*k1*(sum([a_u[u]*t1[u] for u in range(0,r)]))
        h_2  =k2**2+l*(sum([t2[u]**2    for u in range(0,r)]))-2*k2*(sum([a_u[u]*t2[u] for u in range(0,r)]))
        h_div=k1*k2+l*(sum([t1[u]*t2[u] for u in range(0,r)]))-k1*sum([a_u[u]*t2[u] for u in range(0,r)])-k2*sum([a_u[u]*t1[u] for u in range(0,r)])
        h_den=h_1+h_2+h_div
        coeff=[lmd1_lpow_pow[h_1]*lmd2_lpow_pow[h_2]*lmd_div_lpow_pow[h_div],den_pow[h_den]]
        for i in range(0,4):
            plus_term_num=coeff[0]*prod([lincom[(s1[u],s2[u])].numer[(a_u[u]%2)*i] for u in range(0,r)])
            plus_term_den=coeff[1]*prod([lincom[(s1[u],s2[u])].denom               for u in range(0,r)])
            pre_tc_f0[i]=Frac_sum(pre_tc_f0[i],[2*plus_term_num,plus_term_den])
    proj_tc_f0=Projective_Theta(pre_tc_f0)
    tc_f0=NullCoord(proj_tc_f0,1)
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_f0
    




#l is sum of sq. 
def Codomain_C3(tc_0:NullCoord,basis:list):
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    ld=(l-1)//2    
    den_pow         =Power_simultaneously(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
    lmd2_lpow_pow   =Power_simultaneously(lmd2_lpow[0]   ,(l-1)**2)
    lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],(ld+1)*ld)
    a_u=Sum_of_square(l)
    r=len(a_u)
    pre_tc_f0=[[prod([tc_0.numer[(a_u[u]%2)*i] for u in range(0,r)]),(tc_0.denom)**r] for i in range(0,4)]
    for (k1,k2) in h_lincom:
        assert((k1,k2)!=(0,0))
        coeff_num=lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2]
        coeff_den=den_pow[k1**2+k2**2+k1*k2]
        ai_k1k2=[(tc_0.Mult(h_lincom[(k1,k2)],a_u[u]))   for u in range(0,r)]
        plus_term_den=coeff_den*prod([(ai_k1k2[u].denom) for u in range(0,r)])
        for i in range(0,4):
            plus_term_num=2*coeff_num*prod([ai_k1k2[u].numer[(a_u[u]%2)*i] for u in range(0,r)])
            pre_tc_f0[i]=Frac_sum(pre_tc_f0[i],[plus_term_num,plus_term_den])
    proj_tc_f0=Projective_Theta(pre_tc_f0)
    tc_f0=NullCoord(proj_tc_f0,tc_0.field(1))
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_f0



#============================================================================


def Evaluation_common(tc_0:Coord,basis:list,pt_list:list,lm_1_lsq,lm_2_lsq):
    assert(len(basis)==3)
    [tc_e1,tc_e2,tc_e12]=basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    l=tc_e1.order
    assert(is_odd(l))
    assert(is_prime(l))
    assert(tc_e2.order==l)
    assert(tc_e12.order==l)
    ml1_lpow =tc_0.MudivLm_lpow(tc_e1,tc_x,tc_xpe1,lm_1_lsq)
    ml2_lpow =tc_0.MudivLm_lpow(tc_e2,tc_x,tc_xpe2,lm_2_lsq)
    [ml1_lpow,ml2_lpow],m_den=Common_denom_frac([ml1_lpow,ml2_lpow])
    assert(ml1_lpow[1]==m_den)
    assert(ml2_lpow[1]==m_den)
    ml1_lpow_pow  =Power_simultaneously(ml1_lpow[0],(l-1))
    ml2_lpow_pow  =Power_simultaneously(ml2_lpow[0],(l-1))
    muden_lpow_pow=Power_simultaneously(m_den      ,2*(l-1))
    return l,ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow




def Evaluation_E1(tc_0:Coord,basis:list,pt_list:list,lmd_pow_product:dict):
    l=basis[0].order
    l,ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Evaluation_common(tc_0,basis,pt_list,lmd_pow_product[(l,0)],lmd_pow_product[(0,l)])
    assert(len(basis)  ==3)
    assert(len(pt_list)==3)
    [tc_e1,tc_e2,tc_e12]  =basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    assert(len(lmd_pow_product)==l**2+2)
    xplincom=XpLinCom(tc_0,tc_e1,tc_e2,tc_e12,tc_x,tc_xpe1,tc_xpe2)
    assert(len(xplincom)==l**2)
    assert(tc_x==xplincom[(0,0)])
    xplincom,den=Dict_common_denom(xplincom)
    assert(len(xplincom)==l**2)
    den_lpow=den**l
    xplincom_lpow={key:Coord([(xplincom[key].numer[i])**l for i in range(0,4)],den_lpow) for key in xplincom.keys()}
    assert(len(xplincom_lpow)==l**2)
    for (k1,k2) in xplincom:
        coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
        xplincom_lpow[(k1,k2)]=xplincom_lpow[(k1,k2)].Mult_frac(coeff)
    xplincom_lpow,den_lpow=Dict_common_denom(xplincom_lpow)
    assert(len(xplincom)==l**2)
    pre_tc_fx=[0,0,0,0]
    for k1 in range(0,l):
        for k2 in range(0,l):
            assert(xplincom_lpow[(k1,k2)].denom==den_lpow)
            for i in range(0,4):
                pre_tc_fx[i]+=xplincom_lpow[(k1,k2)].numer[i]
    tc_fx=Coord(pre_tc_fx,den_lpow)
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_fx




#2
#avoid to take l-th square roots.
def Evaluation_E2(tc_0:NullCoord,basis:list,pt_list:list,lmd_pow_product):
    l=basis[0].order
    l,ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Evaluation_common(tc_0,basis,pt_list,lmd_pow_product[(l,0)],lmd_pow_product[(0,l)])
    [tc_e1,tc_e2,tc_e12]=basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    a_u=Sum_of_square(l)
    r=len(a_u)    
    xplincom=[]
    for u in range(0,r):
        au=a_u[u]
        au_e1  =tc_0.Mult(tc_e1  ,au)
        au_e2  =tc_0.Mult(tc_e2  ,au)
        au_e12 =tc_0.Mult(tc_e12 ,au)
        au_e1.order =l
        au_e2.order =l
        au_e12.order=l
        au_x   =tc_0.Mult(tc_x   ,au)
        au_xpe1=tc_0.Mult(tc_xpe1,au)
        au_xpe2=tc_0.Mult(tc_xpe2,au)
        aux_plus_lincom=XpLinCom(tc_0,au_e1,au_e2,au_e12,au_x,au_xpe1,au_xpe2)
        assert(len(aux_plus_lincom.keys())==l**2)
        xplincom.append(aux_plus_lincom)
    pre_tc_fx=[[0,1],[0,1],[0,1],[0,1]]
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
            for i in range(0,4):
                plus_term_num=coeff[0]*prod([xplincom[u][(k1,k2)].numer[(a_u[u]%2)*i] for u in range(0,r)])
                plus_term_den=coeff[1]*prod([xplincom[u][(k1,k2)].denom               for u in range(0,r)])
                pre_tc_fx[i]=Frac_sum(pre_tc_fx[i],[plus_term_num,plus_term_den])
    pre_tc_fx,den=Common_denom_frac(pre_tc_fx)
    tc_fx=Coord([pre_tc_fx[i][0] for i in range(0,4)],den)
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_fx
    




def Evaluation_E3(tc_0:NullCoord,basis:list,pt_list:list,lmd_pow_product):
    l=basis[0].order
    l,ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Evaluation_common(tc_0,basis,pt_list,lmd_pow_product[(l,0)],lmd_pow_product[(0,l)])
    [tc_e1,tc_e2,tc_e12]=basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    #ct_mlt=(tc_0.field).n_mul
    #ct_sqr=(tc_0.field).n_sqr
    a_u=Sum_of_square(l)
    r=len(a_u)
    xplincom=XpLinCom(tc_0,tc_e1,tc_e2,tc_e12,tc_x,tc_xpe1,tc_xpe2)
    assert(len(xplincom)==l**2)
    assert(xplincom[(0,0)]==tc_x)
    pre_tc_fx=[[0,1],[0,1],[0,1],[0,1]]
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
            au_k1k2=[tc_0.Mult(xplincom[(k1,k2)],a_u[u]) for u in range(0,r)]
            for i in range(0,4):
                plus_term_num=coeff[0]*prod([au_k1k2[u].numer[(a_u[u]%2)*i] for u in range(0,r)])
                plus_term_den=coeff[1]*prod([au_k1k2[u].denom               for u in range(0,r)])
                pre_tc_fx[i]=Frac_sum(pre_tc_fx[i],[plus_term_num,plus_term_den])
    pre_tc_fx,den=Common_denom_frac(pre_tc_fx)
    tc_fx=Coord([pre_tc_fx[i][0] for i in range(0,4)],den)
    #print("M",(tc_0.field).n_mul-ct_mlt)
    #print("S",(tc_0.field).n_sqr-ct_sqr)
    #print("M+S",(tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
    return tc_fx
    


