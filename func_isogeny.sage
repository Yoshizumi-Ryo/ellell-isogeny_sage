


def Sum_of_square(l:int):
    assert(is_prime(l))
    assert(l!=2)
    if l==3:
        return (1,1,1)
    if l%4==1:
        for a in range(0,l):
            if is_square(l-a**2):
                b=sqrt(l-a**2)
                return (a,b)
    if l%8==3:
        for a in range(1,l):
            if is_square(l-2*a**2):
                b=sqrt(l-2*a**2)
                return (a,a,b)
    if l%8==7 and l%3==1:
        for a in range(1,l):
            if is_square(l-3*a**2):
                b=sqrt(l-3*a**2)
                return (a,a,a,b)
      
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




def Sum_of_square_mod24(l:int):
    assert(is_prime(l))
    assert(l!=2)
    if l==3:
        return (1,1,1)
    if l%4==1:
        for a in range(0,l):
            if is_square(l-a**2):
                b=sqrt(l-a**2)
                return (a,b)
    if l%8==3:
        for a in range(1,l):
            if is_square(l-2*a**2):
                b=sqrt(l-2*a**2)
                return (a,a,b)
    if l%8==7 and l%3==1:
        for a in range(1,l):
            if is_square(l-3*a**2):
                b=sqrt(l-3*a**2)
                return (a,a,a,b)
    if l==23:
        return (1,1,1,1,1,3,3)
    if l==47:
        return (1,1,3,3,3,3,3)
    if l==71:
        return (2,3,3,7)
    if l==167:
        return (1,1,1,1,1,9,9)
    if l==191:
        return (4,5,5,5,5,5,5,5)               
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




def All_Sum_of_square(l:int):
    ct=0
    result=[]
    assert(is_prime(l))
    assert(l!=2)
    if l%4==1:
        for a in range(0,l):
            if is_square(l-a**2):
                b=sqrt(l-a**2)
                if a<=b:
                    result.append((a,b))
                    ct+=1
    else:
        for a in range(0,l):
            for b in range(a,l):
                for c in range(b,l):
                    if is_square(l-a**2-b**2-c**2):
                        d=sqrt(l-a**2-b**2-c**2)
                        if c<=d:
                            if a==0:
                                result.append((b,c,d))
                            else:
                                result.append((a,b,c,d))
    return result
                          




def Half_coeff_without0(l:int):
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
def Half_LinCom(tc_0:NullCoord,tc_e1:Coord,tc_e2:Coord,tc_e12:Coord,l:int):
    assert(is_odd(l))
    assert(is_prime(l))
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






def Remain_Half_coeff_without0(l:int):
    remianed_set={(k1,k2) for k1 in range(0,l) for k2 in range(0,l)}
    remianed_set=remianed_set.difference(Half_coeff_without0(l))
    remianed_set=remianed_set.difference({(0,0)})
    assert(2*len(remianed_set)+1==l**2)
    return remianed_set
    



#construct k_1*e_1+k_2*e_2 for 0<=k_1,k_2<l.
def XpLinCom(tc_0:NullCoord,basis:list,pt_list:list):
    [tc_e1,tc_e2,tc_e12]=basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    l=tc_e1.order
    assert(is_odd(l))
    tc_xpe12=tc_0.Extended_Addition(tc_x,tc_e1,tc_e2,tc_xpe1,tc_e12,tc_xpe2)
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
    for key in xplincom:
        type(xplincom[key]==Coord)
    return xplincom



#============================================================================


def Lmd_lpow(tc_0: NullCoord,tc_e:Coord,tc_lde,tc_ldm1e):
    try:
        return tc_e.lmd_lpow_value
    except AttributeError:
        assert(tc_lde!=0)
        assert(tc_ldm1e!=0)
        tc_ldp1e=tc_0.Diff_Add(tc_lde,tc_e,tc_ldm1e)
        tc_e.lmd_lpow_value=[tc_lde.numer[0]*tc_ldp1e.denom,tc_ldp1e.numer[0]*tc_lde.denom]
        return tc_e.lmd_lpow_value



def Codomain_common(tc_0:NullCoord,basis:list):
    [tc_e1,tc_e2,tc_e12]=basis
    l=tc_e1.order
    h_lincom=Half_LinCom(tc_0,tc_e1,tc_e2,tc_e12,l) 
    assert(is_prime(l))
    assert(len(h_lincom)==(l**2-1)//2)
    ld=(l-1)//2
    h_lincom[(0,0)]=tc_0
    lmd1_lpow =Lmd_lpow(tc_0,tc_e1 ,h_lincom[(ld,0)] ,h_lincom[(ld-1,0)]   )
    lmd2_lpow =Lmd_lpow(tc_0,tc_e2 ,h_lincom[(0,ld)] ,h_lincom[(0,ld-1)]   )
    lmd12_lpow=Lmd_lpow(tc_0,tc_e12,h_lincom[(ld,ld)],h_lincom[(ld-1,ld-1)])
    del h_lincom[(0,0)]
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    return l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den





"""
def Codomain_C1(tc_0:NullCoord,basis:list):
    [tc_e1,tc_e2,tc_e12]=basis
    l=tc_e1.order
    ld=(l-1)//2
    assert(is_prime(l))
    a_u=Sum_of_square_mod24(l)
    r=len(a_u)
    max_au=max(set(a_u))
    mult_e1  =tc_0.Mult_Mult(tc_e1  ,max_au)
    mult_e2  =tc_0.Mult_Mult(tc_e2  ,max_au)
    mult_e12 =tc_0.Mult_Mult(tc_e12 ,max_au)
    lincom_au=[Half_LinCom(tc_0,mult_e1[a_u[0]],mult_e2[a_u[0]],mult_e12[a_u[0]],l)]
    for u in range(1,r):
        if a_u[u]==a_u[u-1]:
            lincom_au.append(lincom_au[u-1])
        else:
            lincom_au.append(Half_LinCom(tc_0,mult_e1[a_u[u]],mult_e2[a_u[u]],mult_e12[a_u[u]],l))
    lmd1_lpow =Lmd_lpow(tc_0,tc_e1 ,tc_0.Mult(tc_e1,ld) ,tc_0.Mult(tc_e1,ld-1))
    lmd2_lpow =Lmd_lpow(tc_0,tc_e2 ,tc_0.Mult(tc_e2,ld) ,tc_0.Mult(tc_e2,ld-1))
    lmd12_lpow=Lmd_lpow(tc_0,tc_e12,tc_0.Mult(tc_e12,ld),tc_0.Mult(tc_e12,ld-1))
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    den_pow         =Multpower_straight(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_sq(lmd1_lpow[0]   ,l)
    lmd2_lpow_pow   =Multpower_sq(lmd2_lpow[0]   ,l)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(ld+1)*ld)
    r=len(a_u)
    pre_tc_f0=[[prod([tc_0.numer[(a_u[u]%2)*i] for u in range(0,r)]),(tc_0.denom)**r] for i in range(0,4)]
    for (k1,k2) in Half_coeff_without0(l):
        assert((k1,k2)!=(0,0))
        auk1e1_auk2e2=[lincom_au[u][(k1,k2)] for u in range(0,r)]
        coeff_num=lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2]
        coeff_den=den_pow[k1**2+k2**2+k1*k2]
        plus_term_den=coeff_den*prod([(auk1e1_auk2e2[u].denom) for u in range(0,r)])
        for i in range(0,4):
            plus_term_num=2*coeff_num*prod([auk1e1_auk2e2[u].numer[(a_u[u]%2)*i] for u in range(0,r)])
            pre_tc_f0[i]=Frac_add(pre_tc_f0[i],[plus_term_num,plus_term_den])
    proj_tc_f0=Projective_Theta(pre_tc_f0)
    tc_f0=NullCoord(proj_tc_f0,tc_0.field(1))
    return tc_f0
"""




#l is sum of sq. 
def CodSq(tc_0:NullCoord,basis:list):
    l,lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    a_u=Sum_of_square(l)
    r=len(a_u)    
    ss1=[(a_u[u]*(l-1))%l         for u in range(0,r)]
    tt1=[(a_u[u]*(l-1)-ss1[u])//l for u in range(0,r)]
    max_exp=(l-1)**2+l*(sum([tt1[u]**2 for u in range(0,r)]))-2*(l-1)*(sum([a_u[u]*tt1[u] for u in range(0,r)]))
    assert(max_exp<=r*(l-1))
    lmd1_lpow_pow   =Multpower_straight(lmd1_lpow[0]   ,max_exp)
    lmd2_lpow_pow   =Multpower_straight(lmd2_lpow[0]   ,max_exp)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],max_exp)
    den_pow         =Multpower_straight(den            ,max_exp*3) 
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
            pre_tc_f0[i]=Frac_add(pre_tc_f0[i],[2*plus_term_num,plus_term_den])
    proj_tc_f0=Projective_Theta(pre_tc_f0)
    tc_f0=NullCoord(proj_tc_f0,1)
    return tc_f0
    





#take l-th power.
def CodOne(tc_0:NullCoord,basis:list):
    l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    den_pow         =Multpower_straight(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_sq(lmd1_lpow[0],l)
    lmd2_lpow_pow   =Multpower_sq(lmd2_lpow[0],l)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(l-1)**2)                                  
    h_lincom_lpow={key:Coord([h_lincom[key].numer[i]**l for i in range(0,4)],h_lincom[key].denom**l) for key in h_lincom.keys()}
    assert(len(h_lincom_lpow)==(l**2-1)//2)
    for (k1,k2) in h_lincom:
        coeff=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow[k1**2+k2**2+k1*k2]]
        h_lincom_lpow[(k1,k2)]=h_lincom_lpow[(k1,k2)].Mult_frac(coeff)
    assert(len(h_lincom_lpow)==(l**2-1)//2)
    h_lincom_lpow[(0,0)]=Coord([tc_0.numer[i]**l for i in range(0,4)],tc_0.denom**l)
    h_lincom_lpow,_=Dict_common_denom(h_lincom_lpow)
    assert(len(h_lincom_lpow)==(l**2+1)//2)
    tc_f0=NullCoord(h_lincom_lpow[(0,0)].numer,1)
    assert(len(Half_coeff_without0(l))==(l**2-1)//2)
    for key in Half_coeff_without0(l):
        for i in range(0,4):
            tc_f0.numer[i]+=2*h_lincom_lpow[key].numer[i]
    return tc_f0




"""
#l is sum of sq. 
def Codomain_C1(tc_0:NullCoord,basis:list):
    l,h_lincom,lmd1_lpow,lmd2_lpow,lmd_div_lpow,den=Codomain_common(tc_0,basis)
    ld=(l-1)//2    
    den_pow         =Multpower_straight(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_sq(lmd1_lpow[0]   ,l)
    lmd2_lpow_pow   =Multpower_sq(lmd2_lpow[0]   ,l)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(ld+1)*ld)
    a_u=Sum_of_square(l)
    max_au=max(set(a_u))
    r=len(a_u)
    pre_tc_f0=[[prod([tc_0.numer[(a_u[u]%2)*i] for u in range(0,r)]),(tc_0.denom)**r] for i in range(0,4)]
    for (k1,k2) in h_lincom:
        assert((k1,k2)!=(0,0))
        coeff_num=lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2]
        coeff_den=den_pow[k1**2+k2**2+k1*k2]
        mult_k1e1k2e2=tc_0.Mult_Mult(h_lincom[(k1,k2)],max_au)
        au_k1k2=[mult_k1e1k2e2[a_u[u]]   for u in range(0,r)]
        plus_term_den=coeff_den*prod([(au_k1k2[u].denom) for u in range(0,r)])
        for i in range(0,4):
            plus_term_num=2*coeff_num*prod([au_k1k2[u].numer[(a_u[u]%2)*i] for u in range(0,r)])
            pre_tc_f0[i]=Frac_add(pre_tc_f0[i],[plus_term_num,plus_term_den])
    proj_tc_f0=Projective_Theta(pre_tc_f0)
    tc_f0=NullCoord(proj_tc_f0,tc_0.field(1))
    return tc_f0
"""


#============================================================================




def Product_power_lambda(basis:list):
    assert(len(basis)==3)
    [tc_e1,tc_e2,tc_e12]=basis
    l=tc_e1.order
    assert(is_prime(l))
    lmd1_lpow =tc_e1.lmd_lpow_value
    lmd2_lpow =tc_e2.lmd_lpow_value
    lmd12_lpow=tc_e12.lmd_lpow_value
    lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
    [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
    den_pow         =Multpower_straight(den            ,3*(l-1)**2)
    lmd1_lpow_pow   =Multpower_straight(lmd1_lpow[0]   ,l**2)
    lmd2_lpow_pow   =Multpower_straight(lmd2_lpow[0]   ,l**2)
    lmd_div_lpow_pow=Multpower_straight(lmd_div_lpow[0],(l-1)**2)
    lmd_pow_product=dict()
    for k1 in range(0,l):
        for k2 in range(0,l):
            lmd_pow_product[(k1,k2)]=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow[k1**2+k2**2+k1*k2]]
    lm_1_lsq=[lmd1_lpow_pow[l],den_pow[l]]
    lm_2_lsq=[lmd2_lpow_pow[l],den_pow[l]]
    assert(len(lmd_pow_product)==l**2)
    return [lmd_pow_product,lm_1_lsq,lm_2_lsq]






def Take_lmd_power(tc_x:Coord,xple1:Coord,xple2:Coord,lm_1_lsq:Coord,lm_2_lsq:Coord,l:int):
    ml1_lpow=[tc_x.numer[0]*(lm_1_lsq[1])*xple1.denom,tc_x.denom*(lm_1_lsq[0])*xple1.numer[0]]
    ml2_lpow=[tc_x.numer[0]*(lm_2_lsq[1])*xple2.denom,tc_x.denom*(lm_2_lsq[0])*xple2.numer[0]]
    [ml1_lpow,ml2_lpow],m_den=Common_denom_frac([ml1_lpow,ml2_lpow])
    ml1_lpow_pow  =Multpower_straight(ml1_lpow[0],(l-1))
    ml2_lpow_pow  =Multpower_straight(ml2_lpow[0],(l-1))
    muden_lpow_pow=Multpower_straight(m_den      ,2*(l-1))
    return ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow




def EvalOne(tc_0:Coord,basis:list,pt_list:list,lmd_data:list):
    l=basis[0].order
    [lmd_pow_product,lm_1_lsq,lm_2_lsq]=lmd_data
    assert(type(lmd_pow_product)==dict)
    assert(type(lm_1_lsq)==list)
    assert(type(lm_2_lsq)==list)
    assert(len(basis)  ==3)
    assert(len(pt_list)==3)
    assert(len(lmd_pow_product)==l**2)
    xplincom=XpLinCom(tc_0,basis,pt_list)
    assert(len(xplincom)==l**2)
    [tc_e1,tc_e2,_]=basis
    tc_x=pt_list[0]
    xple1=tc_0.Diff_Add(xplincom[(l-1,0)],tc_e1,xplincom[(l-2,0)])
    xple2=tc_0.Diff_Add(xplincom[(0,l-1)],tc_e2,xplincom[(0,l-2)])
    ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Take_lmd_power(tc_x,xple1,xple2,lm_1_lsq,lm_2_lsq,l)
    xplincom,den=Dict_common_denom(xplincom)
    den_lpow=den**l
    xplincom_lpow={key:Coord([(xplincom[key].numer[i])**l for i in range(0,4)],den_lpow) for key in xplincom.keys()}
    for (k1,k2) in xplincom:
        coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
        xplincom_lpow[(k1,k2)]=xplincom_lpow[(k1,k2)].Mult_frac(coeff)
    xplincom_lpow,den_lpow=Dict_common_denom(xplincom_lpow)
    K=tc_0.field
    pre_tc_fx=[K(0),K(0),K(0),K(0)]
    for k1 in range(0,l):
        for k2 in range(0,l):
            assert(xplincom_lpow[(k1,k2)].denom==den_lpow)
            for i in range(0,4):
                pre_tc_fx[i]+=xplincom_lpow[(k1,k2)].numer[i]
    tc_fx=Coord(pre_tc_fx,den_lpow)
    return tc_fx





def EvalOne_test(tc_0:Coord,basis:list,pt_list:list,lmd_data:list):
    l=basis[0].order
    [lmd_pow_product,lm_1_lsq,lm_2_lsq]=lmd_data
    assert(type(lmd_pow_product)==dict)
    assert(type(lm_1_lsq)==list)
    assert(type(lm_2_lsq)==list)
    assert(len(basis)  ==3)
    assert(len(pt_list)==3)
    assert(len(lmd_pow_product)==l**2)
    xplincom=XpLinCom(tc_0,basis,pt_list)
    assert(len(xplincom)==l**2)
    [tc_e1,tc_e2,_]=basis
    tc_x=pt_list[0]
    xple1=tc_0.Diff_Add(xplincom[(l-1,0)],tc_e1,xplincom[(l-2,0)])
    xple2=tc_0.Diff_Add(xplincom[(0,l-1)],tc_e2,xplincom[(0,l-2)])
    ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Take_lmd_power(tc_x,xple1,xple2,lm_1_lsq,lm_2_lsq,l)
    xplincom,den=Dict_common_denom(xplincom)
    den_lpow=den**l
    xplincom_lpow={key:Coord([(xplincom[key].numer[i])**l for i in range(0,4)],den_lpow) for key in xplincom.keys()}
    for (k1,k2) in xplincom:
        coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
        xplincom_lpow[(k1,k2)]=xplincom_lpow[(k1,k2)].Mult_frac(coeff)
    xplincom_lpow,den_lpow=Dict_common_denom(xplincom_lpow)
    K=tc_0.field
    pre_tc_fx=[K(0),K(0),K(0),K(0)]
    for k1 in range(0,l):
        for k2 in range(0,l):
            assert(xplincom_lpow[(k1,k2)].denom==den_lpow)
            for i in range(0,4):
                pre_tc_fx[i]+=xplincom_lpow[(k1,k2)].numer[i]
    tc_fx=Coord(pre_tc_fx,den_lpow)
    return xplincom_lpow




#2
#avoid to take l-th square roots.
def EvalSq(tc_0:NullCoord,basis:list,pt_list:list,lmd_data:list):
    l=basis[0].order
    [lmd_pow_product,lm_1_lsq,lm_2_lsq]=lmd_data
    assert(type(lmd_pow_product)==dict)
    assert(type(lm_1_lsq)==list)
    assert(type(lm_2_lsq)==list)
    [tc_e1,tc_e2,tc_e12]=basis
    [tc_x,tc_xpe1,tc_xpe2]=pt_list
    a_u=Sum_of_square_mod24(l)
    max_au=max(set(a_u))
    r=len(a_u)    
    auxplincom_list=[] #a_u*x+linear combination for 0<=u<r.
    mult_e1  =tc_0.Mult_Mult(tc_e1  ,max_au)
    mult_e2  =tc_0.Mult_Mult(tc_e2  ,max_au)
    mult_e12 =tc_0.Mult_Mult(tc_e12 ,max_au)
    mult_x   =tc_0.Mult_Mult(tc_x   ,max_au)
    mult_xpe1=tc_0.Mult_Mult(tc_xpe1,max_au)
    mult_xpe2=tc_0.Mult_Mult(tc_xpe2,max_au)
    for u in range(0,r):
        au=a_u[u]
        au_e1  =mult_e1  [au]
        au_e2  =mult_e2  [au]
        au_e12 =mult_e12 [au]
        au_x   =mult_x   [au]
        au_xpe1=mult_xpe1[au]
        au_xpe2=mult_xpe2[au]
        au_e1.order =l
        au_e2.order =l
        au_e12.order=l
        if u!=0 and au==a_u[u-1]:
            auxplincom=auxplincom_list[u-1]
        else:
            auxplincom=XpLinCom(tc_0,[au_e1,au_e2,au_e12],[au_x,au_xpe1,au_xpe2])
        assert(len(auxplincom.keys())==l**2)
        if au==1:
            xplincom=auxplincom
        auxplincom_list.append(auxplincom)
    if all([a_u[u]!=1 for u in range(0,r)]):
        xple1=tc_0.Kxpy_xpy(l,tc_e1,tc_x,tc_xpe1)
        xple2=tc_0.Kxpy_xpy(l,tc_e2,tc_x,tc_xpe2)
    else:
        xple1=tc_0.Diff_Add(xplincom[(l-1,0)],tc_e1,xplincom[(l-2,0)])
        xple2=tc_0.Diff_Add(xplincom[(0,l-1)],tc_e2,xplincom[(0,l-2)])
    ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Take_lmd_power(tc_x,xple1,xple2,lm_1_lsq,lm_2_lsq,l)
    K=tc_x.field
    pre_tc_fx=[[K(0),K(1)],[K(0),K(1)],[K(0),K(1)],[K(0),K(1)]]
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
            for i in range(0,4):
                plus_term_num=coeff[0]*prod([auxplincom_list[u][(k1,k2)].numer[(a_u[u]%2)*i] for u in range(0,r)])
                plus_term_den=coeff[1]*prod([auxplincom_list[u][(k1,k2)].denom               for u in range(0,r)])
                pre_tc_fx[i]=Frac_add(pre_tc_fx[i],[plus_term_num,plus_term_den])
    pre_tc_fx,den=Common_denom_frac(pre_tc_fx)
    tc_fx=Coord([pre_tc_fx[i][0] for i in range(0,4)],den)
    return tc_fx
    



"""
def Evaluation_E3(tc_0:NullCoord,basis:list,pt_list:list,lmd_data:list):
    l=basis[0].order
    [lmd_pow_product,lm_1_lsq,lm_2_lsq]=lmd_data
    assert(type(lmd_pow_product)==dict)
    assert(type(lm_1_lsq)==list)
    assert(type(lm_2_lsq)==list)
    a_u=Sum_of_square(l)
    max_au=max(set(a_u))
    r=len(a_u)
    xplincom=XpLinCom(tc_0,basis,pt_list)
    assert(len(xplincom)==l**2)
    [tc_e1,tc_e2,_]=basis
    tc_x=pt_list[0]    
    xple1=tc_0.Diff_Add(xplincom[(l-1,0)],tc_e1,xplincom[(l-2,0)])
    xple2=tc_0.Diff_Add(xplincom[(0,l-1)],tc_e2,xplincom[(0,l-2)])
    ml1_lpow_pow,ml2_lpow_pow,muden_lpow_pow=Take_lmd_power(tc_x,xple1,xple2,lm_1_lsq,lm_2_lsq,l) 
    K=tc_0.field
    pre_tc_fx=[[K(0),K(1)],[K(0),K(1)],[K(0),K(1)],[K(0),K(1)]]
    for k1 in range(0,l):
        for k2 in range(0,l):
            coeff=[lmd_pow_product[(k1,k2)][0]*ml1_lpow_pow[k1]*ml2_lpow_pow[k2],lmd_pow_product[(k1,k2)][1]*muden_lpow_pow[k1+k2]]
            mult_xpk1e1k2e2=tc_0.Mult_Mult(xplincom[(k1,k2)],max_au)
            au_k1k2=[mult_xpk1e1k2e2[a_u[u]] for u in range(0,r)]
            for i in range(0,4):
                plus_term_num=coeff[0]*prod([au_k1k2[u].numer[(a_u[u]%2)*i] for u in range(0,r)])
                plus_term_den=coeff[1]*prod([au_k1k2[u].denom               for u in range(0,r)])
                pre_tc_fx[i]=Frac_add(pre_tc_fx[i],[plus_term_num,plus_term_den])
    pre_tc_fx,den=Common_denom_frac(pre_tc_fx)
    tc_fx=Coord([pre_tc_fx[i][0] for i in range(0,4)],den)
    return tc_fx
"""


