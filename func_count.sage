


#p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859



#p=321767516318489570622094452273532938297547



def Count_operation_of_isogeny_former(p,count_from,count_max:int):
    for l in range(count_from,count_max+1):
        if is_prime(l):
            assert((p+1)%l==0)
    #print("p+1=", factor(p+1))
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(count_from,count_max+1) if is_prime(l)])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    E_0m=EllipticCurve(field_K,[field_K(-1),field_K(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    _,_,_,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,_,_,_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    return count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K




def Count_operation_of_isogeny_former_1mod4(p,count_from,count_max:int):
    for l in range(count_from,count_max+1):
        if is_prime(l) and (l%4==1):
            assert((p+1)%l==0)
    #print("p+1=", factor(p+1))
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(count_from,count_max+1) if (is_prime(l)and (l%4==1))])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    E_0m=EllipticCurve(field_K,[field_K(-1),field_K(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    _,_,_,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,_,_,_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    return count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K






def Count_later_codomain(count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,kind):
    ell_list=[]
    c_list=[]
    for l in range(3,count_max+1):
        if is_prime(l):
            l=ZZ(l)
            ell_list.append(l)
            k=N_A//l
            tc_e1 =tc_0.Mult(tc_f1 ,k)
            tc_e2 =tc_0.Mult(tc_f2 ,k)
            tc_e12=tc_0.Mult(tc_f12,k)
            tc_e1.order =l
            tc_e2.order =l
            tc_e12.order=l
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="CodSq":
                tc_f0=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            if kind=="CodOne":
                tc_f0=CodOne(tc_0,[tc_e1,tc_e2,tc_e12])
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            c_list.append(3*ct_mlt_1+2*ct_sqr_1-3*ct_mlt_0-2*ct_sqr_0)    
            #c_list.append(ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2)     
    return ell_list,c_list






def Count_later_evaluation(count_from,count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,kind):
    ell_list=[]
    e_list=[]
    for l in range(count_from,count_max+1):
        if is_prime(l):
            l=ZZ(l)
            ell_list.append(l)
            k=N_A//l
            tc_e1 =tc_0.Mult(tc_f1 ,k)
            tc_e2 =tc_0.Mult(tc_f2 ,k)
            tc_e12=tc_0.Mult(tc_f12,k)
            tc_e1.order =l
            tc_e2.order =l
            tc_e12.order=l
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            print("r=",len(Sum_of_square(l)))
            #---------------------------------------------------------------
            tc_f0=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="EvalSq":
                tc_fx=EvalSq(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            if kind=="EvalOne":
                tc_fx=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            e_list.append(ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2)       
    return ell_list,e_list






def Count_later_evaluation_1mod4(count_from,count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,kind):
    ell_list=[]
    e_list=[]
    for l in range(count_from,count_max+1):
        if is_prime(l) and (l%4==1):
            l=ZZ(l)
            ell_list.append(l)
            k=N_A//l
            tc_e1 =tc_0.Mult(tc_f1 ,k)
            tc_e2 =tc_0.Mult(tc_f2 ,k)
            tc_e12=tc_0.Mult(tc_f12,k)
            tc_e1.order =l
            tc_e2.order =l
            tc_e12.order=l
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            print("r=",len(Sum_of_square(l)))
            #---------------------------------------------------------------
            tc_f0=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="EvalSq":
                tc_fx=EvalSq(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            if kind=="EvalOne":
                tc_fx=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            e_list.append(ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2)       
    return ell_list,e_list










