







def Count_operation_of_isogeny(count_max:int):
    assert(count_max<=199)
    p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859
    print("p+1=", factor(p+1))
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(3,count_max+1) if is_prime(l)])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    Fp2=field_K.subfield(2)
    E_0m=EllipticCurve(Fp2,[Fp2(-1),Fp2(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    P_B,Q_B=E_0m.torsion_basis(N_B)
    E_0m=E_0m.base_extend(field_K)
    P_A=E_0m(P_A)
    Q_A=E_0m(Q_A)
    P_B=E_0m(P_B)
    Q_B=E_0m(Q_B)
    _,_,_,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,_,_,_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    assert(tc_0.Is_order(tc_f1  ,N_A))
    assert(tc_0.Is_order(tc_f2  ,N_A))
    assert(tc_0.Is_order(tc_f12 ,N_A))
    #-------------------------------------------------------------------------
    ell_list=[]
    c1_list=[]
    c2_list=[]
    c3_list=[]
    e1_list=[]
    e2_list=[]
    e3_list=[]
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
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            print("r=",len(Sum_of_square(l)))
            #C1----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C1(tc_0,[tc_e1,tc_e2,tc_e12])
            c1_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C2----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C2(tc_0,[tc_e1,tc_e2,tc_e12])
            c2_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C3----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C3(tc_0,[tc_e1,tc_e2,tc_e12])
            c3_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #----------------------------------------------------------------
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #E1----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E1(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e1_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E2----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E2(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e2_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E3----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E3(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e3_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)    
            #----------------------------------------------------------------        
    assert(len(ell_list)==len(c1_list)==len(c2_list)==len(c3_list)==len(e1_list)==len(e2_list)==len(e3_list))
    return ell_list,c1_list,c2_list,c3_list,e1_list,e2_list,e3_list










def Count_operation_of_isogeny_former(count_max:int):
    assert(count_max<=199)
    p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859
    print("p+1=", factor(p+1))
    print("Please wait 2 minutes.")
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(3,count_max+1) if is_prime(l)])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    Fp2=field_K.subfield(2)
    E_0m=EllipticCurve(Fp2,[Fp2(-1),Fp2(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    P_B,Q_B=E_0m.torsion_basis(N_B)
    E_0m=E_0m.base_extend(field_K)
    P_A=E_0m(P_A)
    Q_A=E_0m(Q_A)
    P_B=E_0m(P_B)
    Q_B=E_0m(Q_B)
    _,_,_,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,_,_,_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    assert(tc_0.Is_order(tc_f1  ,N_A))
    assert(tc_0.Is_order(tc_f2  ,N_A))
    assert(tc_0.Is_order(tc_f12 ,N_A))
    return count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K





def Count_operation_of_isogeny_latter(count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K):
    ell_list=[]
    c1_list=[]
    c2_list=[]
    c3_list=[]
    e1_list=[]
    e2_list=[]
    e3_list=[]
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
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            print("r=",len(Sum_of_square(l)))
            #C1----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C1(tc_0,[tc_e1,tc_e2,tc_e12])
            c1_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C2----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C2(tc_0,[tc_e1,tc_e2,tc_e12])
            c2_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C3----------------------------------------------------------------
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_f0=Codomain_C3(tc_0,[tc_e1,tc_e2,tc_e12])
            c3_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #----------------------------------------------------------------
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #E1----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E1(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e1_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E2----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E2(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e2_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E3----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt=(tc_0.field).n_mul
            ct_sqr=(tc_0.field).n_sqr
            tc_fx=Evaluation_E3(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            e3_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)    
            #----------------------------------------------------------------        
    assert(len(ell_list)==len(c1_list)==len(c2_list)==len(c3_list)==len(e1_list)==len(e2_list)==len(e3_list))
    return ell_list,c1_list,c2_list,c3_list,e1_list,e2_list,e3_list








def Count_operation_of_isogeny_latter_codomain(count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,number):
    assert(number==1 or number==2 or number==3)
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
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            #C1----------------------------------------------------------------
            if number==1:
                tc_0.Reset_data()
                tc_e1.Reset_data()
                tc_e2.Reset_data()
                tc_e12.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_f0=Codomain_C1(tc_0,[tc_e1,tc_e2,tc_e12])
                c_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C2----------------------------------------------------------------
            if number==2:
                tc_0.Reset_data()
                tc_e1.Reset_data()
                tc_e2.Reset_data()
                tc_e12.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_f0=Codomain_C2(tc_0,[tc_e1,tc_e2,tc_e12])
                c_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #C3----------------------------------------------------------------
            if number==3:
                tc_0.Reset_data()
                tc_e1.Reset_data()
                tc_e2.Reset_data()
                tc_e12.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_f0=Codomain_C3(tc_0,[tc_e1,tc_e2,tc_e12])
                c_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #----------------------------------------------------------------        
    return ell_list,c_list






def Count_operation_of_isogeny_latter_evaluation(count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,number):
    assert(number==1 or number==2 or number==3)
    ell_list=[]
    e_list=[]
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
            tc_xpe1 =tc_0.Kxpy_xpy(k,tc_f1,tc_x,tc_xpf1)#x+e_1
            tc_xpe2 =tc_0.Kxpy_xpy(k,tc_f2,tc_x,tc_xpf2)#x+e_2
            #--------------------
            K.reset_count()
            print("")
            print("l=",l) 
            print("r=",len(Sum_of_square(l)))
            #C1---------------------------------------------------------------
            tc_f0=Codomain_C1(tc_0,[tc_e1,tc_e2,tc_e12])
            #----------------------------------------------------------------
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #E1----------------------------------------------------------------
            if number==1:
                tc_x.Reset_data()
                tc_xpe1.Reset_data()
                tc_xpe2.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_fx=Evaluation_E1(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
                e_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E2----------------------------------------------------------------
            if number==2:
                tc_x.Reset_data()
                tc_xpe1.Reset_data()
                tc_xpe2.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_fx=Evaluation_E2(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
                e_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)
            #E3----------------------------------------------------------------
            if number==3:
                tc_x.Reset_data()
                tc_xpe1.Reset_data()
                tc_xpe2.Reset_data()  
                ct_mlt=(tc_0.field).n_mul
                ct_sqr=(tc_0.field).n_sqr
                tc_fx=Evaluation_E3(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
                e_list.append((tc_0.field).n_mul+(tc_0.field).n_sqr-ct_mlt-ct_sqr)    
            #----------------------------------------------------------------        
    return ell_list,e_list











