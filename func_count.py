

##################################################################################
# Functions for counting the number of arithmetic operations of (ell,ell)-isogeny.
##################################################################################

from sage.all import *
from class_count import Finite_field_with_count
from func_for_attack import Attack_prepare
from func_isogeny import CodOne,CodSq,Product_power_lambda,EvalSq,EvalOne
from class_theta import NullCoord,Coord



def Count_prepare(p:int,count_from:int,count_max:int):
    for l in range(count_from,count_max+1):
        if is_prime(l):
            assert((p+1)%l==0)
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(count_from,count_max+1) if is_prime(l)])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    E_0m=EllipticCurve(field_K,[field_K(-1),field_K(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    _,_,_,tc_0,(tc_f1,tc_f2,tc_f12),(tc_x,tc_xpf1,tc_xpf2),_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    pt_data=[tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2]
    return count_max,N_A,pt_data,K




def Count_prepare_1mod4(p:int,count_from:int,count_max:int):
    for l in range(count_from,count_max+1):
        if is_prime(l) and (l%4==1):
            assert((p+1)%l==0)
    #field with counter.
    K=Finite_field_with_count(p)
    #setting.----------------------------------------------------------------
    N_A=prod([l for l in range(count_from,count_max+1) if (is_prime(l)and (l%4==1))])
    N_B=2  #order of the point.
    #y^2=x(x-1)(x+1)=x^3-x.
    field_K=K.F
    E_0m=EllipticCurve(field_K,[field_K(-1),field_K(0)])
    P_A,Q_A=E_0m.torsion_basis(N_A)
    _,_,_,tc_0,(tc_f1,tc_f2,tc_f12),(tc_x,tc_xpf1,tc_xpf2),_=Attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
    return count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K






def Count_Cod(count_max:int,N_A:int,pt_data:list,K,kind:str):
    """ 
    count the number of arithmetic operations to compute codomain of (ell,ell)-isogeny.
    """
    [tc_0,tc_f1,tc_f2,tc_f12,_,_,_]=pt_data
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
            tc_0.Reset_data()
            tc_e1.Reset_data()
            tc_e2.Reset_data()
            tc_e12.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="CodSq":
                _=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            if kind=="CodOne":
                _=CodOne(tc_0,[tc_e1,tc_e2,tc_e12])
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            result=3*ct_mlt_1+2*ct_sqr_1-3*ct_mlt_0-2*ct_sqr_0
            #result_2=ct_mlt_1+ct_sqr_1-ct_mlt_0-ct_sqr_0
            c_list.append(result) 
            print(kind,"ell=",l," ",result)
            #print(kind,"ell=",l," ",ct_sqr_1-ct_sqr_0,ct_mlt_1-ct_mlt_0) 
            #c_list.append(ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2)     
    return ell_list,c_list






def Count_Eval(count_from:int,count_max:int,N_A:int,pt_data:list,K,kind:str):
    """ 
    count the number of arithmetic operations to compute evaluation of (ell,ell)-isogeny.
    """
    [tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2]=pt_data
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
            #print("")
            #print("l=",l) 
            #print("r=",len(Sum_of_square(l)))
            #---------------------------------------------------------------
            _=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="EvalSq":
                _=EvalSq(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            if kind=="EvalOne":
                _=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            result=ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2
            print(kind,"ell=",l," ",result)  
            #print(kind,"ell=",l," ",ct_sqr_1-ct_sqr_0,ct_mlt_1-ct_mlt_0)  
            e_list.append(result)       
    return ell_list,e_list






def Count_Eval_1mod4(count_from:int,count_max:int,N_A:int,tc_0:NullCoord,tc_f1:Coord,tc_f2:Coord,tc_f12:Coord,tc_x:Coord,tc_xpf1:Coord,tc_xpf2:Coord,K,kind:str):
    """ 
    count the number of arithmetic operations to compute evaluation of (ell,ell)-isogeny.
    For the case of l=1 (mod 4).
    """
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
            _=CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
            lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])
            #----------------------------------------------------------------
            tc_x.Reset_data()
            tc_xpe1.Reset_data()
            tc_xpe2.Reset_data()  
            ct_mlt_0=(tc_0.field).n_mul
            ct_sqr_0=(tc_0.field).n_sqr
            if kind=="EvalSq":
                _=EvalSq(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            if kind=="EvalOne":
                _=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
            ct_mlt_1=(tc_0.field).n_mul
            ct_sqr_1=(tc_0.field).n_sqr
            e_list.append(ct_mlt_1*3+ct_sqr_1*2-ct_mlt_0*3-ct_sqr_0*2)       
    return ell_list,e_list

