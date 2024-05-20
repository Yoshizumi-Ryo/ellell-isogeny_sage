


p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859


#field with counter.
K=finite_field_with_count(p)


#setting.====================
N_A=prod([l for l in range(3,55) if is_prime(l)])
N_B=101  #order of the point.
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
tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpe1,tc_xpe2,tc_y,tc_ype1,tc_ype2=Sample_construct(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A,K)
#-----------------------


for l in range(3,55):
    if is_prime(l):
        l=ZZ(l)
        r_1=N_A//l
        tc_e1 =tc_0.Mult(tc_f1 ,r_1)
        tc_e2 =tc_0.Mult(tc_f2 ,r_1)
        tc_e12=tc_0.Mult(tc_f12,r_1)
        tc_e1.order =l
        tc_e2.order =l
        tc_e12.order=l
        #--------------------
        #reset.
        K.reset_count()
        print("l=",l)
        print("r=",l%4)
        tc_e1.Reset_data()
        tc_e2.Reset_data()
        tc_e12.Reset_data()  
        print("a")
        tc_f0_a=Codomain_al_sumsq_a(tc_0,tc_e1,tc_e2,tc_e12)
        tc_e1.Reset_data()
        tc_e2.Reset_data()
        tc_e12.Reset_data() 
        print("b") 
        tc_f0_b=Codomain_al_sumsq_b(tc_0,tc_e1,tc_e2,tc_e12)
        tc_e1.Reset_data()
        tc_e2.Reset_data()
        tc_e12.Reset_data()  
        print("c")    
        tc_f0_c=Codomain_al_lpow_c (tc_0,tc_e1,tc_e2,tc_e12)
        print("")






l_list=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]

A=[305,1442,2980,11193,17008,23684,49583,79919,113234,134008,161644,274652,404621,374963,376976]

B=[376,1186,2628,5727,6963,11567,15949,25055,32141,44520,51963,63407,77275,101051,105035]

C=[275,1149,2266,5860,7959,13376,17627,26731,42106,50320,69425,83753,96698,120784,146433]


fig,ax=plt.subplots()


ax.plot(l_list,A,color='r')
ax.plot(l_list,B,color='b')
ax.plot(l_list,C,color='g')
plt.show()



#-----------------------
#codomain.

tc_f0_a=Codomain_al_sumsq_a(tc_0,tc_e1,tc_e2,tc_e12)

tc_f0_b=Codomain_al_sumsq_b(tc_0,tc_e1,tc_e2,tc_e12)

tc_f0_c=Codomain_al_lpow_c (tc_0,tc_e1,tc_e2,tc_e12)


#check.

tc_f0=tc_f0_a
tc_f0_a.Is_same_proj(tc_f0)
tc_f0_b.Is_same_proj(tc_f0)
tc_f0_c.Is_same_proj(tc_f0)


#--------------
#eveluation.

tc_f1pe1=tc_0.Mult(tc_f1,r_1+1)
tc_f1pe2=tc_0.Kxpy_xpy(r_1,tc_f2,tc_f1,tc_f12)




#================================================================
lmd1_lpow =tc_0.Lmd_lpow(tc_e1)
lmd2_lpow =tc_0.Lmd_lpow(tc_e2)
lmd12_lpow=tc_0.Lmd_lpow(tc_e12)
lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
[lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
den_pow         =Power_simultaneously(den            ,3*(l-1)**2)
lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
lmd2_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],(l-1)**2)
lmd_pow_product=dict()
for k1 in range(0,l):
    for k2 in range(0,l):
        lmd_pow_product[(k1,k2)]=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow [k1**2+k2**2+k1*k2]]
#================================================================










tc_fx_a=Evaluation_al_sumsq_a(tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2,lmd_pow_product)

tc_fx_b=Evaluation_al_sumsq_b(tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2)

tc_fx_c=Evaluation_al_lpow_c (tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2)

#check.
tc_f0.Is_order(tc_fx_a,r_1)
tc_f0.Is_order(tc_fx_b,r_1)
tc_f0.Is_order(tc_fx_c,r_1)
#----------------------------------------------------------------



for l in range(120,140):
    if is_prime(l):
        if (l%4==1):
            print(l)
            l=ZZ(l)
            r_1=N_A//l
            tc_e1 =tc_0.Mult(tc_f1 ,r_1)
            tc_e2 =tc_0.Mult(tc_f2 ,r_1)
            tc_e12=tc_0.Mult(tc_f12,r_1)
            tc_e1.order =l
            tc_e2.order =l
            tc_e12.order=l                
            tc_f1pe1=tc_0.Mult(tc_f1,r_1+1)
            tc_f1pe2=tc_0.Kxpy_xpy(r_1,tc_f2,tc_f1,tc_f12)
            #--------------------
            #reset.
            K.reset_count()
            print("l=",l)
            print("r=",l%4)  
            #================================================================
            lmd1_lpow =tc_0.Lmd_lpow(tc_e1)
            lmd2_lpow =tc_0.Lmd_lpow(tc_e2)
            lmd12_lpow=tc_0.Lmd_lpow(tc_e12)
            lmd_div_lpow=[lmd12_lpow[0]*lmd1_lpow[1]*lmd2_lpow[1],lmd12_lpow[1]*lmd1_lpow[0]*lmd2_lpow[0]]
            [lmd1_lpow,lmd2_lpow,lmd_div_lpow],den=Common_denom_frac([lmd1_lpow,lmd2_lpow,lmd_div_lpow])
            den_pow         =Power_simultaneously(den            ,3*(l-1)**2)
            lmd1_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
            lmd2_lpow_pow   =Power_simultaneously(lmd1_lpow[0]   ,(l-1)**2)
            lmd_div_lpow_pow=Power_simultaneously(lmd_div_lpow[0],(l-1)**2)
            lmd_pow_product=dict()
            for k1 in range(0,l):
                for k2 in range(0,l):
                    lmd_pow_product[(k1,k2)]=[lmd1_lpow_pow[k1**2]*lmd2_lpow_pow[k2**2]*lmd_div_lpow_pow[k1*k2],den_pow [k1**2+k2**2+k1*k2]]
            #===============================================================
            #tc_fx_a=Evaluation_al_sumsq_a(tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2)
            tc_fx_b=Evaluation_al_sumsq_b(tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2,lmd_pow_product)
            tc_fx_c=Evaluation_al_lpow_c (tc_0,tc_e1,tc_e2,tc_e12,tc_f1,tc_f1pe1,tc_f1pe2,lmd_pow_product)
            print("")




l_list=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]



A=[1451,3861,7422,24831,36879,52331,106536,172869,243204,293460,352480,586647,856603, 814525,821610]
B=[1853,3864,12004,24994,27605,48523,78728,141910,152538,264197,256886,320296,440168,632064,554419]
C=[929,2484,4857,12564,17531,30375,39933,61766,99400,118888,164249,203038,233269,291191,358809]



fig,ax=plt.subplots()


ax.plot(l_list,A,color='r')
ax.plot(l_list,B,color='b')
ax.plot(l_list,C,color='g')
plt.show()

