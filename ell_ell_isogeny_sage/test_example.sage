



#setting========================================================

#input paremeter here. -----------------
p =69504748411397252246297776661471 #characteristic.
l=11  #degree.
#----------------------------------------

assert(l//(p+1)==0 or l//(p-1)==0)
K=GF(p**4)
assert(p%4==3)
N_A=l
N_B=2  #order of the point.
assert((p+1)%N_A==0)
assert((p+1)%N_B==0)
assert(gcd(N_A,N_B)==1)
#y^2=x(x-1)(x+1)=x^3-x.
E_0=EllipticCurve(K,[K(-1),K(0)])
P_A,Q_A=E_0.torsion_basis(N_A)
_,_,_,tc_0,tc_e1,tc_e2,tc_e12,tc_x,tc_xpe1,tc_xpe2,tc_y,tc_ype1,tc_ype2=Attack_prepare(E_0,E_0,N_A,N_B,P_A,Q_A,Q_A,P_A,K)
tc_e1.order =l
tc_e2.order =l
tc_e12.order=l
#================================================================



#compute codomain.
tc_f0_c1   =Codomain_C1    (tc_0,[tc_e1,tc_e2,tc_e12])
tc_f0_c1red=Codomain_C1_red(tc_0,[tc_e1,tc_e2,tc_e12])
tc_f0_c2   =Codomain_C2    (tc_0,[tc_e1,tc_e2,tc_e12])
tc_f0=tc_f0_c1


#there are the same theta coordinates.
assert(tc_f0_c1.Is_same_proj(tc_f0_c1red))
assert(tc_f0_c2.Is_same_proj(tc_f0_c1red))
assert(tc_f0_c1.Is_same_proj(tc_f0_c2))


lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])


#compute evaluation.
tc_fx_E1=Evaluation_E1(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
tc_fx_E2=Evaluation_E2(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)


#there should be the same theta coordinates.
assert(tc_fx_E1.Is_same_proj(tc_fx_E2))



#check that their orders are correct.
assert(tc_f0.Is_order(tc_fx_E1,N_B))
assert(tc_f0.Is_order(tc_fx_E2,N_B))




