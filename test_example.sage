



#setting========================================================
p =69504748411397252246297776661471
l=11
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
_,_,_,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,tc_y,tc_ypf1,tc_ypf2=Attack_prepare(E_0,E_0,N_A,N_B,P_A,Q_A,Q_A,P_A,K)
tc_f1.order =l
tc_f2.order =l
tc_f12.order=l
#================================================================



#compute codomain.
tc_f0_c1=Codomain_C1(tc_0,[tc_f1,tc_f2,tc_f12])
tc_f0_c2=Codomain_C2(tc_0,[tc_f1,tc_f2,tc_f12])
tc_f0_c3=Codomain_C3(tc_0,[tc_f1,tc_f2,tc_f12])
tc_f0=tc_f0_c1


#there are the same theta coordinates.
assert(tc_f0_c1.Is_same_proj(tc_f0_c2))
assert(tc_f0_c2.Is_same_proj(tc_f0_c3))
assert(tc_f0_c3.Is_same_proj(tc_f0_c1))


lmd_data=Product_power_lambda([tc_f1,tc_f2,tc_f12])

tc_0.Reset_data()
tc_f1.Reset_data()
tc_f2.Reset_data()
tc_f12.Reset_data()  

#compute evaluation.
tc_fx_E1=Evaluation_E1(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_data)
tc_fx_E2=Evaluation_E2(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_data)
tc_fx_E3=Evaluation_E3(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_data)


#there should be the same theta coordinates.
assert(tc_fx_E1.Is_same_proj(tc_fx_E2))
assert(tc_fx_E2.Is_same_proj(tc_fx_E3))
assert(tc_fx_E3.Is_same_proj(tc_fx_E1))



#check that their orders are correct.
assert(tc_f0.Is_order(tc_fx_E1,N_B))
assert(tc_f0.Is_order(tc_fx_E2,N_B))
assert(tc_f0.Is_order(tc_fx_E3,N_B))



