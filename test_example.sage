
#examples of isogeny calculations


#setting========================================================

#input parameter here. -----------------
p =69504748411397252246297776661471 #characteristic.

l=11  #degree.
#----------------------------------------

assert((p+1)%l==0)
K=GF(p**2)
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
tc_f0_sq =CodSq (tc_0,[tc_e1,tc_e2,tc_e12])
tc_f0_one=CodOne(tc_0,[tc_e1,tc_e2,tc_e12])
tc_f0=tc_f0_sq

#there are the same theta coordinates.
assert(tc_f0_sq.Is_same_proj(tc_f0_one))


lmd_data=Product_power_lambda([tc_e1,tc_e2,tc_e12])


#compute evaluation.
tc_fx_sq =EvalSq (tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)
tc_fx_one=EvalOne(tc_0,[tc_e1,tc_e2,tc_e12],[tc_x,tc_xpe1,tc_xpe2],lmd_data)


#there should be the same theta coordinates.
assert(tc_fx_sq.Is_same_proj(tc_fx_one))



#check that their orders are correct.
assert(tc_f0.Is_order(tc_fx_sq ,N_B))
assert(tc_f0.Is_order(tc_fx_one,N_B))








