

#setting==============

assert(l//(p+1)==0 or l//(p-1)==0)

K=GF(p**4)
assert(p%4==3)

N_A=l
N_B=2  #order of the point.
assert((p+1)%N_A==0)
assert((p+1)%N_B==0)
assert(gcd(N_A,N_B)==1)

#y^2=x(x-1)(x+1)=x^3-x.
#field_K=K.F
#Fp2=K.subfield(2)

E_0=EllipticCurve(K,[K(-1),K(0)])
P_A,Q_A=E_0.torsion_basis(N_A)
#P_B,Q_B=E_0.torsion_basis(N_B)
#E_0m=E_0m.base_extend(K)
#P_A=E_0m(P_A)
#Q_A=E_0m(Q_A)
#P_B=E_0m(P_B)
#Q_B=E_0m(Q_B)



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


lmd_pow_product=Product_power_lambda(tc_0,[tc_f1,tc_f2,tc_f12])

#compute evaluation.
tc_fx_E1=Evaluation_E1(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_pow_product)
tc_fx_E2=Evaluation_E2(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_pow_product)
tc_fx_E3=Evaluation_E3(tc_0,[tc_f1,tc_f2,tc_f12],[tc_x,tc_xpf1,tc_xpf2],lmd_pow_product)


#there should be the same theta coordinates, but...
print(tc_fx_E1.Is_same_proj(tc_fx_E2))
print(tc_fx_E2.Is_same_proj(tc_fx_E3))
print(tc_fx_E3.Is_same_proj(tc_fx_E1))



#check that their orders are correct.
print(tc_f0.Is_order(tc_fx_E1,N_B))
print(tc_f0.Is_order(tc_fx_E2,N_B))
print(tc_f0.Is_order(tc_fx_E3,N_B))



