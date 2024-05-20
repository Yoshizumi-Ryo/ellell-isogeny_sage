



#setting.====================
N_A=3*l
N_B=5
#y^2=x(x-1)(x+1)=x^3-x.
K=GFp4pow(p)
E_0m=EllipticCurve(K,[K(-1),K(0)])
P_A,Q_A=E_0m.torsion_basis(N_A)
P_B,Q_B=E_0m.torsion_basis(N_B)
E_0m,PA_EB,QA_EB,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2,lv2_ZS2,lv2_ZS2pf1,lv2_ZS2pf2=attack_prepare(E_0m,E_0m,N_A,N_B,Q_A,P_A,P_A,Q_A)
lv2tnp,lv2_e1,lv2_e2,lv2_x=attack_isogeny_for_time(N_A,lv2tnp,lv2_f1,lv2_f2,lv2_f12,lv2_ZS1,lv2_ZS1pf1,lv2_ZS1pf2)
assert(not Is_Ellipitc_product(lv2tnp)[0])
#=================================


print("al_sumsq"),
isogeny_time(lv2tnp,lv2_e1,lv2_e2,lv2_x,l,"al_sumsq") #(i)

print("al_lpow"),
isogeny_time(lv2tnp,lv2_e1,lv2_e2,lv2_x,l,"al_lpow")  #(ii)

#print("tl_sumsq"),
#isogeny_time(lv2tnp,lv2_e1,lv2_e2,lv2_x,l,"tl_sumsq") #(iii)

#print("tl_lpow"),
#isogeny_time(lv2tnp,lv2_e1,lv2_e2,lv2_x,l,"tl_lpow")  #(iv)





#---------------



for l in range(3,300):
    if is_prime(l):
        l=ZZ(l)
        a=sum_of_square(l)
        r=len(a)
        set1=Set([])
        for i in range(0,r):
            set1={(((a[i]*k_1)%l),((a[i]*k_2)%l)) for k_1 in range(0,l) for k_2 in range(0,l)}.union(set1)
        len(set1)==l**2




for l in range(5,100):
    type(l)


a=sum_of_square(l)
r=len(a)
set1=Set([])
for i in range(0,r):
    set1={((a[i]*k_1)%l,(a[i]*k_2)%l) for k_1 in range(0,l) for k_2 in range(0,l)}.union(set1)
len(set1)==l**2









