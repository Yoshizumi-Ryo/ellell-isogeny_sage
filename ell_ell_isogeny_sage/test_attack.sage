  

#input parameter here. ---------------------------------------------

p=276154505650672190920223
N_A=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101
N_B=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79


'''
#propsed parameter.
p=0x1E409D8D53CF3BEB65B5F41FB53B25EBEAF37761CD8BA996684150A40FFFFFFFF
N_A=3^(56)* 31*43*59*271*311*353*461*593*607*647*691*743*769*877*1549
N_B=2^(32)*5^(21)*7*11*163*1181*2389*5233*8353*10139*11939*22003*25391*41843
'''
#--------------------------------------------------------------------


print("p=",p)
print("bit length of p=",len(ZZ(p).digits(2)))
print("N_A=",N_A)
print("N_B=",N_B)

assert((p-1)%N_A==0)
assert((p+1)%N_B==0)



#the 4th primitive root of 1=========
K,_,Fp2=GFp4pow(p)
X=gen(K['X'])
f=X**2+1
zeta_4=f.roots()[0][0]
assert(zeta_4**2==-1)
#===================================


#public setting.====================
print("public setting")
#y^2=x(x-1)(x+1)=x^3-x.
E_0m=EllipticCurve(K,[K(-1),K(0)])
P_A,Q_A=E_0m.torsion_basis(N_A)
P_B,Q_B=E_0m.torsion_basis(N_B)
#==================================



#Bob calculates secretly.===========
print("Bob calculates secretly.")
coff_B=randint(0,(N_B-1))
R_B=P_B+coff_B*Q_B
assert(order(R_B)==N_B)
E_B,PA_EB,QA_EB=Elliptic_Cyclic(E_0m,R_B,P_A,Q_A)
assert(order(PA_EB)==N_A)
assert(order(QA_EB)==N_A)
#===================================

time_start=time.time()

#Construction auxiliary path.=========
print("Construction auxiliary path.")
E_0p=EllipticCurve(K,[K(1),K(0)])
iso_mtop=E_0m.isomorphisms(E_0p)[1]
E_pr,alpha_PA,alpha_QA=Auxiliary_path(E_0p,N_A,N_B,iso_mtop(P_A),iso_mtop(Q_A))
#=====================================


#Main attack.=========================
print("Main attack.")
Attacker_ker=Attack_total(E_0m,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_PA,alpha_QA,zeta_4,K)
#=====================================


time_end=time.time()
time_diff=time_end-time_start


#Check if this attacke succeed.========
print("This attack is a success.")
assert(order(Attacker_ker)==N_B)
print(Attacker_ker.weil_pairing(R_B,N_B)==1)
#=====================================



print("Time for this attack",time_diff,"seconds.")



