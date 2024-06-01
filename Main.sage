


#First, import the following.
import time
import itertools
import matplotlib.pyplot as plt
import numpy as np
from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
from sage.functions.log import logb
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic


#Second, load the following.
load("func_fraction.sage")
load("class_theta.sage")
load("func_elliptic.sage")
load("func_isogeny.sage")
load("class_count.sage")
load("func_for_attack.sage")
load("func_count.sage")




#If you attack for B-SIDH, load the following.
#attack for B-SIDH.-------------------------------------------
p=276154505650672190920223
N_A=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101
N_B=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79

load("test_attack.sage")
#-------------------------------------------------------------





#If you count the operation of isogeny.
#count opertions of isogeny------------------------------------------------------------

#construct data.
ell_list,c1_list,c2_list,c3_list,e1_list,e2_list,e3_list=Count_operation_of_isogeny(199)


#represent graph of the count.
fig,ax=plt.subplots()
ax.plot(ell_list,c1_list,color='r')
ax.plot(ell_list,c2_list,color='b')
ax.plot(ell_list,c3_list,color='g')
plt.show()

#represent graph of the count.
fig,ax=plt.subplots()
ax.plot(ell_list,e1_list,color='r')
ax.plot(ell_list,e2_list,color='b')
ax.plot(ell_list,e3_list,color='g')
plt.show()
#--------------------------------------------------------------------------------------

#If you want to count one by one.------------------------------------------------
#first.
count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K=Count_operation_of_isogeny_former(199)
#second.
ell_list,c1_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,c2_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)
ell_list,c3_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3)
ell_list,e1_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,e2_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)
ell_list,e3_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3)

#--------------------------------------------------------------------------------------



