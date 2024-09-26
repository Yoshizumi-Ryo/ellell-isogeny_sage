
#First, import the following.
import time
import itertools
import sys
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




if __name__ == "__main__":
    args = sys.argv
    assert(len(args) == 3)
    if args[1]=="attack":
        security_bits=int(args[2])
        assert(security_bits==30 or security_bits==128)
        if security_bits==30:
            p=276154505650672190920223
            N_A=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101
            N_B=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79
        if security_bits==128:
            print("Remar that it takes about 11 hours. ")
            p=0x1E409D8D53CF3BEB65B5F41FB53B25EBEAF37761CD8BA996684150A40FFFFFFFF
            N_A=3^(56)* 31*43*59*271*311*353*461*593*607*647*691*743*769*877*1549
            N_B=2^(32)*5^(21)*7*11*163*1181*2389*5233*8353*10139*11939*22003*25391*41843 
        attack_time,attack_result=BSIDH_construct_attack(p,N_A,N_B)
        print("Attack successful:",attack_result)
        print("Attack time (sec):",attack_time)

    if args[1]=="count":
        print("Please wait 20 seconds.")
        max_L=int(args[2])
        assert(max_L<=200)
        p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859
        count_max,N_A,pt_data,K=Count_prepare(p,3,max_L)
        Count_Cod   (max_L,N_A,pt_data,K,"CodSq") 
        print("")
        Count_Cod   (max_L,N_A,pt_data,K,"CodOne") 
        print("")
        Count_Eval(3,max_L,N_A,pt_data,K,"EvalSq") 
        print("")
        Count_Eval(3,max_L,N_A,pt_data,K,"EvalOne")
        print("")


