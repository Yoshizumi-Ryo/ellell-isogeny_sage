



import time
import itertools
import matplotlib.pyplot as plt
import numpy as np
from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
from sage.functions.log import logb
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic


load("func_fraction.py")
load("class_theta.py")
load("func_elliptic.py")
load("func_isogeny.py")
load("class_count.py")
load("func_for_attack.py")



#attack for B-SIDH.-------------------------------------------
p=276154505650672190920223
N_A=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101
N_B=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79

load("test_attack.py")
#-------------------------------------------------------------



#calculation example-----------------------------------------
p =69504748411397252246297776661471
l=11

load("test_example.py")
#-------------------------------------------------------------




#count the arithemetic operations-----------------------------
p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859
#not prepaed. 

load("test_count.py")
#-------------------------------------------------------------


