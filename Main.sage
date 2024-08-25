


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




