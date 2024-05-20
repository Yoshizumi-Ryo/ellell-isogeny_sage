

import time
import itertools
from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
from sage.schemes.elliptic_curves.ell_curve_isogeny import isogeny_codomain_from_kernel
from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
from sage.functions.log import logb
#============================


load("func_elliptic.py")
load("func_additions.py")
load("func_isogeny.py")
load("func_attack.py")

#---------------------

