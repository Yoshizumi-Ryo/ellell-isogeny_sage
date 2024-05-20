
"""
from sage.all import (
    GF,
    PolynomialRing,
    ZZ
)
"""



#construct big finite field-------


#return finite field F_{p^4}.
def GFp4pow(p):
    Fp2=GF(p**2)
    z=Fp2(0)
    while z.is_square():
        z=Fp2.random_element()
    t=ZZ(z + z**p)
    n=ZZ(z**(p+1))
    R = PolynomialRing(GF(p),'x')
    x=R.gen()
    f=x**4-t*x**2+n
    Fp4=GF(p**4,modulus=f,name='a')
    return Fp4








# finite field with operations count
def finite_field_with_count(p):
    assert(is_prime(p))
    class Field():
        def __init__(self, x):
            if type(x) == Field:
                self.v = x.v
            else:
                self.v = Field.F(x)

        def __add__(self, x):
            Field.n_add += 1
            if not type(x) == Field:
                return Field(self.v + x)
            else:
                return Field(self.v + x.v)
        
        def __iadd__(self, x):
            return self + x
        
        def __neg__(self):
            return Field(-self.v)

        def __sub__(self, x):
            return self + (-x)

        def __isub__(self, x):
            return self - x

        def __mul__(self, x):
            if not type(x) == Field:
                return x * self
            elif x==1:
                return self
            else:
                Field.n_mul += 1
                return Field(self.v * x.v)

        # multiplication from right, i.e., x*self
        def __rmul__(self, x):
            if not type(x) == Field:
                if x == 0:
                    return Field(0)
                elif x < 1:
                    return (-x)*(-self)
                else:
                    ret = Field(self.v)
                    for b in bin(x)[3:]:
                        ret = ret + ret
                        if b == "1":
                            ret = ret + self
                return ret
            else:
                return self*x

        def __imul__(self, x):
            return self * x

        def inv(self):
            Field.n_inv += 1
            return Field(self.v ** -1)
        
        def __truediv__(self, x):
            return self * x.inv()

        def __rtruediv__(self, x):
            return x * self.inv()

        def __pow__(self, e):
            if e == 0:
                return Field(1)
            elif e < 0:
                return self.inv() ** (-e)
            elif e == 2:
                Field.n_sqr += 1
                return Field(self.v**2)
            else:
                ret = Field(self.v)
                for b in bin(e)[3:]:
                    ret = ret**2
                    if b == "1":
                        ret = ret * self
            return ret

        def __eq__(self, x) -> bool:
            if type(x) == Field:
                return self.v == x.v
            else:
                return self.v == x

        def __str__(self):
            return str(self.v)

        def random_element():
            return Field(Field.F.random_element())

        def reset_count():
            Field.n_add = 0
            Field.n_mul = 0
            Field.n_sqr = 0
            Field.n_inv = 0
        def count():
            return ("M",Field.n_mul),("S",Field.n_sqr),("I",Field.n_inv),("M+S",Field.n_mul+Field.n_sqr)
        
    K=GFp4pow(p)
    Field.F=K
    Field.n_add = 0
    Field.n_mul = 0
    Field.n_sqr = 0
    Field.n_inv = 0
    return Field


    
    






