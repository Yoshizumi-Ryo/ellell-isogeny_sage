
# ell_ell_isogeny_sage

In this file, you can count operations of calculating $(\ell,\ell)$-isogeny and implement the SIDH attack on B-SIDH.

This is writen in [SageMath](https://www.sagemath.org).

## 1.Counting the number of arithmetic operations.

You can count the number of arithmetic operations for four alogithms $\mathtt{CodSq}, \mathtt{CodOne}, \mathtt{EvalSq}, \mathtt{EvalOne}$ for prime numbers $3\le \ell\le L$ where 
$3\le L\le 200$. 

Write the following command:

```
sage main.sage "count" {the above L}
```

For example, 

```
$ sage main.sage "count" 20
Please wait 20 seconds.
CodSq ell= 3   1071
CodSq ell= 5   2711
CodSq ell= 7   6740
CodSq ell= 11   14924
CodSq ell= 13   18579
CodSq ell= 17   31829
CodSq ell= 19   44376

CodOne ell= 3   771
CodOne ell= 5   2452

...

EvalOne ell= 17   63336
EvalOne ell= 19   83431
```

## 2.The SIDH Attack on B-SIDH.

You can implement the SIDH attack on B-SIDH with two types of parameter: one with 30 security bits and the other with 128.
Remark that for 128 security bits parameter, it takes about 11 hours.

Write the following command:
```
sage main.sage "attack" {30 or 128}
```

For example, 
```
$ sage main.sage "attack" 30
isogeny chain: [13, 23, 37, 43, 47, 47, 59, 61, 71, 73, 97, 101]
ell= 13
ell= 23

...

ell= 97
ell= 101
Attack successful: True
Attack time (sec): 60.96083307266235
```




