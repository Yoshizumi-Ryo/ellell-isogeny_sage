
# ellell_isogeny

In this file, you can count operations of calculating $(\ell,\ell)$-isogeny and implement an attack for B-SIDH.

This is writen in [SageMath](https://www.sagemath.org).

## How to use

First, please load `Main.sage`: 

```
load("Main.sage")
```


### 1.Example

When you want to calculate  an $(\ell,\ell)$-isogeny, please look at `test_example.sage`.

In the file, first you implement "setting", then compute  $(\ell,\ell)$-isogeny, for example, 

```
CodSq(tc_0,[tc_e1,tc_e2,tc_e12])
```

### 2.Attack for B-SIDH

When you want to implement attack for B-SIDH, please load `test_attack.sage`.

```
load("test_attack.sage")
```

Then, key exchange between Alice and Bob and attack for it begin. 

If you want to implement other parameters, please rewrite parameter of the above of `test_attack.sage`.


### 3.Count operation to calculate $(\ell,\ell)$-isogeny 

If you want to count the operation of $(\ell,\ell)$-isogeny for 4 ways $\mathtt{CodSq}, \mathtt{CodOne}, \mathtt{EvalSq}, \mathtt{EvalOne}$ for $3\le \ell\le L$, 
please use function `Count_operation_of_isogeny(L)`. You can implement for $L < 200$, for example

```
Count_operation_of_isogeny(30)
```
More concletely, please see `test_count.sage`.



