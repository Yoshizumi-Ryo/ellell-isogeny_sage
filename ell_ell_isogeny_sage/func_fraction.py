
#give bianry representation of N with i-digits.
def Binary_Exp(N:int,i:int):
    N_list=ZZ(N).digits(2)
    N_tup=tuple(N_list)
    assert(len(N_tup)<=i)
    for m in range(0,i-len(N_tup)):
        N_tup=N_tup+tuple([0])
    return N_tup
        
    
def B_val(t:tuple):
    return sum([t[j]*(2**j) for j in range(0,len(t))])
    



def Common_Denom(denoms:list):
    assert(all(denoms[i]!=0 for i in range(0,len(denoms))))
    N=len(denoms)    
    d_list=ZZ(N-1).digits(2)
    n=len(d_list)
    a={}
    b={}
    a={Binary_Exp(m,n):denoms[m] for m in range(0,N)}
    for i in range(1,n+1):
        d_i_list=tuple([d_list[j] for j in range(i,n)])
        assert(len(d_i_list)==n-i)
        d_i_val=B_val(d_i_list)
        #(d_i_list==Binary_Exp(d_i_val,n-i))
        for e_val in range(0,d_i_val+1):
            e_list=Binary_Exp(e_val,n-i)
            if e_list==d_i_list and d_list[i-1]==0:
                a[e_list]=a[tuple([0])+e_list]
            else:
                a[e_list]=a[tuple([0])+e_list]*a[tuple([1])+e_list]
    total_product=a[()]
    b[tuple([0])]=a[tuple([1])]
    b[tuple([1])]=a[tuple([0])]
    for i in reversed(range(0,n-1)):
        d_i_list=tuple([d_list[j] for j in range(i,n)])
        d_i_val=B_val(d_i_list)
        for e_val in range(0,d_i_val+1):
            e_list=Binary_Exp(e_val,n-i)
            if e_list==d_i_list and e_list[0]==0 and d_i_list[0]==0:
                r_e_list=tuple([e_list[i] for i in range(1,len(e_list))])
                b[e_list]=b[r_e_list]
            else:
                r_e_list=tuple([e_list[i] for i in range(1,len(e_list))])
                c=(e_list[0]+1)%2
                c_e_list=tuple([c])+r_e_list
                b[e_list]=b[r_e_list]*a[c_e_list]
    b_list=[b[Binary_Exp(i,n)] for i in range(0,N)]
    assert(total_product!=0)
    return b_list,total_product




def Common_denom_frac(fracs:list):
    assert(all(len(fracs[i])==2 for i in range(0,len(fracs))))
    assert(all(fracs[i][1]!=0   for i in range(0,len(fracs))))
    denoms=[fracs[i][1] for i in range(0,len(fracs))]
    numes,den=Common_Denom(denoms)
    new_fracs=[[(numes[i]*fracs[i][0]),den] for i in range(0,len(fracs))]
    assert(den!=0)
    return new_fracs,den




def Dict_common_denom(lincom:dict):
    list_lincom=sorted(lincom.items())
    k=len(list_lincom)
    list_val=[(list_lincom[i][1]).denom for i in range(0,k)]
    list_num,den=Common_Denom(list_val)
    dict_num={list_lincom[i][0]:list_num[i] for i in range(0,k)}
    for k in lincom.keys():
        assert(type(lincom[k])==Coord)
        lincom[k].numer=[dict_num[k]*lincom[k].numer[j] for j in range(0,4)]
        lincom[k].denom=den
    assert(den!=0)
    return lincom,den




def Frac_sum(frac1:list,frac2:list):
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    return [frac1[0]*frac2[1]+frac1[1]*frac2[0],frac1[1]*frac2[1]]




def Frac_minus(frac1:list,frac2:list):
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    if frac1[0]==0:
        return [-frac2[0],frac2[1]]
    elif frac2[0]==0:
        return frac1
    elif frac1[1]==frac2[1]:
        return [frac1[0]-frac2[0],frac1[1]]
    else:
        return [frac1[0]*frac2[1]-frac1[1]*frac2[0],frac1[1]*frac2[1]]



def Frac_mult(frac1:list,frac2:list):
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    K=parent(frac1[0])
    if frac1[0]==0:
        return [K(0),K(1)]
    elif frac2[0]==0:
        return [K(0),K(1)]
    elif frac1[0]==frac2[0] and frac1[1]==frac2[1]:
        return [frac1[0]**2,frac1[1]**2]
    elif frac1[0]==frac2[0]:
        return [frac1[0]**2,frac1[1]*frac2[1]]
    elif frac1[1]==frac2[1]:
        return [frac1[0]*frac2[0],frac1[1]**2]
    else:
        return [frac1[0]*frac2[0],frac1[1]*frac2[1]]



def Frac_div(frac1:list,frac2:list):
    assert(len(frac1)==2)
    assert(len(frac2)==2)
    assert(frac1[1]!=0)
    assert(frac2[1]!=0)
    return [frac1[0]*frac2[1],frac1[1]*frac2[0]]




def Projective_Theta(frac_coord:list):
    assert(len(frac_coord)==4)
    [[n1,d1],[n2,d2],[n3,d3],[n4,d4]]=frac_coord
    d12=d1*d2
    d34=d3*d4
    return [n1*d2*d34,d1*n2*d34,d12*n3*d4,d12*d3*n4]




def Power_simultaneously(value,max_exp:int):
    pow_list=[parent(value)(1),value]
    value_pow=value
    for exp in range(2,max_exp+1):
        value_pow*=value
        pow_list.append(value_pow)   
    assert(len(pow_list)==max_exp+1) 
    return pow_list



