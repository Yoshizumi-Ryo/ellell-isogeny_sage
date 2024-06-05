
#If you count the operation of isogeny.
#count opertions of isogeny------------------------------------------------------------



p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859


#first.
count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K=Count_operation_of_isogeny_former(p,3,199)

#second.
ell_list,c1_list=Count_later_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1) #(C1)
ell_list,c2_list=Count_later_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2) #(C2)
ell_list,c3_list=Count_later_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3) #(C3)

ell_list,e1_list=Count_later_evaluation(3,39,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1) #(E1)
ell_list,e2_list=Count_later_evaluation(3,199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2) #(E2)
ell_list,e3_list=Count_later_evaluation(3,199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3) #(E3)



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






#--------------------------------------------------------------------------------------
#for 211<=ell<=307.

p=321767516318489570622094452273532938297547

count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K=Count_operation_of_isogeny_former(p,211,307)
#second.
ell_list,e1_list=Count_later_evaluation(250,307,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,e2_list=Count_later_evaluation(250,307,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)

#--------------------------------------------------------------------------------------



p=547*2-1


count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K=Count_operation_of_isogeny_former(p,547,547)




ell_list,e1_list=Count_later_evaluation(503,503,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,e2_list=Count_later_evaluation(503,503,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)



for k in range(2,100):
    if is_prime(547*k-1):
        print(k)
    
