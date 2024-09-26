
#If you count the operation of isogeny.
#count opertions of isogeny------------------------------------------------------------



p=826791736418446924644415105270960270928927659729776400179861442336062222833458285859


#first.
count_max,N_A,pt_data,K=Count_prepare(p,3,max_L)

#second.
ell_list,csq_list =Count_Cod(max_L,N_A,pt_data,K,"CodSq") 
ell_list,cone_list=Count_Cod(max_L,N_A,pt_data,K,"CodOne") 

ell_list,esq_list =Count_Eval(3,max_L,N_A,pt_data,K,"EvalSq") 
ell_list,eone_list=Count_Eval(3,max_L,N_A,pt_data,K,"EvalOne")




#represent graph of the count.
fig,ax=plt.subplots()
ax.plot(ell_list,csq_list ,color='r')
ax.plot(ell_list,cone_list,color='g')
plt.show()

#represent graph of the count.
fig,ax=plt.subplots()
ax.plot(ell_list,esq_list ,color='r')
ax.plot(ell_list,eone_list,color='b')
plt.show()
#--------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------
#for 211<=ell<=307.

p=321767516318489570622094452273532938297547

count_max,N_A,pt_data,K=Count_prepare(p,211,307)
#second.
ell_list,e1_list=Count_Eval(250,307,N_A,pt_data,K,"EvalSq")

ell_list,e2_list=Count_Eval(250,307,N_A,pt_data,K,"EvalOne")

#--------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------
#for 307<=l<600 s.t. l==1 mod 4. 


p=252*10440816211446376383645271540652922930448181682727450202777-1


count_max,N_A,pt_data,K=Count_prepare_1mod4(p,307,600)


ell_list,e1_list=Count_Eval_1mod4(307,600,N_A,pt_data,K,"EvalSq")

ell_list,e2_list=Count_Eval_1mod4(307,600,N_A,pt_data,K,"EvalOne")

#--------------------------------------------------------------------------------------