
#If you count the operation of isogeny.
#count opertions of isogeny------------------------------------------------------------

#construct data.
ell_list,c1_list,c2_list,c3_list,e1_list,e2_list,e3_list=Count_operation_of_isogeny(199)


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



#If you want to count one by one.------------------------------------------------
#first.
count_max,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K=Count_operation_of_isogeny_former(199)
#second.
ell_list,c1_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,c2_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)
ell_list,c3_list=Count_operation_of_isogeny_latter_codomain(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3)
ell_list,e1_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,1)
ell_list,e2_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,2)
ell_list,e3_list=Count_operation_of_isogeny_latter_evaluation(199,N_A,tc_0,tc_f1,tc_f2,tc_f12,tc_x,tc_xpf1,tc_xpf2,K,3)

#--------------------------------------------------------------------------------------



