

#0,1,2,3=[0,0],[0,1],[1,0],[1,1]
class Coord:
    def __init__(self,numer:list,denom):
        assert(len(numer)==4)
        assert(denom!=0)
        self.numer=numer
        self.order=0
        self.field=parent(numer[0])
        if type(numer[0])!=type(denom):
            denom=self.field(denom)
        self.denom=denom
        assert(type(numer[0])==type(denom))

    
    def Hadamard(self):
        tc=self.numer
        h_num=[tc[0]+tc[1]+tc[2]+tc[3],
               tc[0]-tc[1]+tc[2]-tc[3],
               tc[0]+tc[1]-tc[2]-tc[3],
               tc[0]-tc[1]-tc[2]+tc[3]]
        return Coord(h_num,self.denom)
    
    def Square(self):
        try:
            return self.square
        except AttributeError:
            self.square=Coord([self.numer[i]**2 for i in range(0,4)],self.denom**2)
            return self.square

    #Hadamard of sauare.
    def Tc_Hsq(self):
        try:
            return self.hsq
        except AttributeError:
            self.hsq=(self.Square()).Hadamard()
            return self.hsq
    

    def Zichi(self):
        try:
            return self.zichi
        except AttributeError:
            prod_00_01=self.numer[0]*self.numer[1]
            prod_00_10=self.numer[0]*self.numer[2]
            prod_00_11=self.numer[0]*self.numer[3]
            prod_01_10=self.numer[1]*self.numer[2]
            prod_01_11=self.numer[1]*self.numer[3]
            prod_10_11=self.numer[2]*self.numer[3]
            self.zichi={}
            self.zichi["_10^00"]=2*prod_00_10+2*prod_01_11
            self.zichi["_10^01"]=2*prod_00_10-2*prod_01_11
            self.zichi["_01^00"]=2*prod_00_01+2*prod_10_11
            self.zichi["_01^10"]=2*prod_00_01-2*prod_10_11
            self.zichi["_11^00"]=2*prod_00_11+2*prod_01_10
            self.zichi["_11^11"]=2*prod_00_11-2*prod_01_10
            return self.zichi
   

    def Mult_const(self,lmd):
        return Coord([self.numer[i]*lmd for i in range(0,4)],self.denom)

    def Mult_frac(self,frac:list):
        return Coord([self.numer[i]*frac[0] for i in range(0,4)],self.denom*frac[1])
 
    def Is_same_affine(self,tc_x):
        return all(self.numer[i]*tc_x.denom==tc_x.numer[i]*self.denom for i in range(0,4))
    
    def Is_same_proj(self,tc_x):
        assert(tc_x.numer[0]!=0)
        assert(self.numer[0]!=0)
        assert(type(self.numer[0])==type(tc_x.numer[0]))
        assert(type(self.numer[1])==type(tc_x.numer[1]))
        assert(type(self.numer[2])==type(tc_x.numer[2]))
        assert(type(self.numer[3])==type(tc_x.numer[3]))
        assert(parent(self.numer[0])==parent(tc_x.numer[0]))
        assert(parent(self.numer[1])==parent(tc_x.numer[1]))
        assert(parent(self.numer[2])==parent(tc_x.numer[2]))
        assert(parent(self.numer[3])==parent(tc_x.numer[3]))
        ratio=self.numer[0]/tc_x.numer[0]
        return all(self.numer[i]==ratio*tc_x.numer[i] for i in range(0,4))
            
    def Common_denom(self,tc_2):
        d=self.denom*tc_2.denom
        new_tc_1=Coord([self.numer[i]*tc_2.denom for i in range(0,4)],d)
        new_tc_2=Coord([tc_2.numer[i]*self.denom for i in range(0,4)],d)
        return new_tc_1,new_tc_2

    def Product(self,tc_2):
        prod_tc=Coord([self.numer[i]*tc_2.numer[i] for i in range(0,4)],self.denom*tc_2.denom)
        return prod_tc
           
    def Divide(self,tc_2): #self[i]/tc_2[i]
        assert(all(tc_2.numer[i]!=0 for i in range(0,4)))
        num_frac=[[self.numer[i],tc_2.numer[i]] for i in range(0,4)]
        num_frac,den=Common_denom_frac(num_frac)
        tc_1div2=Coord([num_frac[i][0]*tc_2.denom for i in range(0,4)],den*self.denom)
        return tc_1div2
     
    #term appeard in Riemann relation.
    def RR_product_term(self,tc_y,chi:int,i:int,j:int):
        sum=(self.field)(0)
        for t in range(0,4):
            chi_t=(self.field)((-1)**((chi//2)*(t//2)+(chi%2)*(t%2)))
            ipt=2*((i//2+t//2)%2)+(i%2+t%2)%2
            jpt=2*((j//2+t//2)%2)+(j%2+t%2)%2
            sum+=chi_t*self.numer[ipt]*tc_y.numer[jpt]
        return [sum,self.denom*tc_y.denom]
    
    def Reset_data(self):
        self.square=0
        self.hsq=0
        self.zichi=0
        self.lmd_lpow_value=0
        self.normalized_tc=0
        del self.square, self.hsq, self.zichi, self.lmd_lpow_value, self.normalized_tc

        

#-----------------------



class NullCoord(Coord):
    #for differential addition.
    def Kappa_1(self,tc_x:Coord,tc_y:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        #assert((tc_0_Hsq.numer)!=[0,0,0,0])
        tc_x_Hsq=tc_x.Tc_Hsq()
        tc_y_Hsq=tc_y.Tc_Hsq()
        prod_tc=tc_x_Hsq.Product(tc_y_Hsq)
        z0_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z0_chi.Hadamard()
        kappa_ii.denom*=4
        return kappa_ii

    #for doubling.
    def Kappa_3(self,tc_x:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        prod_tc=(tc_x.Tc_Hsq()).Square()
        z0_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z0_chi.Hadamard()
        kappa_ii.denom*=4
        return kappa_ii
    
    
    def Kappa_2(self,tc_x:Coord,tc_y:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        tc_x_Hsq=tc_x.Tc_Hsq()
        tc_y_Hsq=tc_y.Tc_Hsq()
        prod_tc=tc_x_Hsq.Product(tc_y_Hsq)
        z0_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z0_chi.Hadamard()
        kappa_ii.denom*=4
        zichi_0=self.Zichi()
        zichi_x=tc_x.Zichi()
        zichi_y=tc_y.Zichi()
        assert(zichi_0[key]!=0 for key in zichi_0)
        z_10_00=[zichi_x["_10^00"]*zichi_y["_10^00"],zichi_0["_10^00"]]
        z_10_01=[zichi_x["_10^01"]*zichi_y["_10^01"],zichi_0["_10^01"]]
        z_01_00=[zichi_x["_01^00"]*zichi_y["_01^00"],zichi_0["_01^00"]]
        z_01_10=[zichi_x["_01^10"]*zichi_y["_01^00"],zichi_0["_01^00"]]
        z_11_00=[zichi_x["_11^00"]*zichi_y["_11^00"],zichi_0["_11^00"]]
        z_11_11=[zichi_x["_11^11"]*zichi_y["_11^11"],zichi_0["_11^11"]]
        [z_10_00,z_10_01],den_10=Common_denom_frac([z_10_00,z_10_01])
        [z_01_00,z_01_10],den_01=Common_denom_frac([z_01_00,z_01_10])
        [z_11_00,z_11_11],den_11=Common_denom_frac([z_11_00,z_11_11])
        d_0sq=tc_0_Hsq.denom
        d_xsqysq=prod_tc.denom
        den_10*=d_xsqysq
        den_01*=d_xsqysq
        den_11*=d_xsqysq
        z_10_00=[z_10_00[0]*d_0sq,den_10]
        z_10_01=[z_10_01[0]*d_0sq,den_10]
        z_01_00=[z_01_00[0]*d_0sq,den_01]
        z_01_10=[z_01_10[0]*d_0sq,den_01]
        z_11_00=[z_11_00[0]*d_0sq,den_11]
        z_11_11=[z_11_11[0]*d_0sq,den_11]
        kappa_ij={(i,i):[kappa_ii.numer[i],kappa_ii.denom] for i in range(0,4)}
        kappa_ij[(0,2)]=[z_10_00[0]+z_10_01[0],4*den_10]
        kappa_ij[(1,3)]=[z_10_00[0]-z_10_01[0],4*den_10]
        kappa_ij[(0,1)]=[z_01_00[0]+z_01_10[0],4*den_01]
        kappa_ij[(2,3)]=[z_01_00[0]-z_01_10[0],4*den_01]
        kappa_ij[(0,3)]=[z_11_00[0]+z_11_11[0],4*den_11]
        kappa_ij[(1,2)]=[z_11_00[0]-z_11_11[0],4*den_11]
        return kappa_ij
    


    #[cf.LR16]section5.
    #from x,y,we compute set {x+y,x-y}. We can't distinguish x+y,x-y.
    def Kappa_2_2(self,tc_x,tc_y):
        lv2tnp=[self.numer[i]/self.denom for i in range(0,4)]
        lv2_x =[tc_x.numer[i]/tc_x.denom for i in range(0,4)]
        lv2_y =[tc_y.numer[i]/tc_y.denom for i in range(0,4)]    
        assert(not Is_Ellipitc_product(lv2tnp)[0])
        z_i_chi={i_chi:(product_term(lv2_x,i_chi[0],i_chi[1])*product_term(lv2_y,i_chi[0],i_chi[1])/product_term(lv2tnp,i_chi[0],i_chi[1])) for i_chi in even_theta_lv22}
        #kappa_{i,j}
        k_ij={(0,0):(z_i_chi[(0,0)]+z_i_chi[(0,1)]+z_i_chi[(0,2)]+z_i_chi[(0,3)]),
            (1,1):(z_i_chi[(0,0)]-z_i_chi[(0,1)]+z_i_chi[(0,2)]-z_i_chi[(0,3)]),
            (2,2):(z_i_chi[(0,0)]+z_i_chi[(0,1)]-z_i_chi[(0,2)]-z_i_chi[(0,3)]),
            (3,3):(z_i_chi[(0,0)]-z_i_chi[(0,1)]-z_i_chi[(0,2)]+z_i_chi[(0,3)]),
            (2,0):(z_i_chi[(2,0)]+z_i_chi[(2,1)]),
            (3,1):(z_i_chi[(2,0)]-z_i_chi[(2,1)]),
            (1,0):(z_i_chi[(1,0)]+z_i_chi[(1,2)]),
            (3,2):(z_i_chi[(1,0)]-z_i_chi[(1,2)]),
            (3,0):(z_i_chi[(3,0)]+z_i_chi[(3,3)]),
            (1,2):(z_i_chi[(3,0)]-z_i_chi[(3,3)])}
        for A in {(2,0),(3,1),(1,0),(3,2),(3,0),(1,2)}:
            k_ij[(A[1],A[0])]=k_ij[A]
        for key in k_ij:
            k_ij[key]=[k_ij[key],4]
        return k_ij




    def Normal_Add(self,tc_x:Coord,tc_y:Coord,number:int):
        assert(number==1 or number==2)
        #kappa_ij=self.Kappa_2(tc_x,tc_y)
        kappa_ij=self.Kappa_2(tc_x,tc_y)
        #print((kappa_ij[(0,2)][0]/kappa_ij[(0,2)][1])/k_ij[(0,2)])
        alpha=1
        
        root_in=(kappa_ij[(0,alpha)][0]/kappa_ij[(0,alpha)][1])**2-(kappa_ij[(0,alpha)][0]/kappa_ij[(0,alpha)][1])*(kappa_ij[(0,0)][0]/kappa_ij[(0,alpha)][1])
        print(100,sqrt(root_in))
        
        kappa_a0sq=[kappa_ij[(0,alpha)][0]**2,kappa_ij[(0,0)][1]**2]          
        kappa_aa00=[kappa_ij[(alpha,alpha)][0]*kappa_ij[(0,0)][0],kappa_ij[(alpha,alpha)][1]*kappa_ij[(0,0)][1]]
        D=Frac_minus(kappa_a0sq,kappa_aa00)
        print(sqrt(D[0]/D[1]))
        root_D=[sqrt((D[0]*D[1])),D[1]]
        k00_rootD=[kappa_ij[(0,0)][i]*root_D[i] for i in range(0,2)]
        X0=kappa_ij[(0,0)]
        X1=Frac_sum(kappa_ij[(0,alpha)],root_D)
        X2_num=Frac_minus(Frac_mult(X1,kappa_ij[(0,2)]),Frac_mult(kappa_ij[(0,0)],kappa_ij[(alpha,2)]))
        X2=Frac_div(X2_num,k00_rootD)
        X3_num=Frac_minus(Frac_mult(X1,kappa_ij[(0,3)]),Frac_mult(kappa_ij[(0,0)],kappa_ij[(alpha,3)]))
        X3=Frac_div(X3_num,k00_rootD)
        [X0,X1,X2,X3],den=Common_denom_frac([X0,X1,X2,X3])
        tc_xpy=Coord([X0[0],X1[0],X2[0],X3[0]],den)
        if number==1:
            return tc_xpy
        Y0=[self.field(1),self.field(1)]
        ka0_min_rootD=Frac_minus(kappa_ij[(0,alpha)],root_D)
        Y1_num=Frac_minus(Frac_mult(ka0_min_rootD,kappa_ij[(0,1)]),Frac_mult(kappa_ij[(0,0)],kappa_ij[(alpha,1)]))
        Y1=Frac_div(Y1_num,k00_rootD)
        Y2_num=Frac_minus(Frac_mult(ka0_min_rootD,kappa_ij[(0,2)]),Frac_mult(kappa_ij[(0,0)],kappa_ij[(alpha,2)]))
        Y2=Frac_div(Y2_num,k00_rootD)
        Y3_num=Frac_minus(Frac_mult(ka0_min_rootD,kappa_ij[(0,3)]),Frac_mult(kappa_ij[(0,0)],kappa_ij[(alpha,3)]))
        Y3=Frac_div(Y3_num,k00_rootD)
        [Y0,Y1,Y2,Y3],den_y=Common_denom_frac([Y0,Y1,Y2,Y3])
        tc_xmy=Coord([Y0[0],Y1[0],Y2[0],Y3[0]],den_y)
        return tc_xpy,tc_xmy





    #calculate x+z from x,z,x+y,y+z.
    def Compatible_Add(self,tc_x:Coord,tc_z:Coord,tc_xpy:Coord,tc_ypz:Coord):
        tc_xpz_1,tc_xpz_2=self.Normal_Add(tc_x,tc_z,2)
        tc_1,tc_2=self.Normal_Add(tc_xpy,tc_xpz_1,2)
        tc_3,tc_4=self.Normal_Add(tc_xpy,tc_xpz_2,2)
        if   tc_ypz.Is_same_proj(tc_1) or tc_ypz.Is_same_proj(tc_2):
            return tc_xpz_2
        elif tc_ypz.Is_same_proj(tc_3) or tc_ypz.Is_same_proj(tc_4):
            return tc_xpz_1
        else :
            assert(False)

    def Double(self,tc_x:Coord):
        kappa_ii=self.Kappa_3(tc_x)
        tc_2x=kappa_ii.Divide(self)
        return tc_2x
        
    def Diff_Add(self,tc_x:Coord,tc_y:Coord,tc_xmy:Coord):
        if (self==tc_xmy):
            return self.Double(tc_x)
        kappa_ii=self.Kappa_1(tc_x,tc_y)
        tc_xpy=kappa_ii.Divide(tc_xmy)
        return tc_xpy


    #calculate kx+y from x,y,x+y.
    def Kxpy_xpy(self,k:int,tc_x:Coord,tc_y:Coord,tc_xpy:Coord):
        if k==0:
            return tc_y
        elif k==1:
            return tc_xpy
        elif k==2:
            return self.Diff_Add(tc_xpy,tc_x,tc_y) #2x+y=(x+y)+x
        else:
            bit_km1=(ZZ(k-1)).digits(2)
            n=len(bit_km1) #bit lingth of k-1.
            assert(bit_km1[n-1]==1)
            X=tc_x
            Y=tc_xpy
            Z=tc_y
            for i in range(0,n-1):
                if bit_km1[i]==1:
                    dX=self.Double(X)
                    XpY=self.Diff_Add(X,Y,Z)
                    X=dX
                    Y=XpY
                else:
                    assert(bit_km1[i]==0)
                    dX=self.Double(X)
                    XpZ=self.Diff_Add(X,Z,Y)
                    X=dX
                    Z=XpZ
            result=self.Diff_Add(X,Y,Z)
            return result



    #calculate kx+y from x,y,x+y for plural k.
    def Plural_kxpy_xpy(self,coeff:list,tc_x:Coord,tc_y:Coord,tc_xpy:Coord):
        k=max(coeff)
        bit_km1=(ZZ(k-1)).digits(2)
        n=len(bit_km1) #bit lingth of k-1
        assert(bit_km1[n-1]==1)
        #X[i]=(2^i)*tc_x.
        X=[tc_x]
        for i in range(0,n-1):
            X.append(self.Double(X[i]))
        assert(len(X)==n)
        bit_coeffm1=[(ZZ(coeff[t]-1)).digits(2) for t in range(0,len(coeff))]
        kpxy_list=[]
        
        for t in range(0,len(coeff)):
            Y=tc_xpy
            Z=tc_y
            for i in range(0,len(bit_coeffm1[t])):
                if bit_coeffm1[t][i]==1:
                    Y=self.Diff_Add(X[i],Y,Z)
                else:
                    Z=self.Diff_Add(X[i],Z,Y)
            kpxy_list.append(Y)
        assert(len(kpxy_list)==len(coeff))
        return kpxy_list
        
        
        
    def Mult(self,tc_x:Coord,k:int):
        return self.Kxpy_xpy(k,tc_x,self,tc_x)
    
    
    def Plural_Mult(self,tc_x:Coord,coeff:list):
        return self.Plural_kxpy_xpy(coeff,tc_x,self,tc_x)
    
    
        
    def Is_order(self,tc_x:Coord,k:int):
        assert(k>=1)
        if k==1:
            return self.Is_same_proj(tc_x)
        else:
            if self.Is_same_proj(tc_x):
                print("it is 0")
                return False
            for i in factor(k):
                s=k//(i[0])
                if self.Is_same_proj(self.Mult(tc_x,s)):
                    return False
            if self.Is_same_proj(self.Mult(tc_x,k)):
                return True
            print("the order is over.")
            return False

    
        
    #[numer.denom]
    def Lmd_lpow(self,tc_e:Coord):
        try:
            return tc_e.lmd_lpow_value
        except AttributeError:
            l=tc_e.order
            assert(l!=0)
            ld=(l-1)//2
            [num,den]=self.Plural_Mult(tc_e,[ld,ld+1])
            tc_e.lmd_lpow_value=[num.numer[0]*den.denom,den.numer[0]*num.denom]
            return tc_e.lmd_lpow_value

    def Is_excellent(self,tc_e:Coord):
        assert(tc_e.order!=0)
        return self.Is_same_affine(self.Mult(tc_e,tc_e.order))

   
    #[numer.denom]
    def MudivLm_lpow(self,tc_e:Coord,tc_x:Coord,tc_epx:Coord,lm_l_sq):
        l=tc_e.order
        assert(l!=0)
        lepx=self.Kxpy_xpy(l,tc_e,tc_x,tc_epx)
        assert(type(lepx)==Coord)
        mudivlm_lpow=[tc_x.numer[0]*(lm_l_sq[1])*lepx.denom,tc_x.denom*(lm_l_sq[0])*lepx.numer[0]]
        return mudivlm_lpow


    #Three way addition.
    def Extended_Addition(self,tc_x:Coord,tc_y:Coord,tc_z:Coord,tc_xpy:Coord,tc_ypz:Coord,tc_zpx:Coord):
        prod_0_ypz  =[self.numer  [i]*tc_ypz.numer[i] for i in range(0,4)]
        prod_zpx_xpy=[tc_zpx.numer[i]*tc_xpy.numer[i] for i in range(0,4)]
        prod_y_z    =[tc_y.numer  [i]*tc_z.numer  [i] for i in range(0,4)]
        sum_prod_0_ypz  =[0,0,0,0]
        sum_prod_zpx_xpy=[0,0,0,0]
        sum_prod_y_z    =[0,0,0,0]
        for chi in range(0,4):
            chi_t=[((-1)**((chi//2)*(t//2)+(chi%2)*(t%2))) for t in range(0,4)]
            sum_prod_0_ypz  [chi]=sum([chi_t[t]*prod_0_ypz  [t] for t in range(0,4)])
            sum_prod_zpx_xpy[chi]=sum([chi_t[t]*prod_zpx_xpy[t] for t in range(0,4)])
            sum_prod_y_z    [chi]=sum([chi_t[t]*prod_y_z    [t] for t in range(0,4)])        
        E_chi=[[(sum_prod_0_ypz[chi]*sum_prod_zpx_xpy[chi]),sum_prod_y_z[chi]] for chi in range(0,4)]
        E_chi,den1=Common_denom_frac(E_chi)
        pre_tc_xpypz=[[0,1],[0,1],[0,1],[0,1]]
        for i in range(0,4):
            chi_i=[((-1)**((chi//2)*(i//2)+(chi%2)*(i%2))) for chi in range(0,4)]
            sum_part_i=sum([chi_i[chi]*E_chi[chi][0] for chi in range(0,4)])
            pre_tc_xpypz[i]=[sum_part_i,4*tc_x.numer[i]]
        pre_tc_xpypz,den_xpypz=Common_denom_frac(pre_tc_xpypz)
        den_xpypz*=tc_xpy.denom*tc_ypz.denom*tc_zpx.denom*den1
        dx_dy_dz=tc_x.denom*tc_y.denom*tc_z.denom*self.denom
        for i in range(0,4):
            pre_tc_xpypz[i][0]*=dx_dy_dz
        tc_xpypz=Coord([pre_tc_xpypz[i][0] for i in range(0,4)],den_xpypz)
        return tc_xpypz
    

