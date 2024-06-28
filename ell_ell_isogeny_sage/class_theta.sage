

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
    

    def Sum_chit_thetaipt_thetat(self):
        try:
            return self.sctt_ichi
        except AttributeError:
            prod_00_01=self.numer[0]*self.numer[1]
            prod_00_10=self.numer[0]*self.numer[2]
            prod_00_11=self.numer[0]*self.numer[3]
            prod_01_10=self.numer[1]*self.numer[2]
            prod_01_11=self.numer[1]*self.numer[3]
            prod_10_11=self.numer[2]*self.numer[3]
            self.sctt_ichi={}
            self.sctt_ichi[(1,0)]=2*prod_00_01+2*prod_10_11
            self.sctt_ichi[(1,2)]=2*prod_00_01-2*prod_10_11
            self.sctt_ichi[(2,0)]=2*prod_00_10+2*prod_01_11
            self.sctt_ichi[(2,1)]=2*prod_00_10-2*prod_01_11
            self.sctt_ichi[(3,0)]=2*prod_00_11+2*prod_01_10
            self.sctt_ichi[(3,3)]=2*prod_00_11-2*prod_01_10
            return self.sctt_ichi


    def Mult_const(self,lmd):
        return Coord([self.numer[i]*lmd for i in range(0,4)],self.denom)

    def Mult_frac(self,frac:list):
        return Coord([self.numer[i]*frac[0] for i in range(0,4)],self.denom*frac[1])
 
    def Is_same_affine(self,tc_x):
        return all(self.numer[i]*tc_x.denom==tc_x.numer[i]*self.denom for i in range(0,4))
    
    def Is_same_proj(self,tc_x):
        K=self.field
        ratio=self.numer[0]/tc_x.numer[0]
        return all(self.numer[i]==ratio*tc_x.numer[i] for i in range(0,4))
            
    def Common_denom(self,tc_2):
        d=self.denom*tc_2.denom
        new_tc_1=Coord([self.numer[i]*tc_2.denom for i in range(0,4)],d)
        new_tc_2=Coord([tc_2.numer[i]*self.denom for i in range(0,4)],d)
        return new_tc_1,new_tc_2

    def Multiply(self,tc_2):
        prod_tc=Coord([self.numer[i]*tc_2.numer[i] for i in range(0,4)],self.denom*tc_2.denom)
        return prod_tc
           
    def Divide(self,tc_2): #self[i]/tc_2[i]
        assert(all(tc_2.numer[i]!=0 for i in range(0,4)))
        num_frac=[[self.numer[i],tc_2.numer[i]] for i in range(0,4)]
        num_frac,den=Common_denom_frac(num_frac)
        tc_1div2=Coord([num_frac[i][0]*tc_2.denom for i in range(0,4)],den*self.denom)
        return tc_1div2
         
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
    def Kappa_ii(self,tc_x:Coord,tc_y:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        #assert((tc_0_Hsq.numer)!=[0,0,0,0])
        tc_x_Hsq=tc_x.Tc_Hsq()
        tc_y_Hsq=tc_y.Tc_Hsq()
        prod_tc=tc_x_Hsq.Multiply(tc_y_Hsq)
        z0_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z0_chi.Hadamard()
        kappa_ii.denom*=4
        return kappa_ii

    #for doubling.
    def Kappa_double(self,tc_x:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        prod_tc=(tc_x.Tc_Hsq()).Square()
        z0_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z0_chi.Hadamard()
        kappa_ii.denom*=4
        return kappa_ii
    
    

    def Kappa_ij(self,tc_x:Coord,tc_y:Coord):
        tc_0_Hsq=self.Tc_Hsq()
        tc_x_Hsq=tc_x.Tc_Hsq()
        tc_y_Hsq=tc_y.Tc_Hsq()
        prod_tc=tc_x_Hsq.Multiply(tc_y_Hsq)
        z_00_chi=prod_tc.Divide(tc_0_Hsq)
        kappa_ii=z_00_chi.Hadamard()
        kappa_ii.denom*=4
        sctt_ichi_0=self.Sum_chit_thetaipt_thetat()
        sctt_ichi_x=tc_x.Sum_chit_thetaipt_thetat()
        sctt_ichi_y=tc_y.Sum_chit_thetaipt_thetat()
        z_10_00=[sctt_ichi_x[(2,0)]*sctt_ichi_y[(2,0)],sctt_ichi_0[(2,0)]]
        z_10_01=[sctt_ichi_x[(2,1)]*sctt_ichi_y[(2,1)],sctt_ichi_0[(2,1)]]
        z_01_00=[sctt_ichi_x[(1,0)]*sctt_ichi_y[(1,0)],sctt_ichi_0[(1,0)]]
        z_01_10=[sctt_ichi_x[(1,2)]*sctt_ichi_y[(1,2)],sctt_ichi_0[(1,2)]]
        z_11_00=[sctt_ichi_x[(3,0)]*sctt_ichi_y[(3,0)],sctt_ichi_0[(3,0)]]
        z_11_11=[sctt_ichi_x[(3,3)]*sctt_ichi_y[(3,3)],sctt_ichi_0[(3,3)]]
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
    




    def Normal_Add(self,tc_x:Coord,tc_y:Coord,number:int):
        assert(number==1 or number==2)
        kappa_ij=self.Kappa_ij(tc_x,tc_y)
        alpha=1
        kappa_a0sq=[kappa_ij[(0,alpha)][0]**2,kappa_ij[(0,alpha)][1]**2]          
        kappa_aa00=[kappa_ij[(alpha,alpha)][0]*kappa_ij[(0,0)][0],kappa_ij[(alpha,alpha)][1]*kappa_ij[(0,0)][1]]
        D=Frac_sub(kappa_a0sq,kappa_aa00)
        root_D=[sqrt(D[0]*D[1]),D[1]]
        X1_num=Frac_add(kappa_ij[(0,alpha)],root_D)
        X0=[1,1]
        X1=Frac_div(X1_num,kappa_ij[(0,0)])
        X2_num=Frac_sub(Frac_mult(X1,kappa_ij[(0,2)]),kappa_ij[(1,2)])
        X3_num=Frac_sub(Frac_mult(X1,kappa_ij[(0,3)]),kappa_ij[(1,3)])
        X2=Frac_div(X2_num,root_D)
        X3=Frac_div(X3_num,root_D)
        [X0,X1,X2,X3],den_X=Common_denom_frac([X0,X1,X2,X3])
        tc_xpy=Coord([X0[0],X1[0],X2[0],X3[0]],den_X)
        if number==1:
            return tc_xpy
        Y0=Frac_div(kappa_ij[0,0],X0)
        Y1=Frac_div(kappa_ij[1,1],X1)
        Y2=Frac_div(kappa_ij[2,2],X2)
        Y3=Frac_div(kappa_ij[3,3],X3)
        [Y0,Y1,Y2,Y3],den_Y=Common_denom_frac([Y0,Y1,Y2,Y3])
        tc_xmy=Coord([Y0[0],Y1[0],Y2[0],Y3[0]],den_Y)
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
        else:
            assert(False)


    def Double(self,tc_x:Coord):
        kappa_ii=self.Kappa_double(tc_x)
        tc_2x=kappa_ii.Divide(self)
        return tc_2x
        
    def Diff_Add(self,tc_x:Coord,tc_y:Coord,tc_xmy:Coord):
        if (self==tc_xmy):
            return self.Double(tc_x)
        kappa_ii=self.Kappa_ii(tc_x,tc_y)
        tc_xpy=kappa_ii.Divide(tc_xmy)
        return tc_xpy



    #calculate kx+y from x,y,x+y.
    def Kxpy_xpy(self,k:int,tc_x:Coord,tc_y:Coord,tc_xpy:Coord):
        assert(k>=0)
        if k==0 or k==1 or k==2 or k==3 or k==5:
            tc_0xpy=tc_y
            if k==0:
                return tc_0xpy
            tc_1xpy=tc_xpy
            if k==1:
                return tc_1xpy
            tc_2xpy=self.Diff_Add(tc_1xpy,tc_x,tc_0xpy)
            if k==2:
                return tc_2xpy
            tc_3xpy=self.Diff_Add(tc_2xpy,tc_x,tc_1xpy)
            if k==3:
                return tc_3xpy
            tc_4xpy=self.Diff_Add(tc_3xpy,tc_x,tc_2xpy)
            tc_5xpy=self.Diff_Add(tc_4xpy,tc_x,tc_3xpy)
            if k==5:
                return tc_5xpy
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



    def Mult(self,tc_x:Coord,k:int):
        assert(k>=0)
        if k<=5:
            tc_0x=self
            if k==0:
                return tc_0x
            tc_1x=tc_x
            if k==1:
                return tc_1x
            tc_2x=self.Double(tc_1x)
            if k==2:
                return tc_2x
            tc_3x=self.Diff_Add(tc_2x,tc_x,tc_x)
            if k==3:
                return tc_3x
            tc_4x=self.Diff_Add(tc_3x,tc_x,tc_2x)
            if k==4:
                return tc_4x
            tc_5x=self.Diff_Add(tc_4x,tc_x,tc_3x)
            if k==5:
                return tc_5x
        else:
            return self.Kxpy_xpy(k,tc_x,self,tc_x)
    
    
    
    
    
    def Mult_Mult(self,tc_x:Coord,max_au:int):
        assert(max_au>=1)
        kx_list=[self,tc_x]
        if max_au==1:
            return kx_list
        else:
            for k in range(2,max_au+1):
                kx=self.Diff_Add(kx_list[k-1],tc_x,kx_list[k-2])
                kx_list.append(kx)
            assert(len(kx_list)==max_au+1)
            return kx_list
                

    
        
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

        

    def Is_excellent(self,tc_e:Coord,l:int):
        return self.Is_same_affine(self.Mult(tc_e,l))
    

    #Three way addition.
    def Extended_Addition(self,tc_x:Coord,tc_y:Coord,tc_z:Coord,tc_xpy:Coord,tc_ypz:Coord,tc_zpx:Coord):
        prod_0_ypz  =[self.numer  [i]*tc_ypz.numer[i] for i in range(0,4)]
        prod_zpx_xpy=[tc_zpx.numer[i]*tc_xpy.numer[i] for i in range(0,4)]
        prod_y_z    =[tc_y.numer  [i]*tc_z.numer  [i] for i in range(0,4)]
        K=self.field
        sum_prod_0_ypz  =[K(0),K(0),K(0),K(0)]
        sum_prod_zpx_xpy=[K(0),K(0),K(0),K(0)]
        sum_prod_y_z    =[K(0),K(0),K(0),K(0)]
        for chi in range(0,4):
            chi_t=[K((-1)**((chi//2)*(t//2)+(chi%2)*(t%2))) for t in range(0,4)]
            for t in range(0,4):
                sum_prod_0_ypz  [chi]+=chi_t[t]*prod_0_ypz  [t]
                sum_prod_zpx_xpy[chi]+=chi_t[t]*prod_zpx_xpy[t]
                sum_prod_y_z    [chi]+=chi_t[t]*prod_y_z    [t]   
        E_chi=[[(sum_prod_0_ypz[chi]*sum_prod_zpx_xpy[chi]),sum_prod_y_z[chi]] for chi in range(0,4)]
        E_chi,den1=Common_denom_frac(E_chi)
        pre_tc_xpypz=[[0,1],[0,1],[0,1],[0,1]]
        for i in range(0,4):
            chi_i=[K((-1)**((chi//2)*(i//2)+(chi%2)*(i%2))) for chi in range(0,4)]
            sum_part_i=K(0)
            for chi in range(0,4):
                sum_part_i+=K(chi_i[chi]*E_chi[chi][0])
            pre_tc_xpypz[i]=[sum_part_i,4*tc_x.numer[i]]
        pre_tc_xpypz,den_xpypz=Common_denom_frac(pre_tc_xpypz)
        den_xpypz*=tc_xpy.denom*tc_ypz.denom*tc_zpx.denom*den1
        dx_dy_dz=tc_x.denom*tc_y.denom*tc_z.denom*self.denom
        for i in range(0,4):
            pre_tc_xpypz[i][0]*=dx_dy_dz
        tc_xpypz=Coord([pre_tc_xpypz[i][0] for i in range(0,4)],den_xpypz)
        return tc_xpypz
    

