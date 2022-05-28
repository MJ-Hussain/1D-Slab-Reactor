####################################################################
###                  1_Dimensional Slab reactor                  ###
####################################################################
'''
Netron Diffusion equation is given by
-d/dx D(x).d/dx phi(x) + Sigma_Abs(x).phi(x) = 1/K.Nu_Sigma_f.phi(x)
where 
D          => Neutron Diffusion Coefficient
Sigma_Abs  => Neutron Absorption Coefficient
Nu_Sigma_f => Neutron fission Coefficient

Finite Difference discretization results a tridiagonal System 
=>     aL(i) phi(i-1) + aC(i) phi(i) + aR(i) phi(i+1) = 1/K aD(i) phi(i)

K is criticality eigen value for which Power Method is employed
With initial guess of K and phi the resulting linear system has form
=>          Ax=B
Which is solved using Successive Over-relaxation method in each outer iteration of Power Method

'''
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------#
#Successive over-relaxation method for Solution of tridiagonal system
#------------------------------------------------------------------------#
def SOR(a,b,c,B,N):
    phi_prev=np.zeros(N)
    phi=np.zeros(N)
    #initial guess
    phi_prev[:]=1.0
    phi[:]=1.0
    #Over-relaxation parameter
    w=1.4
    #Tolerence for covergence 
    tol=1E-03
    #Tolerence in each iteration to be compared
    diff=999.0
    
    while diff>tol:
        
        phi_prev=phi.copy()
        
        phi[0]=(B[0]-c[0]*phi_prev[1])*w/b[0]+(1-w)*phi_prev[0]
        for i in range(1,N-1):
            phi[i]=(B[i]-a[i]*phi[i-1]-c[i]*phi_prev[i+1])*w/b[i]+(1-w)*phi_prev[i]
        phi[N-1]=(B[N-1]-a[N-1]*phi[N-2])*w/b[N-1]+(1-w)*phi_prev[N-1]
        diff=abs(sum(phi)-sum(phi_prev))/sum(phi_prev)
    return phi 
#--------------------------------------------------------------------------------------#

##--------##
#Given data
##--------##
#Netron Diffusion Coefficient
D1=0.9   
#Netron Absorption Coefficient  
SigA1=0.066
#Netron Fission Coefficient 
NSgf1=0.07 


#Reactor length vector
#For each length corresponding K value to be determined
L=np.arange(30,60,1)
length=len(L)
#.
K_vec=np.zeros(length)
#Number of Mesh points
N=50
#Initialize Coefficient Matrix Tridiagonls
aL=np.zeros(N)
aC=np.zeros(N)
aR=np.zeros(N)
d=np.zeros(N)

#Initial Guess for K
K=1.0

#Loop for each L value
for m in range(length):
 #Mesh size
    dx=L[m]/(N-1)
    
    #Coefficient Matrix Tridiagonls Calculation
    for i in range(1,N-1):
        aL[i]=-D1/dx
        aR[i]=-D1/dx
        aC[i]=SigA1*dx-aR[i]-aL[i]
        d[i]=NSgf1*dx
    
    #Left boundary condition reflective
    rL=0.0
    aR[0]=-D1/dx
    aC[0]=SigA1*dx-aR[0]+rL*D1/(D1+dx*rL/2)
    d[0]=NSgf1*dx
    #Right boundary condition Vacuum
    rR=0.5
    aL[N-1]=-D1/dx
    aC[N-1]=SigA1*dx-aL[N-1]+D1*rR/(D1+dx*rR/2)
    d[N-1]=NSgf1*dx
    
##-------------------##    
#     Power Method    #
##-------------------##
    #  Source term initial guess => Nu_Sigma_f*phi
    Source=np.ones(N)
    #To store previous source value
    S_old=np.zeros(N)
    #Left hand side of System Ax=B 
    B=np.zeros(N)
    
    #K tolerence
    K_tol=1E-04
    #K tolerence in each iteration initialize
    K_diff=999
    #Store previous value of K
    K_prev=0.0
    
    
    while K_diff>K_tol:
        K_prev=K
        S_old=Source.copy()
        B=Source/K
        #Call for Successive Over-relaxation method to solve Ax=B
        Flux=SOR(a=aL,b=aC,c=aR,B=B,N=N)
        
        Source[:]=0.0
        #update Source term
        Source=d*Flux
        #Update K value
        K=K_prev*sum(Source)/sum(S_old)
        #K tolerence calculation
        K_diff=abs(K-K_prev)/K_prev
#-----Power method loop ends here------------#
    #Store Converged K value  
    K_vec[m]=K    
#----L value loop ends here--------------#

#Finding K value near to 1 => criticality
idx = (np.abs(K_vec-1)).argmin()
print('The Length of Reactor at which it becomes critical is ',"{0:.2f}".format(L[idx]))

#-----------------------------#
#           Plots             #
#-----------------------------#
plt.figure(1)
plt.plot(L,K_vec)
plt.plot(L[idx],K_vec[idx],'ro',label='K=1')
plt.xlabel('Length (cm)')
plt.ylabel('Criticality value')
plt.title('Length vs criticality of a 1-Dimensional Slab Reactor')
plt.legend()
plt.figure(2)
x=np.linspace(0,30,N)
plt.plot(x,Flux)
plt.xlabel('Dimension in X (a.u.)')
plt.ylabel('Flux (a.u.)')
plt.title('Flux Profile in 1-Dimensional Slab Reactor')
plt.show()

