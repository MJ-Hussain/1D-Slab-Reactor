%________________1D 2 Group Reactor_______________%
clc; clear;
%Given data
%Cross-sections of given materials 
   CXs=[1/(3*0.2375) 0.002 0 0.0293         %water G1
       1/(3*1.2467) 0.0286 0 0              %water G2
       1.5118 9.0075e-3 6.8547e-3 1.4282e-2 %Assembley Ao G1
       0.403 8.0567e-2 1.3903e-1 0          %Assembley Ao G2
       1.5094 8.7859e-3 6.4038e-3 1.4522e-2 %Assembley Bo G1
       0.4033 7.4984e-2 1.2704e-1 0         %Assembley Bo G2
       1.4965 7.9306e-03 4.6046e-3 1.5703e-02 %Assembley Co G1
       0.40253 5.1306e-2 7.602e-2 0];         %Assembley Co G2
%_____________________________________________%

n_reg=4;  %number of regions

len=20; %cm %region length

dx=0.5; %mesh size

N=n_reg*len/dx; %number of mesh points

ng=2;  %number of groups
%____________________________________________%
nm=0;
st=0;
sp=0;
for i=1:n_reg
    st=sp+1;
    sp=sp+len/dx;
    for g=1:ng
    nm=nm+1;
    D(g,st:sp)=CXs(nm,1);
    
    Sig_a(g,st:sp)=CXs(nm,2);
    
    Sig_f(g,st:sp)=CXs(nm,3);
    
    Sig_s12(g,st:sp)=CXs(nm,4);
    
    end
end
%___________________________________________%

%Tri_diagonal matrix elements of Coefficient matrix
Lhs=zeros(ng,N);
Cs=zeros(ng,N);
Rhs=zeros(ng,N);
nF=zeros(ng,N);
S_12=zeros(ng,N);
for g=1:ng
    for i=2:N-1
       Lhs(g,i)= -2*D(g,i-1)*D(g,i)/(D(g,i)*dx+D(g,i-1)*dx);
       Rhs(g,i)= -2*D(g,i)*D(g,i+1)/(D(g,i+1)*dx+D(g,i)*dx);
       Cs(g,i)=(Sig_a(g,i)+Sig_s12(g,i))*dx-Lhs(g,i)-Rhs(g,i);
       nF(g,i)=Sig_f(g,i)*dx;
       S_12(g,i)=Sig_s12(g,i)*dx;
    end
    
    %Left side Boundary
    gL=1E5;
    Rhs(g,1)= -2*D(g,1)*D(g,2)/(D(g,2)*dx+D(g,1)*dx);
    Cs(g,1)=(Sig_a(g,1)+Sig_s12(g,1))*dx-Rhs(g,1)+gL*D(g,1)/(D(g,1)+gL*0.5*dx);
    nF(g,1)=Sig_f(g,1)*dx;
    S_12(g,1)=Sig_s12(g,1)*dx;
    
    %Right side Boundary
    gR=0;
    Lhs(g,N)=-2*D(g,N-1)*D(g,N)/(D(g,N)*dx+D(g,N-1)*dx);
    Cs(g,N)=(Sig_a(g,N)+Sig_s12(g,N))*dx-Lhs(g,N)+gR*D(g,N)/(D(g,N)+gR*0.5*dx);
    nF(g,N)=Sig_f(g,N)*dx;
    S_12(g,N)=Sig_s12(g,N)*dx;
end

%_________________________________________________________________%
K_tolerence=1e-06;
Flux_tolerence=1e-08;

%Call for solution of K and Flux using Power Method
[Flux,K]=PM(Rhs,Cs,Lhs,nF,S_12,K_tolerence,Flux_tolerence,ng,N);
% %________________________________________________________%

%Flux Normalization
mx=max(Flux(1,:));
Norm_F1=Flux(1,:)/mx*1E12;
Norm_F2=Flux(2,:)/mx*1E12;

%Power peaking factor
%Fast flux 
PPf=max(Flux(1,:))/mean(Flux(1,:));
PPt=max(Flux(2,:))/mean(Flux(2,:));

%Print of solution
fprintf('The max K value is %6.4f\n',K)
fprintf('Power peaking Factor for Fast Flux is %4.2f\nPower Peaking factor for thermal flux is %4.2f\n',PPf,PPt)

%Plot of Fast and thermal flux
X=linspace(0,n_reg*len,N);
plot(X,Norm_F1,X,Norm_F2)
xlabel('Length (cm)'); ylabel('Normalized Flux')
legend('Fast Flux', 'Thermal Flux')


%______________________________________________________________________%
%                      Power method
%______________________________________________________________________%
function [Phi,K]=PM(R,C,L,F,S12,K_tol,Flux_tol,ng,N)
%Initial guesses
Phi=ones(ng,N);
K=1.0;
Src(1:N)=(F(1,:).*Phi(1,:)+F(2,:).*Phi(2,:));

%Tolerence
%K_tol=1e-6;
K_diff=1E3;
Src_tol=1e-04;
Src_diff=1E3;

%Number of iterations
iter=0;
    while (K_diff>K_tol) && (Src_diff>Src_tol)
        
        iter=iter+1;
        K_prev=K;
        S_prev=Src;

        for g=1:ng
            %Fast group Flux
           if g==1
              S=(F(1,:).*Phi(1,:)+F(2,:).*Phi(2,:))./K;
              %Thermal Group Flux
           else
               S=S12(1,:).*Phi(1,:);
           end
            
            %Flux_tol=1E-08; %tolerence
            Flux_diff=1E3;
            Flux_prev=ones(1,N);
            w=1.5;  %Relaxation parameter

            while Flux_diff>Flux_tol
                 
                 Phi(g,1)=(S(1)-R(g,1)*Flux_prev(2))*w/C(g,1)+(1-w)*Flux_prev(1);
                 for i=2:N-1
                    Phi(g,i)=(S(i)-L(g,i)*Phi(g,i-1)-R(g,i)*Flux_prev(i+1))*w/C(g,i)+(1-w)*Flux_prev(i);
                 end
                 Phi(g,N)=(S(N)-L(g,N)*Phi(g,N-1))*w/C(g,N)+(1-w)*Flux_prev(N);
                 Flux_diff=abs(sum(Phi(g,:))-sum(Flux_prev));
                 Flux_prev=Phi(g,:);
            end
             
        end
        Src=0.0;
        Src=F(1,:).*Phi(1,:)+F(2,:).*Phi(2,:);
        K=K_prev*sum(Src)/sum(S_prev);
        
        Src_diff=abs(sum(Src)-sum(S_prev))/sum(S_prev);
        K_diff=abs(K-K_prev)/K_prev;
  
    end
end
