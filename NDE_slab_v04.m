%-------------------------------%
% Slab Reactor in One-Dimension
%-------------------------------%
clear; clc;

%___________%
%Given data
%-----------%
%Diffusion Coefficient
D=0.90; %cm
%Absorption Coefficient
Sig_Abs=0.066; %cm^-1
%Fission Coefficient
Nu_Sgf=0.07;


%Number of mesh points
N=100;
%Length vector for which criticality to be determined
L=1:0.1:50;
l=length(L);
%K value vector to store each K value against L value
K_vec=zeros(1,l);

%Loop for each length value
for m=1:l

%delta x => Mesh size
dx=L(m)/(N-1);

% %___________________________________________________________%
%Finite Difference matrix coefficients of Diffusion Equation
aL=zeros(1,N);
aC=zeros(1,N);
aR=zeros(1,N);
aD=zeros(1,N);

for i=2:N-1
  aL(i)=-D/dx;
  aR(i)=-D/dx;
  aC(i)=Sig_Abs*dx-aL(i)-aR(i);
  aD(i)=Nu_Sgf*dx;
end

%Left Side boundary condition (Reflective)
Lb=0.5;
aR(1)=-D/dx;
aC(1)=Sig_Abs*dx+D*Lb/(1+dx*Lb/2)-aR(1);
aD(1)=Nu_Sgf*dx;
%Right side boundary condition (Vacuum)
Rb=0.5;
aL(N)=-D/dx;
aC(N)=Sig_Abs*dx+D*Rb/(1+dx*Rb/2)-aL(N);
aD(N)=Nu_Sgf*dx;

%Initial guesses for K and Source
K=1.0;
Source(1:N)=1.0;

%Tolerence for K
K_tol=1e-4;
K_inst_tol=999;

%_____________%
%Inverse Power method
%-------------%

while K_inst_tol>K_tol
    %Old value update
    K_old=K;
    S_old=Source;
    B=Source/K;
%Calling fuction of Successive over-relaxation for Flux solution
    Flux_calc=SOR(aL,aC,aR,B,N);

    %Source matrix update
    Source(1:N)=0.0;
    Source=aD.*Flux_calc;
    %K value update
    K=K_old*sum(Source)/sum(S_old);
    %K tolerence calculation
    K_inst_tol=abs(K-K_old)/K_old;
end

%K value matrix
K_vec(m)=K;
end

%Critical value from criticality matrix
[value,index]=min(abs(K_vec-1));
K_value=K_vec(index);
%Printing Crtical value of length
fprintf('Length for reactor for which k=1 is %4.2f\n',L(index))


%___________________________________%
%Plot of K value vs Lenngth of the Slab reactor
%---------------------------------------------%
plot(L,K_vec)
xlabel ("Length (cm)"); ylabel ("K");
title('Criticality value vs Length of 1-D Slab reactor');
hold on
plot(L(index),K_vec(index),'.g', 'MarkerSize',20)
% hold on
% 
% dim = [0.75 0.7 0.09 0.17];
% annotation('textbox',dim,'String','K=1','FitBoxToText','on');
hold off


%____________________________________________%
%Function of Successive over relaxation method
%--------------------------------------------%
function Flux_out=SOR(aL,aC,aR,B,N)
%Flux initial guess
Flux=ones(1,N);
%Relaxation parameter
w=1.4;
%Flux convergence tolerence
Flux_tol=1E-03;
Flux_diff=999;

while Flux_diff>Flux_tol
 
 Flux_old=Flux;
 Flux(1)=(B(1)-aR(1)*Flux_old(2))*w/aC(1)+(1-w)*Flux_old(1);
for i=2:N-1
    Flux(i)=(B(i)-aL(i)*Flux(i-1)-aR(i)*Flux_old(i+1))*w/aC(i)+(1-w)*Flux_old(i);
end
Flux(N)=(B(N)-aL(N)*Flux(N-1))*w/aC(N)+(1-w)*Flux_old(N);

Flux_diff=abs(sum(Flux)-sum(Flux_old));
end
Flux_out=Flux;
end