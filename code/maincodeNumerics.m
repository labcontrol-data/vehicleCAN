% paper: ``Stability of asynchronous sampled-data systems with 
%          input delay: application to an automotive throttle valve''
% Matlab 2017a
% Written by Alessandro N. Vargas
% email: vargas.alessandro@gmail.com
% Last updated: Abr 21, 2021

clear all, clc, close all, format short, format long,


%##############################################
%#  paper's numerical example
%##############################################
Ac=[-2 1;
   0  -0.9];
Bc=[1 0]';
K = [1 -1];

H=Ac+Bc*K;

n=max(size(Ac));

if max(real(eig(Ac+0.5*eye(n))))<0
    disp('The pair (Ac,Bc) is stabilizable')
else
    disp('!!! The pair (Ac,Bc) is NOT stabilizable !!!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute maximum delay according to the condition in the paper [*]
% [*] F. Mazenc and D. Normand-Cyrot, IEEE Trans Aut Control,   
%     Reduction Model Approach for Linear Systems With Sampled 
%     Delayed Inputs,  2013, vol. 58,  p. 1263-1268
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this command calculates the matrix Qr that satisfies
% (A+BK)'*Qr + Qr*(A+BK) == -I < 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=max(size(Ac));
Q = lyap(H', eye(n));
residual_value = norm( H'*Q + Q*H +  eye(n) );

tau=1;

DeltaM_old = inv(4*sqrt(6)*norm(Bc)*norm(K))...
              *min( [inv(norm(H)*norm(Q*expm(Ac*tau)))  inv(norm(expm(Ac*tau)))    ] )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New method described in the novel paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = expm(Ac*tau)*Bc*K;

beta = 0.56;

val1 = (1 - 1/(2*beta))*inv(6*norm(Q*expm(Ac*tau)*Bc*K)^2*norm(H)^2) ;

a = 2*norm(expm(Ac*tau)*Bc*K)^2;
b = 2*beta/3;
c = -1;
rootsVal = [(-b+sqrt(b^2 - 4*a*c))/(2*a)  (-b-sqrt(b^2 - 4*a*c))/(2*a)];
val2 = max(rootsVal);

DeltaM_new = min(val1,val2)

delta = DeltaM_new ;

a1 = -1 + 1/(2*beta) + (3*delta*norm(Q*expm(Ac*tau)*Bc*K)^2*norm(H)^2)
a2 = -1 + 2*beta*delta/3 + delta^2*norm(expm(Ac*tau)*Bc*K)^2

if (a1<0)&&(a2<0)
    sprintf('Note that both a1 and a2 are negative under delta = %.8f',delta)
end
