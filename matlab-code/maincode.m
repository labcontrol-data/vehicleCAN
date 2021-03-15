% this code shows the identified mode, simulation, and graphics
% Matlab 2017a
% Written by Alessandro N. Vargas
% email: vargas.alessandro@gmail.com
% Last updated: March 12, 2021

clear all, clc, close all, format short, format long,


%##############################################
%#  simple numerical example
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
Qr = lyap(H', eye(n));
residual_value = norm( H'*Qr + Qr*H +  eye(n) );

tau=1;

DeltaM_old = inv(4*sqrt(6)*norm(Bc)*norm(K))...
              *min( [inv(norm(H)*norm(Qr*expm(Ac*tau)))  inv(norm(expm(Ac*tau)))    ] )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New method described in the novel paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = 1;

val1 = (1 - 1/(2*beta))*inv(6*norm(Qr*expm(Ac*tau)*Bc*K)^2*norm(H)^2);

a = 2*norm(expm(Ac*tau)*Bc*K)^2;
b = 2*beta/3;
c = -1;
rootsVal = [(-b+sqrt(b^2 - 4*a*c))/(2*a)  (-b-sqrt(b^2 - 4*a*c))/(2*a)];
val2 = max(rootsVal);
DeltaM_new = min(val1,val2)

varEPS = 1e-6;
delta = -varEPS + DeltaM_new ;

a1 = -1 + 1/(2*beta) + (6*delta*norm(Qr*expm(Ac*tau)*Bc*K)^2*norm(H)^2)
a2 = -1 + 2*beta*delta/3 + 2*delta^2*norm(expm(Ac*tau)*Bc*K)^2

if (a1<0)&&(a2<0)
    sprintf('Note that both a1 and a2 are negative under delta = %.8f',delta)
end
 

format short,
 

%##############################################
% code that process experimental data
%##############################################


%##############################################
%#  code to check the theoretical condition
%##############################################

vectau=[]; vecDeltaM=[];

% Ac =[-0.0421    4.3891         0;
%    -8.6703   -5.7236    8.2085;
%   -24.2994  -42.7586  -48.5546];
% Bc =[0; 0;12.7132];

Ac =[-0.0716    4.4868         0;
   -6.8028   -5.0147    6.3637;
  -20.2376  -18.0076  -45.3408];
Bc =[0; 0; 11.6605];

% Ac =[-0.0741    4.4873         0;
%      -6.9432   -4.9586    6.6801;
%      -25.6087  -20.1949  -63.4326];
% Bc =[ 0; 0; 15.8059];

%K = [-1.33 -0.1 -0.2];
%Ac=[-2 1;
%    0  -0.9];
%Bc=[1 0]';
%K = [1 -1];

% K = [-2 -0.78 -0.1 7.4];
K = [-1.3 -0.1 -0.2 0.5];
Ac = [Ac zeros(3,1);
      -1 0 0 -1];
Bc = [Bc; 0];

% 
% A=Amf;
% B=Bext;

H=Ac+Bc*K;


n=max(size(Ac));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this command calculates the matrix Qr that satisfies
% (A+BK)'*Qr + Qr*(A+BK) == -I < 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qr = lyap(H',eye(n));

N=3600;    % number of point for each second interval
TS = 1/N;  % sampling time

beta = 2;  %old value = 1.2

delta = TS;
ExpectedValue=0;

lambMinQr = min(eig(Qr));
lambMaxQr = max(eig(Qr));
for i=1:N
    tau = TS*i;
    %a1 = -1 + 1/(2*beta) + (6*norm(sqrt(Qr)*expm(Ac*tau)*Bc*K)^2*delta*norm(H)^2);
    %a2 = -1 + 2*beta*delta*inv(3*lambMinQr) + 2*delta^2*norm(expm(Ac*tau)*Bc*K)^2;
    a1 = -1 + 1/(2*beta) + 6*delta*norm(Qr*expm(Ac*tau)*Bc*K)^2*norm(H)^2;
    a2 = -1 + 2*beta*delta/3 + 2*delta^2*norm(expm(Ac*tau)*Bc*K)^2;
    vecAlpha(i) = max(a1,a2); 
    vecTau(i) = tau;
end


figure(1)
plot(vecTau,vecAlpha)
xlabel('tau')
ylabel('alpha(tau)')
grid




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experimental data under control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('dataThrottle_Labcontrol_VehiclePlatoonControl11-Mar-2016_19h34m')
posy{1}=10*simout.signals.values(2950:8000,2);
vely{1}=10*simout.signals.values(2950:8000,3);
current{1}=simout.signals.values(2950:8000,1);
vm{1} = simout.signals.values(2950:8000,5);
reference{1} = 2.4*simout.signals.values(2950:8000,4);

gainGradeToRad = 0.01745329251994;
posy{1} = gainGradeToRad*posy{1}';
vely{1}= gainGradeToRad*vely{1}';


load('dataThrottle_Labcontrol_VehiclePlatoonControl11-Mar-2016_19h44m')
posy{2}=10*simout.signals.values(2950:8000,2);
vely{2}=10*simout.signals.values(2950:8000,3);
current{2}=simout.signals.values(2950:8000,1);
vm{2} = simout.signals.values(2950:8000,5);
reference{2} = 2.4*simout.signals.values(2950:8000,4);

posy{2} = gainGradeToRad*posy{2}';
vely{2}= gainGradeToRad*vely{2}';

%  load('dataThrottle_Labcontrol_VehiclePlatoonControl12-Mar-2016_1h18m26s')
%  posy{1}=10*simout.signals.values(2000:8000,2);
%  vely{1}=10*simout.signals.values(2000:8000,3);
%  current{1}=simout.signals.values(2000:8000,1);
%  vm{1} = simout.signals.values(2000:8000,5);

N=max(size(posy{1}))
t = TS*[0:N-1];

figure(21)
subplot(5,1,1)
hold on
plot(t,posy{1},'b')
plot(t,posy{2},'k')
hold off
grid
ylabel('position(rad)')

subplot(5,1,2)
hold on
plot(t,vely{1},'b')
plot(t,vely{2},'k')
hold off
grid
ylabel('velocity(rad/sec)')

subplot(5,1,3)
hold on
plot(t,current{1},'b')
plot(t,current{2},'k')
hold off
grid
ylabel('current(A)')

subplot(5,1,4)
hold on
plot(t,vm{1},'b')
plot(t,vm{2},'k')
hold off
grid
ylabel('control(V)')

subplot(5,1,5)
plot(t,reference{1},'b')
ylabel('reference(V)')

xlabel('time(sec)')
legend('delay=10msec','delay=180msec')


theta=[posy{1}; posy{2}];
w=[vely{1}; vely{2}];
c=[current{1} current{2}]';
u = [vm{1} vm{2}]';
ref = [reference{1} reference{2}]';

savefile = 'controlDataThrottleVehiclePlatoon.mat';
save(savefile, 'theta', 'w','c','u','t','ref','vecTau','vecAlpha','-v7');
