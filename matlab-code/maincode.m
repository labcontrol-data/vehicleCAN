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

Ac =[-0.0716    4.4868         0;
   -6.8028   -5.0147    6.3637;
  -20.2376  -18.0076  -45.3408];
Bc =[0; 0; 11.6605];

K = [-1.3 -0.1 -0.2 0.5];
Ac = [Ac zeros(3,1);
      -1 0 0 -1];
Bc = [Bc; 0];

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

N=max(size(posy{1}))
t = TS*[0:N-1];

figure(2)
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identification procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load('dataIdentificatinThrottle_Labcontrol_HassaneAissam12-Nov-2015_18h35m')
load('dataIdentificatinThrottle_Labcontrol_HassaneAissam20-Feb-2016_15h8m')

posy=10*simout.signals.values(200:39500,2);
vely=10*simout.signals.values(200:39500,3);
current=simout.signals.values(200:39500,1);
vm = simout.signals.values(200:39500,4);

gainGradeToRad = 0.01745329251994;
posy = gainGradeToRad*posy';
vely= gainGradeToRad*vely';

Nit = max(size(posy));

t=[0:TS:TS*Nit-TS]';

for k=1:max(size(posy))-1
    x{k}=[posy(k);vely(k);current(k)];
end

dataTrajectory = [];
for k=1:max(size(posy))-1
    dataTrajectory = [dataTrajectory; [vely(k) ]];
end

Y1=[]; U1=[];
Y2=[]; U2=[];
Y3=[]; U3=[];
for k=1:max(size(posy))-2
    Y1=[Y1; x{k+1}(1)];
    U1=[U1; x{k}(1) x{k}(2)];
    
    Y2=[Y2; x{k+1}(2)];
    U2=[U2; x{k}(1) x{k}(2) x{k}(3)];
    
    Y3=[Y3; x{k+1}(3)];
    U3=[U3; x{k}(2) x{k}(3) vm(k)];
end
beta1=inv(U1'*U1)*U1'*Y1;
beta2=inv(U2'*U2)*U2'*Y2;
beta3=inv(U3'*U3)*U3'*Y3;

A=[beta1(1) beta1(2) 0;
   beta2(1) beta2(2) beta2(3);
   0   beta3(1)  beta3(2)];
B=[0;0;beta3(3)];
C=eye(3,3);
D=0*B;
K=0*C;
x0=x{1};

%A(1,1)=1;

init_sys=idss(A,B,C,D,K,x0,TS);

data=iddata([posy;vely;current']',vm, TS);
data.InputName={'tensao'};
data.OutputName={'pos-rad','vel-rads','current'};

%figure(1),
%compare(data,init_sys)

init_sys.Structure.a.Free(3,1)=false;
init_sys.Structure.b.Free(1,1)=false;
init_sys.Structure.b.Free(2,1)=false;
init_sys.Structure.c.Free=false;

[modelo,x0]=ssest(data,init_sys,'TS',TS);

%figure(3),
%compare(data,modelo)

%sysc = d2c(modelo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparacao dados apresentacao artigo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=10*simout.signals.values(:,2);
w=10*simout.signals.values(:,3);
c=simout.signals.values(:,1);
u = simout.signals.values(:,4);

gainGradeToRad = 0.01745329251994;
theta = gainGradeToRad*theta';
w= gainGradeToRad*w';

data=iddata([theta;w;c']',u, TS);
data.InputName={'tensao'};
data.OutputName={'pos-rad','vel-rads','current'};

figure(11),
compare(data,init_sys)


[Y,fit,x0]=compare(data,modelo);
thetaSim=Y.OutputData(:,1);
wSim=Y.OutputData(:,2);
cSim=Y.OutputData(:,3);

n=max(size(A));


savefile = 'identificationThrottleVehiclePlatoon.mat';
save(savefile, 'theta', 'w','c','u','thetaSim', 'wSim','cSim','t','-v7');



clear A, clear B,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for Identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

posy=10*simout.signals.values(8000:end,2);
vely=10*simout.signals.values(8000:end,3);
current=simout.signals.values(8000:end,1);
vm = simout.signals.values(8000:end,4);

gainGradeToRad = 0.01745329251994;
posy = gainGradeToRad*posy';
vely= gainGradeToRad*vely';

Nit = max(size(posy));

t=[0:TS:TS*Nit-TS]';

for k=1:max(size(posy))-1
    x{k}=[posy(k);vely(k);current(k)];
end

dataTrajectory = [];
for k=1:max(size(posy))-1
    dataTrajectory = [dataTrajectory; [vely(k) ]];
end

 
px=1;
cluster_number = px; % Set yout cluster number
[idx,centrTraject,sumDist] = kmeans(dataTrajectory, cluster_number,'MaxIter',500); % idx is the index array, for each sample data
sprintf(' cluster_number = %d,  sumaDistancia=%f',px,sum(sumDist))

Y1=[]; U1=[];
Y2=[]; U2=[];
Y3=[]; U3=[];
for j=1:px
    Y1{j}=[]; U1{j}=[]; Y2{j}=[]; U2{j}=[]; Y3{j}=[]; U3{j}=[];
end


for k=1:max(size(posy))-2
    m = idx(k);
    
    Y1{m}=[Y1{m}; x{k+1}(1)];
    U1{m}=[U1{m}; x{k}(1) x{k}(2)];
    
    Y2{m}=[Y2{m}; x{k+1}(2)];
    U2{m}=[U2{m}; x{k}(1) x{k}(2) x{k}(3)];
    
    Y3{m}=[Y3{m}; x{k+1}(3)];
    U3{m}=[U3{m}; x{k}(1) x{k}(2) x{k}(3) vm(k)];
end

for j=1:px
    beta1=inv(U1{j}'*U1{j})*U1{j}'*Y1{j};
    beta2=inv(U2{j}'*U2{j})*U2{j}'*Y2{j};
    beta3=inv(U3{j}'*U3{j})*U3{j}'*Y3{j};
    
    A{j}=[beta1(1) beta1(2) 0;
        beta2(1) beta2(2) beta2(3);
        beta3(1) beta3(2) beta3(3)];
    B{j}=[0;0;beta3(4)];
end

x=x{1};
vecx1=[x(1)];
vecx2=[x(2)];
vecx3=[x(3)];

for k=1:max(size(posy))-2
    m = idx(k);
    x = A{m}*x + B{m}*vm(k);
    vecx1=[vecx1 x(1)];
    vecx2=[vecx2 x(2)];
    vecx3=[vecx3 x(3)];
end

figure(10)
subplot(3,1,1)
hold on
plot(t,posy,'k')
plot(t(1:max(size(vecx1))),vecx1,'b')
hold off
legend('real','simul.')
ylabel('position')

subplot(3,1,2)
hold on
plot(t,vely,'k')
plot(t(1:max(size(vecx2))),vecx2,'b')
axis([0 80 -1 1])
hold off
legend('real','simul.')
ylabel('velocity')

subplot(3,1,3)
hold on
plot(t,current,'k')
plot(t(1:max(size(vecx3))),vecx3,'b')
axis([0 80 -3 3])
hold off
legend('real','simul.')
ylabel('current'),xlabel('time (s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert from discrete- to continuous-time model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=max(size(A{1}));
for i=1:px
    Ad{i}=A{i};
    Bd{i}=B{i};
    A{i} = (Ad{i} - eye(n,n))*inv(TS);
    B{i} = Bd{i}*inv(TS);
end

G = [-2 -0.78 -0.1 7.4];
for i=1:px 
    Aext{i} = [A{i} zeros(3,1);
              -1  zeros(1,3)];
    Bext{i} = [B{i}; 0];
    Amf{i} = Aext{i}+Bext{i}*G;
end

Acontinuous = A{1}
Bcontinuous = B{1}

A=Ad{1};
B=Bd{1};
C=eye(3,3);
D=0*B;
K=0*C;
x0=[posy(1);vely(1);current(1)];
modelo2=idss(A,B,C,D,K,x0,TS);

figure(14),
compare(data,modelo)

