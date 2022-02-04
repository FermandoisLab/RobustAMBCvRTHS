%% Adaptive model-based compensator for RTHS benchmark with non-linear experimental substructure

%Simulation parameters
fs=4096;           %Sampling frequency   
dtsim=1/fs;        %Time step

% Sensors:
%Displacement sensor
rms_noise_desp = 6.25e-14;    %  Noise power
rmsdesp=sqrt(rms_noise_desp/dtsim)*1000;  %rms mm

%Force sensor
rms_noise_F = 1.16e-3; % Noise power
rmsF=sqrt(rms_noise_F/dtsim);  %rms N

% Saturation limits 
sat_limit_upper = +3.8;         % Volts
sat_limit_lower = -3.8;         % Volts
% Quantization interval
quantize_int = 1 / 2^18;        % 18 bit channel

%% Transfer system dynamics and Initial Model
%Transfer system nominal parameters
gain=1;      
a1_b0=2.13e13;      
a2=4.23e6;           
a3=3.3;
b1=425;
b2=1e5;
% mass, damping and stiffness for model without specimen
meLoadSys=0;
keLoadSys=1; 
ceLoadSys=0;

%Load system's transfer function
s=tf('s');
Gest = tf(1,[meLoadSys ceLoadSys keLoadSys]);
Ga = tf(1,[1 a3]);
Gs = tf(a1_b0,[1 b1 b2]);
Gcsi = a2*s;
G0 = feedback(Ga*Gest,Gcsi,-1);
GpLoadSys = gain*feedback(Gs*G0,1,-1);

%Third order model
orden=3;
[num,den]=tfdata(GpLoadSys,'v');
num=num(end);
den=den(end-orden:end);
GpLoadSysred=tf(num,den);    %Initial model
AMB=tf(den/num,1);           
A_amb_i=flip(den/num)        %Initial compensator's parameters a_i

%% Compensator

% Backward difference
u=[1 0 0 0];
dudt=[1 -1 0 0]/dtsim;
dudt2=[1 -2 1 0]/dtsim^2;
dudt3=[1 -3 3 -1]/dtsim^3;
Derivadas=[u;dudt;dudt2;dudt3];  %FIR filter for numerical derivatives

%Adaptive parameters constraints
%a0
amb0max=inf;   
amb0min=0;
%a1
amb1max=inf;
amb1min=0;
%a2
amb2max=inf;
amb2min=0;
%a3
amb3max=inf;
amb3min=0;
%For saturation block
max_amb=[amb0max,amb1max,amb2max,amb3max];
min_amb=[amb0min,amb1min,amb2min,amb3min];

%Filter for adaptation
n = 4;   %order
fc = 20; %cutoff frequency
fs=1/dtsim;
[numfilter,denfilter] = butter(n,fc/(fs/2)); %discrete filter coefficients

% Adaptive gains
adaptivegain=diag(10.^[8.4 6.2 2.1 0.8]);   %obtained from calibration

%% Ground motion

load('KobeAccelNoScaling.mat');
load('ElCentroAccelNoScaling.mat');

escala=50/100;                %scaling factor

%Earthquakes:
Ug=ElCentroAccel';
earthquake= ['El Centro'];

% Ug=KobeAccel';
% earthquake= ['Kobe'];

dtUg=Ug(2,1)-Ug(1,1);  %acceleration record time step
Ugdesde=0;             %initial time for simulation
Ughasta=Ug(end,1);            %final time for simulation
Ug=[(0:dtUg:(Ughasta-Ugdesde))',Ug(Ugdesde/dtUg+1:Ughasta/dtUg+1,2)];
Ug(:,2)=Ug(:,2)*escala;  %scaled acceleration

ttotal=max(Ug(:,1));  %total time of simulation

%% Reference Structure
m=1000;  %mass per floor kg
d=0.05;  %damping ratio per mode

M=eye(3)*m;  %Mass matrix
K=10^7*[2.6055,-2.3134,0.5937;-2.3134,3.2561,-1.4420;0.5937,-1.4420,0.9267];  %Stiffness matrix N/m


[mod,w2]=eig(K,M); %modes and natural frequencies
w=sqrt(w2);        %natural frecuencies (rad/sec)
f=w/2/pi*ones(3,1) %natural frecuencies (Hz)

C=mod'^-1*(2*d*w)*mod^-1;  %Damping matrix

%% Partitioning and experimental substructure model

%Experimental substructure properties
me=29.1;   %mass kg
ce=114.6;  %viscous damping
ke=1.19*10^6;  %Initial stiffness           

%Bouc Wen parameters
asiv=0.1;        % post-yield stiffness ratio [0,1]  ; 0: elastoplastic ; 1: linear
xy=1/1000;  % yield displacement
n1siv=0.1;     % hysteresis shape parameters[0,1]
n2siv=1-n1siv; % 1-n1siv
nsiv=2.5;      % smoothness [1,inf] ; 1 to 10 smooth ; >20 bi-linear ; >> unstable

%degradation
stiffdeg=0.02*0;    % stiffness degradation ; 0: no degradation
strdeg=0.005*0;     % strenght degradation ; 0: no degradation

%cubic hardening spring
khard=6e9*0;     % Fhard=khard*x^3

% state-space model
ssAe=[0,1;-asiv*ke/me,-ce/me];
ssBe=[0,0;1/me,-1/me];
ssCe=[1,0;0,1;asiv*ke,0];
ssDe=[0,0;0,0;0,1];

%Numerical substructure
rest=zeros(3);
rest(1,1)=1;

Mn=M-me*rest;
Cn=C-ce*rest;
Kn=K-ke*rest;

G1=ones(3,1);
G2=[1,0,0]';

ssAn=[zeros(3),eye(3);-Mn\Kn,-Mn\Cn];
ssBn=[zeros(3,2);-Mn\M*G1,-Mn\G2];
ssCn=[eye(3),zeros(3)];
ssDn=zeros(3,2);

%Reference structure
Kmod=Kn;

Ar=[zeros(3),eye(3);-M\Kmod -M\C];
Br=[zeros(3,2);-ones(3,1),-M\G2];
Cr=[eye(6)];
Dr=[zeros(6,2)];

%% Simulation

sim('AMB_NL_2_simulation.slx')

ts=x_m.Time;         %simulation time
x_t=x_t.Data(:,1);   %target displacement
x_m=x_m.Data(:,1);   %measured displacement
x_r1=x_r.Data(:,1);  %reference displacement
Ugs=Ugs.Data(:,1);   %ground motion

F_m=F_m.Data(:,1);   %measured force (total with noise)
rt=rt.Data(:,1);     %experimental force (only restitutive r(t))
F_r=F_r.Data(:,1);   %experimental force in reference structure
Kcur=Kcur.Data(:,1); %experimental stiffness
Fy=Fy.Data(:,1);     %experimental yield force

%Adaptive parameters
amb0=amb.Data(:,1);
amb1=amb.Data(:,2);
amb2=amb.Data(:,3);
amb3=amb.Data(:,4);

%Evaluation
J2=rms(x_t-x_m)/rms(x_t)*100;
J4=rms(x_r1-x_m)/rms(x_r1)*100;
[Amptotal,phitotal,feqtotal,delaytotal] = Freq_Resp_Tong(x_t,x_m,1/dtsim);
%results table
Ji = string({'J2 [%]';'delay [ms]';'J4 [%]'});
Ji = cellstr(Ji);
results=[J2;delaytotal*1000;J4];
table(results,'VariableNames',{'Results'},'RowNames',Ji)
%% Displacement Plots
%measured and reference
figure
subplot(3,1,1)
plot(ts,x_r1*1000,'b')
hold on
plot(ts,x_m*1000,'r--')
legend('x_{r1}','x_m','Orientation','Horizontal')
xlabel('Time [sec]')
ylabel('Disp. [mm]')
grid on

subplot(3,1,2)
plot(ts,abs(x_r1*1000-x_m*1000),'k')
hold on
legend(['NRMSE = ',num2str(J4),' %'])%  ; Delay = ',num2str(delaytotal*1000), ' ms']) 
xlabel('Time [sec]')
ylabel('|error| [mm]')
grid on

subplot(3,1,3)
plot(ts,Ugs/9.81,'b')
legend([earthquake,' ',num2str(escala*100), '%']) 
xlabel('Time [sec]')
ylabel('Earthquake [g]')
grid on

%synchronization
figure
subplot(3,1,1)
plot(ts,x_t*1000,'k')
hold on
plot(ts,x_m*1000,'r--')
legend('x_t','x_m','Orientation','Horizontal','Location','best')
xlabel('Time [sec]')
ylabel('Disp. [mm]')
grid on

subplot(3,1,2)
plot(ts,abs(x_t*1000-x_m*1000),'k')
hold on
legend(['J_2 = ',num2str(J2),' %'])%  ; Delay = ',num2str(delaytotal*1000), ' ms']) 
xlabel('Time [sec]')
ylabel('|error| [mm]')
grid on

subplot(3,1,3)
plot(ts,Ugs/9.81,'b')
legend([earthquake,' ',num2str(escala*100), '%']) 
xlabel('Time [sec]')
ylabel('Earthquake [g]')
grid on

%Zoom
desde=round(9/dtsim);
hasta=round(10/dtsim);
figure
plot(ts(desde:hasta),x_t(desde:hasta)*1000,'k','LineWidth',1)
hold on
plot(ts(desde:hasta),x_m(desde:hasta)*1000,'r--','LineWidth',2)
legend('x_t','x_m','Orientation','Horizontal','Location','best')
xlabel('Time [sec]')
ylabel('Displacement [mm]')
xlim([desde hasta]*dtsim)
grid on

%% Force plots

figure
plot(x_m*1000,rt,'b')
xlabel('Disp [mm]')
ylabel('r(t) [N]')
title('Restitutive force')
grid on

figure
plot(x_r1*1000,F_r,'b')
hold on
plot(x_m*1000,F_m,'r')
legend('Reference','vRTHS','Location','best')
xlabel('Measured Disp. [mm]')
ylabel('Measured Force [N]')
grid on

figure
plot(ts,Kcur/ke,'b')
hold on
plot(ts,Fy/ke/xy,'k--')
legend('K_{cur}/k_{e}','Fy_{cur}/Fy_{i}')
xlabel('Time [sec]')
ylabel('Stiffness/strength degradation')
grid on

%% Adaptive parameters
figure
subplot(2,2,1)
plot(ts,amb0,'k')
xlabel('Time [sec]')
ylabel('a_0')
grid on

subplot(2,2,2)
plot(ts,amb1,'k')
xlabel('Time [sec]')
ylabel('a_1  [sec]')
grid on

subplot(2,2,3)
plot(ts,amb2,'k')
xlabel('Time [sec]')
ylabel('a_2  [sec^2]')
grid on

subplot(2,2,4)
plot(ts,amb3,'k')
xlabel('Time [sec]')
ylabel('a_3   [sec^3]')
grid on

%Identified parameters
A_amb_final=mean([amb0 amb1 amb2 amb3]);

%bode
Gest = tf(1,[me ce ke]);
G0 = feedback(Ga*Gest,Gcsi,-1);
Gp = gain*feedback(Gs*G0,1,-1);

fbode=0.1:0.1:100;
wbode=2*pi*fbode;
[magGp,phaseGp,wout]=bode(Gp,wbode); %Experimental initial bode
[magID,phaseID,wout]=bode(tf(1,flip(A_amb_final)),wbode); %Identified bode
[maginicial,phaseinicial,wout]=bode(tf(1,flip(A_amb_i)),wbode);  %Initial controller

figure
subplot(2,1,1)
plot(fbode,db(squeeze(maginicial)),'r')
hold on
plot(fbode,db(squeeze(magGp)),'b--')
plot(fbode,db(squeeze(magID)),'g--')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on
legend('Initial (compensator)','Exp. sub. (initial)','Identified')
subplot(2,1,2)
plot(fbode,(squeeze(phaseinicial)*pi/180)./wbode'*1000,'r')
hold on
plot(fbode,(squeeze(phaseGp)*pi/180)./wbode'*1000,'b--')
plot(fbode,(squeeze(phaseID)*pi/180)./wbode'*1000,'g--')
xlabel('Frequency [Hz]')
ylabel('Time Delay [msec]')
grid on
