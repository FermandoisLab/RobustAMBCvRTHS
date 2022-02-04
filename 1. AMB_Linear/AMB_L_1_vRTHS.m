%% RTHS benchmark with Adaptive model-based compensator 

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

%Noise filter for adaptation process
n = 4;   %order
fc = 20; %cutoff frequency
fs=1/dtsim;
[numfilter,denfilter] = butter(n,fc/(fs/2)); %discrete filter coefficients

% Adaptive gains
adaptivegain=diag(10.^[8.4 6.2 2.1 0.8]);   %obtained from calibration

%% Ground motion

load('KobeAccelNoScaling.mat');
load('ElCentroAccelNoScaling.mat');
load('MauleAccelNoScaling.mat');

escala=50/100;                %scaling factor

%Earthquakes:
Ug=ElCentroAccel';
earthquake= ['El Centro'];

% Ug=KobeAccel';
% earthquake= ['Kobe'];

% Ug=MauleAccel';
% earthquake= ['Maule'];

dtUg=Ug(2,1)-Ug(1,1);  %acceleration record time step
Ugdesde=0;             %initial time for simulation
Ughasta=Ug(end,1);            %final time for simulation
Ug=[(0:dtUg:(Ughasta-Ugdesde))',Ug(Ugdesde/dtUg+1:Ughasta/dtUg+1,2)];
Ug(:,2)=Ug(:,2)*escala;  %scaled acceleration

ttotal=max(Ug(:,1));  %total time of simulation

%% Reference structure

m=1000;  %mass per floor kg
d=0.05;  %damping ratio per mode

M=eye(3)*m;  %Mass matrix
K=10^7*[2.6055,-2.3134,0.5937;-2.3134,3.2561,-1.4420;0.5937,-1.4420,0.9267];  %Stiffness matrix N/m

[mod,w2]=eig(K,M); %modes and natural frequencies
w=sqrt(w2);        %natural frecuencies (rad/sec)
f=w/2/pi*ones(3,1) %natural frecuencies (Hz)

C=mod'^-1*(2*d*w)*mod^-1;  %Damping matrix

%state-space model
Ar=[zeros(3),eye(3);-M\K -M\C];
Br=[zeros(3,1);-ones(3,1)];
Cr=eye(6);
Dr=zeros(6,1);
REFSYS=ss(Ar,Br,Cr,Dr);

%% Simulations

nsimu=10;   %number of realizations
inc=1;     %inc=0 for nominal case   inc=1 for simulations with uncertainties

%matrices for results and simulated plants
J2results=zeros(1,nsimu);
delayresults=zeros(1,nsimu);
J4results=zeros(1,nsimu);
plants={};

%Bode for initial model and simulated plants
fbode=0.1:0.1:100;  %bode frecuency range
wbode=2*pi*fbode;
[maginicial,phaseinicial,wout]=bode(tf(1,flip(A_amb_i)),wbode);

figure(20)
subplot(2,1,1)
plot(fbode,db(squeeze(maginicial)),'r')
hold on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on
subplot(2,1,2)
plot(fbode,(squeeze(phaseinicial)*pi/180)./wbode'*1000,'r')
hold on
xlabel('Frequency [Hz]')
ylabel('Time Delay [msec]')
grid on

%Simulations
for simu=1:nsimu  
   
disp(simu);
%transfer system uncertainty
a3=3.3+1.3*randn(1,1)*inc;
b1=425+3.3*randn(1,1)*inc;
b2=1e5+3.3e3*randn(1,1)*inc;

%Experimental substructure
%benchmark exp. sub.
me=29.1;
ce=114.6;
ke=1.19*10^6+5e4*randn(1)*inc;

%random ex. sub.
% me=10+(40-10)*rand(1);
% ce=50+(1e3-50)*rand(1);
% ke=1e4+(2e6-1e4)*rand(1);

%Experimental substructure state-space model
ssAe=[0,1;-ke/me,-ce/me];
ssBe=[0;1/me];
ssCe=[1,0;0,1];
ssDe=[0;0];

%Numerical substructure
rest=zeros(3);
rest(1,1)=1;
Mn=M-me*rest;
Cn=C-ce*rest;
Kn=K-ke*rest;
G1=ones(3,1);
G2=[1,0,0]';

%Numerical substructure state-space model
ssAn=[zeros(3),eye(3);-Mn\Kn,-Mn\Cn];
ssBn=[zeros(3,2);-Mn\M*G1,-Mn\G2];
ssCn=[eye(3),zeros(3)];
ssDn=zeros(3,2);
NUMSYS=ss(ssAn,ssBn,ssCn,ssDn);

%Simulation
sim('AMB_L_8_vRTHS_simulation.slx')

%response
ts=x_m.Time;         %time
x_t=x_t.Data(:,1);   %target displacement
x_m=x_m.Data(:,1);   %measured displacement
x_r1=x_r.Data(:,1);  %reference structure response
x_c=x_c.Data(:,1);   %commanded displacement
Ugs=Ugs.Data(:,1);   %ground motion

%Adaptive parameters
amb0=amb.Data(:,1); %a_0
amb1=amb.Data(:,2); %a_1
amb2=amb.Data(:,3); %a_2
amb3=amb.Data(:,4); %a_3

%Error indicators
J2=rms(x_t-x_m)/rms(x_t)*100;  %NRMSE target-measured
J4=rms(x_r1-x_m)/rms(x_r1)*100;  %NRMSE reference-measured
[Amptotal,phitotal,feqtotal,delaytotal] = Freq_Resp_Tong(x_t,x_m,1/dtsim); %FEI indicator

%results
J2results(simu)=J2;
delayresults(simu)=delaytotal;
J4results(simu)=J4;

%simulated plant bode
Ga = tf(1,[1 a3]);
Gs = tf(a1_b0,[1 b1 b2]);
Gcsi = a2*s;
Gest = tf(1,[me ce ke]);
G0 = feedback(Ga*Gest,Gcsi,-1);
Gp = gain*feedback(Gs*G0,1,-1);
plants{simu}=Gp;

[magGp,phaseGp,wout]=bode(Gp,wbode);

figure(20)
subplot(2,1,1)
plot(fbode,db(squeeze(magGp)),'b--')
subplot(2,1,2)
plot(fbode,(squeeze(phaseGp)*pi/180)./wbode'*1000,'b--')


end
subplot(2,1,1)
legend('Initial','Simulated')

%Table with results
Ji = string({'J2 [%]';'delay [ms]';'J4 [%]'});
Ji = cellstr(Ji);
results=[J2results;delayresults*1000;J4results];
table(mean(results,2),std(results,0,2),max(results,[],2),min(results,[],2),'VariableNames',{'Mean','std','Max','Min'},'RowNames',Ji)
%% Time history of last simulation
%reference / measured / earthquake
figure
subplot(3,1,1)
plot(ts,x_r1*1000,'b')
hold on
plot(ts,x_m*1000,'r--')
legend('x_{r1}','x_m','Orientation','Horizontal','Location','best')
xlabel('Time [sec]')
ylabel('Disp. [mm]')
xlim([0 ttotal])
grid on

subplot(3,1,2)
plot(ts,abs(x_r1*1000-x_m*1000),'k')
hold on
legend(['NRMSE = ',num2str(J4),' %'])
xlabel('Time [sec]')
ylabel('|error| [mm]')
xlim([0 ttotal])
grid on

subplot(3,1,3)
plot(ts,Ugs/9.81,'b')
legend([earthquake,' ',num2str(escala*100), '%'],'Location','best')
xlabel('Time [sec]')
ylabel('Earthquake [g]')
xlim([0 ttotal])
grid on

%target / measured
figure
subplot(3,1,1)
plot(ts,x_t*1000,'k')
hold on
plot(ts,x_m*1000,'r--')
legend('x_t','x_m','Orientation','Horizontal','Location','best')
xlabel('Time [sec]')
ylabel('Disp. [mm]')
xlim([0 ttotal])
grid on

subplot(3,1,2)
plot(ts,abs(x_t*1000-x_m*1000),'k')
hold on
legend(['J_2 = ',num2str(J2),' %'])
xlabel('Time [sec]')
ylabel('|error| [mm]')
xlim([0 ttotal])
grid on

subplot(3,1,3)
plot(ts,Ugs/9.81,'b')
legend([earthquake,' ',num2str(escala*100), '%'],'Location','best') 
xlabel('Time [sec]')
ylabel('Earthquake [g]')
xlim([0 ttotal])
grid on

%Zoom figure
desde=round(9/dtsim);  %from
hasta=round(10/dtsim); %to
figure
plot(ts(desde:hasta),x_t(desde:hasta)*1000,'k','LineWidth',1)
hold on
plot(ts(desde:hasta),x_m(desde:hasta)*1000,'r--','LineWidth',2)
legend('x_t','x_m','Orientation','Horizontal','Location','best')
xlabel('Time [sec]')
ylabel('Displacement [mm]')
xlim([desde hasta]*dtsim)
grid on

%% Adaptive parameters at last simulation

%Identified a_i parameters
iddesde=10/dtsim;  %from (sec)
idhasta=35/dtsim;  %to (sec)
A_amb_final=mean([amb0(iddesde:idhasta) amb1(iddesde:idhasta) amb2(iddesde:idhasta) amb3(iddesde:idhasta)]);

%a_i
figure
subplot(2,2,1)
amb0line=plot(ts,amb0,'k');
hold on
plot([ts(iddesde),ts(idhasta)],[A_amb_final(1),A_amb_final(1)],'g:','Linewidth',2);
legend([amb0line],'Parameter a_i','Location','northoutside')
xlabel('Time [sec]')
ylabel('a_0')
xlim([0 ttotal])
grid on

subplot(2,2,2)
plot(ts,amb1,'k')
hold on
amb1id=plot([ts(iddesde),ts(idhasta)],[A_amb_final(2),A_amb_final(2)],'g:','Linewidth',2);
xlabel('Time [sec]')
ylabel('a_1  [sec]')
xlim([0 ttotal])
legend([amb1id],'Identified a_i','Location','northoutside')
grid on

subplot(2,2,3)
plot(ts,amb2,'k')
hold on
plot([ts(iddesde),ts(idhasta)],[A_amb_final(3),A_amb_final(3)],'g:','Linewidth',2);
xlabel('Time [sec]')
ylabel('a_2  [sec^2]')
xlim([0 ttotal])
grid on

subplot(2,2,4)
plot(ts,amb3,'k')
hold on
plot([ts(iddesde),ts(idhasta)],[A_amb_final(4),A_amb_final(4)],'g:','Linewidth',2);
xlabel('Time [sec]')
ylabel('a_3   [sec^3]')
xlim([0 ttotal])
grid on

%Identified bode
Ga = tf(1,[1 a3]);
Gs = tf(a1_b0,[1 b1 b2]);
Gcsi = a2*s;
Gest = tf(1,[me ce ke]);
G0 = feedback(Ga*Gest,Gcsi,-1);
Gp = gain*feedback(Gs*G0,1,-1);

[magGp,phaseGp,wout]=bode(Gp,wbode);   %Actual control plant
[magID,phaseID,wout]=bode(tf(1,flip(A_amb_final)),wbode);  %identified plant
[maginicial,phaseinicial,wout]=bode(tf(1,flip(A_amb_i)),wbode);  %initial plant

figure
subplot(2,1,1)
plot(fbode,db(squeeze(maginicial)),'r')
hold on
plot(fbode,db(squeeze(magGp)),'b--')
plot(fbode,db(squeeze(magID)),'g--')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
grid on
legend('Initial model','Control plant','Identified')
subplot(2,1,2)
plot(fbode,(squeeze(phaseinicial)*pi/180)./wbode'*1000,'r')
hold on
plot(fbode,(squeeze(phaseGp)*pi/180)./wbode'*1000,'b--')
plot(fbode,(squeeze(phaseID)*pi/180)./wbode'*1000,'g--')
xlabel('Frequency [Hz]')
ylabel('Time Delay [msec]')
grid on
