%% Adaptive gains calibration 

%% Calibration plants
%Actuator parameters
a3=3.3;
b1=425;
b2=1e5;
Ga = tf(1,[1 a3]);
Gs = tf(a1_b0,[1 b1 b2]);
Gcsi = a2*s;
%mass, damping and stiffness ranges for calibration plants
mecals=[10,40];
cecals=[50,1e3];
kecals=[1e4,2e6];

nsimu=50;   %Number of realization for visualization of random control plants
[maginicial,phaseinicial,wout]=bode(tf(1,flip(A_amb_i)),wbode);

figure(21)
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
magGps=zeros(length(fbode),nsimu);
delayGps=zeros(length(fbode),nsimu);
for simu=1:nsimu  

%random virtual mass, damping and stiffness
mecal=mecals(1)+(mecals(2)-mecals(1))*rand(1);
cecal=cecals(1)+(cecals(2)-cecals(1))*rand(1);
kecal=kecals(1)+(kecals(2)-kecals(1))*rand(1);

%bode
Gest = tf(1,[mecal cecal kecal]);
G0 = feedback(Ga*Gest,Gcsi,-1);
Gp = gain*feedback(Gs*G0,1,-1);

[magGp,phaseGp,wout]=bode(Gp,wbode);
magGps(:,simu)=db(squeeze(magGp));
delayGps(:,simu)=(squeeze(phaseGp)*pi/180)./wbode'*1000;

figure(21)
subplot(2,1,1)
plot(fbode,db(squeeze(magGp)),'b--')
subplot(2,1,2)
plot(fbode,(squeeze(phaseGp)*pi/180)./wbode'*1000,'b--')

end
subplot(2,1,1)
legend('Initial model','Possible control plants','Location','Best')

%% Calibration earthquake
load('ElCentroAccelNoScaling.mat');

Ugcal=ElCentroAccel';

dtUg=Ugcal(2,1)-Ugcal(1,1);   %record time step
%Portion of the earthquake utilized for calibration to improve efficiency
Ugdesde=5;         %from sec
Ughasta=12;        %to sec
Ugcal=[(0:dtUg:(Ughasta-Ugdesde))',Ugcal(Ugdesde/dtUg+1:Ughasta/dtUg+1,2)];

tcal=max(Ugcal(:,1));  %total time of calibration simulation

escalas=[30,60]/100;  %Scaling factor range for calibration (uniform distribution)

%% SDOF Numerical substructures for calibration

frefcals=[2.8,4];   %natural frequency range (Hz)
drefcals=[3,5]/100; %damping ratio range

%% Calibration simulation example with two different gains

gain1=diag(10.^[3,2,1,0]);  %Adaptive gains 1
gain2=diag(10.^[6,5,2,1]);  %Adaptive gains 2
 
escala=0.4;        %earthquake scaling factor
drefcal=0.04;      %numerical structure damping ratio
wnrefcal=3.5*2*pi; %numerical structure natural frequency
mecal=20;          %experimental mass
cecal=500;         %experimental damping
kecal=1e6;         %experimental stiffness
    Gest = tf(1,[mecal cecal kecal]);
    G0 = feedback(Ga*Gest,Gcsi,-1);
    Gpcal = feedback(Gs*G0,1,-1);
    [numcal,dencal]=tfdata(Gpcal,'v');

%simulation 1
adaptivegain=gain1;
sim('AMB_L_9_calibration_simulation.slx')
%results
ts1=x_m.Time;
x_t1=x_t.Data(:,1);
x_m1=x_m.Data(:,1);
amb01=amb.Data(:,1);
amb11=amb.Data(:,2);
amb21=amb.Data(:,3);
amb31=amb.Data(:,4);
J21=rms(x_t1-x_m1)/rms(x_t1)*100;

%simulation 2
adaptivegain=gain2;
sim('AMB_L_9_calibration_simulation.slx')
%results
ts2=x_m.Time;
x_t2=x_t.Data(:,1);
x_m2=x_m.Data(:,1);
amb02=amb.Data(:,1);
amb12=amb.Data(:,2);
amb22=amb.Data(:,3);
amb32=amb.Data(:,4);
J22=rms(x_t2-x_m2)/rms(x_t2)*100;

%synchronization
figure
subplot(2,2,1)
plot(ts1,x_t1*1000,'k')
hold on
plot(ts1,x_m1*1000,'r--')
legend('x_t','x_m','Orientation','Horizontal','Location','Best')
xlabel('Time [sec]')
ylabel('Disp. [mm]')
title(['log_{10} (\Gamma_a) = diag([',num2str(diag(log10(gain1))'),'])'])
grid on

subplot(2,2,3)
plot(ts1,abs(x_t1*1000-x_m1*1000),'k')
legend(['J_2 = ',num2str(J21),' %'],'Orientation','Horizontal','Location','Best') 
xlabel('Time [sec]')
ylabel('|error| [mm]')
grid on
ylim([0,1])

subplot(2,2,2)
plot(ts2,x_t2*1000,'k')
hold on
plot(ts2,x_m2*1000,'r--')
title(['log_{10} (\Gamma_b) = diag([',num2str(diag(log10(gain2))'),'])'])
xlabel('Time [sec]')
ylabel('Disp. [mm]')
grid on

subplot(2,2,4)
plot(ts2,abs(x_t2*1000-x_m2*1000),'k')
legend(['J_2 = ',num2str(J22),' %'],'Orientation','Horizontal','Location','Best') 
xlabel('Time [sec]')
ylabel('|error| [mm]')
grid on
ylim([0,1])

% Adaptive parameters
figure
subplot(2,2,1)
plot(ts1,amb01,'b')
hold on
plot(ts2,amb02,'r')
legend('\Gamma_a','\Gamma_b','Orientation','Horizontal','Location','Best')
xlabel('Time [sec]')
ylabel('a_0')
grid on
ylim([0.99,1.01])

subplot(2,2,2)
plot(ts1,amb11,'b')
hold on
plot(ts2,amb12,'r')
xlabel('Time [sec]')
ylabel('a_1  [sec]')
grid on
ylim([0.018,0.026])

subplot(2,2,3)
plot(ts1,amb21,'b')
hold on
plot(ts2,amb22,'r')
xlabel('Time [sec]')
ylabel('a_2  [sec^2]')
grid on
ylim([0.00006,0.00012])

subplot(2,2,4)
plot(ts1,amb31,'b')
hold on
plot(ts2,amb32,'r')
xlabel('Time [sec]')
ylabel('a_3   [sec^3]')
grid on
ylim([0,0.000002])

%% Matlab function for calibration

%input: adaptive gain matrix 
%output: mean J2 error (R2) , std of J2 (dR2) , max J2 (maxR2)

%number of realizations
realizations=100;

%function (x=log10(diag(adaptive gain matrix))
fun=@(x) AMB_L_6_R2function(x,escalas,frefcals,drefcals,mecals,cecals,kecals,realizations,Ga,Gs,Gcsi);

%example with Gamma=diag(10.^[7,4,1,-2])
tic
[prom,desv,maxj2]=fun([7,4,1,-2])
toc



