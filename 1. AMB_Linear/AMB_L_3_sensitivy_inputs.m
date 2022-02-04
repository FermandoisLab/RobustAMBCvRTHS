%% Sensitivity analysis for a fixed adaptive gain matrix
% input: adaptive gain matrix
% output: J2 results and simulated random paremeters for N realizations

%Number of realizations
realizations=1000;

%Function
fun2=@(x) AMB_L_7_J2andInputs(x,escalas,frefcals,drefcals,mecals,cecals,kecals,realizations,Ga,Gs,Gcsi);

%Adaptive gains coefficients for evaluation
gains=[9.15 6.2 2.1 -0.8];

%Evaluation
tic
[J2sL,escalaparL,drefcalparL,wnrefcalparL,mecalparL,cecalparL,kecalparL,xtmax]=fun2(gains);
toc
mean(J2sL)
%% Plots

yylim=100;

figure
subplot(2,3,1)
scatter(escalaparL,J2sL)
set(gca,'YScale','log')
grid on
xlabel('Earthquake scale')
ylabel('J_2 [%]')
ylim([0.5 yylim])

subplot(2,3,2)
scatter(drefcalparL*100,J2sL)
set(gca,'YScale','log')
hold on
grid on
xlabel('Ref structure damping %')
ylabel('J_2 [%]')
ylim([0.5 yylim])
title(['log_{10} \Gamma = [',num2str(gains),']'])

subplot(2,3,3)
scatter(wnrefcalparL/2/pi,J2sL)
set(gca,'YScale','log')
hold on
grid on
xlabel('Ref structure freq. Hz')
ylabel('J_2 [%]')
ylim([0.5 yylim])

subplot(2,3,4)
scatter(mecalparL,J2sL)
set(gca,'YScale','log')
hold on
grid on
xlabel('Exp. sub. mass kg')
ylabel('J_2 [%]')
ylim([0.5 yylim])

subplot(2,3,5)
scatter(cecalparL,J2sL)
set(gca,'YScale','log')
hold on
grid on
xlabel('Exp. sub. damping Ns/m')
ylabel('J_2 [%]')
ylim([0.5 yylim])

subplot(2,3,6)
scatter(kecalparL,J2sL)
set(gca,'YScale','log')
hold on
grid on
xlabel('Exp. sub. stiffnes N/m')
ylabel('J_2 [%]')
ylim([0.5 yylim])

%% histograms
J2sL(J2sL>=100)=100;
figure
histogram(J2sL,25)
xlabel('J_2 %')
ylabel('Freq')
title(['log_{10} \Gamma = [',num2str(gains),' ] ;  N = ',num2str(realizations),' realizations'])

figure
subplot(2,3,1)
histogram(escalaparL,10)
xlabel('Earthquake scale')
ylabel('Freq')

subplot(2,3,2)
histogram(drefcalparL*100,10)
xlabel('Ref structure damping %')
ylabel('Freq')
title(['log_{10} \Gamma = [',num2str(gains),' ]'])

subplot(2,3,3)
histogram(wnrefcalparL/2/pi,10)
xlabel('Ref structure freq. Hz')
ylabel('Freq')

subplot(2,3,4)
histogram(mecalparL,10)
xlabel('Exp. sub. mass kg')
ylabel('Freq')

subplot(2,3,5)
histogram(cecalparL,10)
xlabel('Exp. sub. damping Ns/m')
ylabel('Freq')

subplot(2,3,6)
histogram(kecalparL,10)
xlabel('Exp. sub. stiffnes N/m')
ylabel('Freq')

