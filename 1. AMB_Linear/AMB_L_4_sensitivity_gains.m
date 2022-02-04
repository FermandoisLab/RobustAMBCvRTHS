%% Sensitivity analysis for different adaptive gains
%
lb=[2;0;-2;-4]';  %lower bounds for adaptive gains (Gamma_1; Gamma_2; Gamma_3; Gamma_4) 
ub=[10;8;6;4]';   %upper bounds for adaptive gains
puntos = 10;     %Number of realizations
Gammas = lb+(ub-lb).*lhsdesign(puntos,4)  %Realizations

% Gammas(:,2)=6.2  %only for fixing Gamma_2
%% Evaluation
R2sen=zeros(puntos,1);
for i=1:puntos
R2sen(i)=fun(Gammas(i,:));
i
end

% save('SensivityGfix','Gammas','J2sen')

%% 
yylim=1000;

figure
subplot(2,2,1)
scatter(10.^Gammas(:,1),R2sen)
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
grid on
xlabel('\Gamma_0')
ylabel('R_2 [%]')
xticks(10.^(round(lb(1)-0.001):1:round(ub(1))))
xlim(10.^[round(lb(1)-0.001),round(ub(1))])
ylim([0.5 yylim])

subplot(2,2,2)
scatter(10.^Gammas(:,2),R2sen)
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
grid on
xlabel('\Gamma_1')
ylabel('R_2 [%]')
xticks(10.^(round(lb(2)-0.001):1:round(ub(2))))
xlim(10.^[round(lb(2)-0.001),round(ub(2))])
ylim([0.5 yylim])

subplot(2,2,3)
scatter(10.^Gammas(:,3),R2sen)
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
grid on
xlabel('\Gamma_2')
ylabel('R_2 [%]')
xticks(10.^(round(lb(3)-0.001):1:round(ub(3))))
xlim(10.^[round(lb(3)-0.001),round(ub(3))])
ylim([0.5 yylim])

subplot(2,2,4)
scatter(10.^Gammas(:,4),R2sen)
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
grid on
xlabel('\Gamma_3')
ylabel('R_2 [%]')
xticks(10.^(round(lb(4)-0.001):2:round(ub(4))))
xlim(10.^[round(lb(4)-0.001),round(ub(4))])
ylim([0.5 yylim])


