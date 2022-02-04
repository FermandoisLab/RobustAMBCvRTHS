%% Particle Swarm Optimization to obtain Calibrated Adaptive Gains
tic

maxit=12;         %Number of iterations (generations)
population=15;    %Number of evaluations per iteration (swarm size)
lb=[2;0;-2;-4];   %Lower bounds for each gain
ub=[10;8;6;4];    %Upper bound for each gain
inertia=0.5;      %Inertia weight
pacel=0.5;        %weight for best result for each particle
sacel=0.5;        %weight for best result of the swarm

dimension=4;
x0=lb+(ub-lb).*rand(dimension,population); %initial swarm location
v0=-0.5+1*rand(dimension,population);      %initial velocities
func=zeros(1,population);
p_i=x0;
p_i_val=ones(1,population)*inf;
op_location=x0(:,1);
op_value=inf;

%matrices to save the evolution
swarmvalues=zeros(maxit,population);
gain0values=zeros(maxit,population);
gain1values=zeros(maxit,population);
gain2values=zeros(maxit,population);
gain3values=zeros(maxit,population);

%optimization
for i=1:maxit
    for j=1:population
        disp(['Simulation ',num2str((i-1)*population+j),' of ',num2str((maxit+1)*population)])
        [R2,dR2,maxJ2]=fun(x0(:,j));
        func(j)=R2;
        if func(j)<op_value
            op_location=x0(:,j);
            op_value=func(j);
        end
        if p_i_val(j)>func(j)
            p_i(:,j)=x0(:,j);
            p_i_val(j)=func(j);
        end
    end    
swarmvalues(i,:)=func;
gain0values(i,:)=x0(1,:);
gain1values(i,:)=x0(2,:);
gain2values(i,:)=x0(3,:);
gain3values(i,:)=x0(4,:);
v0=inertia*v0+pacel*rand(dimension,population).*(p_i-x0)+sacel*rand(dimension,population).*(op_location-x0);
x0=x0+v0;
for k=1:dimension
    vector=x0(k,:);
    vector(vector<lb(k))=lb(k);
    vector(vector>ub(k))=ub(k);
    x0(k,:)=vector;
end
end

for j=1:population
        disp(['Simulation ',num2str((i)*population+j),' of ',num2str((1+maxit)*population)])
        [R2,dR2,maxJ2]=fun(x0(:,j));
        func(j)=R2;
        if func(j)<op_value
            op_location=x0(:,j);
            op_value=func(j);
        end
    end    
swarmvalues(i+1,:)=func;
gain0values(i+1,:)=x0(1,:);
gain1values(i+1,:)=x0(2,:);
gain2values(i+1,:)=x0(3,:);
gain3values(i+1,:)=x0(4,:);

%Best result
op_location'   %position
op_value       %R2
adaptivegain=diag(10.^op_location');   %adaptivegain defined with the best result
swarmvalues(swarmvalues>100)=100;      %Saturation for plots
toc

%Figures
figure
plot(0:maxit,swarmvalues(:,1))
hold on
for i=2:population
    plot(0:maxit,swarmvalues(:,i))
end
xlabel('Iteration')
ylabel('mean(J_2) [%]')
xticks(0:1:maxit)

figure
subplot(2,2,1)
semilogy(0:maxit,10.^gain0values(:,1))
hold on
for i=2:population
    semilogy(0:maxit,10.^gain0values(:,i))
end
xlabel('Iteration')
ylabel('\Gamma_0')
xticks(1:1:maxit)
yticks(10.^(round(lb(1)-0.001):1:round(ub(1))))
ylim(10.^[round(lb(1)-0.001),round(ub(1))])

subplot(2,2,2)
semilogy(0:maxit,10.^gain1values(:,1))
hold on
for i=2:population
    semilogy(0:maxit,10.^gain1values(:,i))
end
xlabel('Iteration')
ylabel('\Gamma_1')
xticks(1:1:maxit)
yticks(10.^(round(lb(2)-0.001):1:round(ub(2))))
ylim(10.^[round(lb(2)-0.001),round(ub(2))])

subplot(2,2,3)
semilogy(0:maxit,10.^gain2values(:,1))
hold on
for i=2:population
    semilogy(0:maxit,10.^gain2values(:,i))
end
xlabel('Iteration')
ylabel('\Gamma_2')
xticks(1:1:maxit)
yticks(10.^(round(lb(3)-0.001):1:round(ub(3))))
ylim(10.^[round(lb(3)-0.001),round(ub(3))])

subplot(2,2,4)
semilogy(0:maxit,10.^gain3values(:,1))
hold on
for i=2:population
    semilogy(0:maxit,10.^gain3values(:,i))
end
xlabel('Iteration')
ylabel('\Gamma_3')
xticks(0:1:maxit)
yticks(10.^(round(lb(4)-0.001):1:round(ub(4))))
ylim(10.^[round(lb(4)-0.001),round(ub(4))])

x=op_location';

%% R2 Maps around the optimal gain
tic
x=[8.4 , 6.2 , 2.1 , 0.8];    %Optimal gain coeficientes (x=log10(diag(adaptivegain))
puntos=12;   %Grid size puntosxpuntos
nticks=2;
upperbounds=[10,8,6,4];  %Upper and lower bounds for adaptive gains
lowerbounds=[2,0,-2,-4];

% subplot 1/6  Gamma_0 vs Gamma_1
ejex=0;
ejey=1;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_16=zeros(jtotal,itotal);
resultadosdR2_16=zeros(jtotal,itotal);
resultadosmaxJ2_16=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 1/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_16(j,i)=R2;
        resultadosdR2_16(j,i)=dR2;
        resultadosmaxJ2_16(j,i)=maxJ2;
    end
end

% Subplot 2/6  Gamma_0 vs Gamma_2
ejex=0;
ejey=2;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_26=zeros(jtotal,itotal);
resultadosdR2_26=zeros(jtotal,itotal);
resultadosmaxJ2_26=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 2/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_26(j,i)=R2;
        resultadosdR2_26(j,i)=dR2;
        resultadosmaxJ2_26(j,i)=maxJ2;
    end
end

% Subplot 3/6  Gamma_0 vs Gamma_3
ejex=0;
ejey=3;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_36=zeros(jtotal,itotal);
resultadosdR2_36=zeros(jtotal,itotal);
resultadosmaxJ2_36=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 3/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_36(j,i)=R2;
        resultadosdR2_36(j,i)=dR2;
        resultadosmaxJ2_36(j,i)=maxJ2;
    end
end

% Subplot 4/6  Gamma_1 vs Gamma_2
ejex=1;
ejey=2;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_46=zeros(jtotal,itotal);
resultadosdR2_46=zeros(jtotal,itotal);
resultadosmaxJ2_46=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 4/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_46(j,i)=R2;
        resultadosdR2_46(j,i)=dR2;
        resultadosmaxJ2_46(j,i)=maxJ2;
    end
end

% Subplot 5/6  Gamma_1 vs Gamma_3
ejex=1;
ejey=3;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_56=zeros(jtotal,itotal);
resultadosdR2_56=zeros(jtotal,itotal);
resultadosmaxJ2_56=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 5/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_56(j,i)=R2;
        resultadosdR2_56(j,i)=dR2;
        resultadosmaxJ2_56(j,i)=maxJ2;
    end
end

% Subplot 6/6  Gamma_2 vs Gamma_3
ejex=2;
ejey=3;
gainscoef=x;

gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
itotal=length(gains0);
jtotal=length(gains1);
[X,Y]=meshgrid(gains0,gains1);
resultadosR2_66=zeros(jtotal,itotal);
resultadosdR2_66=zeros(jtotal,itotal);
resultadosmaxJ2_66=zeros(jtotal,itotal);
for i=1:itotal
    for j=1:jtotal
        disp(['Evaluation ',num2str((i-1)*itotal+j),'/',num2str(itotal*jtotal), ' | Subplot 6/6'])
        gainscoef(1+ejex)=log10(gains0(i));
        gainscoef(1+ejey)=log10(gains1(j));
        [R2,dR2,maxJ2]=fun(gainscoef);
        resultadosR2_66(j,i)=R2;
        resultadosdR2_66(j,i)=dR2;
        resultadosmaxJ2_66(j,i)=maxJ2;
    end
end
toc

%% Save results 
save('MapsR2dR2maxJ2','resultadosR2_16','resultadosdR2_16','resultadosmaxJ2_16')
save('MapsR2dR2maxJ2','resultadosR2_26','resultadosdR2_26','resultadosmaxJ2_26','-append')
save('MapsR2dR2maxJ2','resultadosR2_36','resultadosdR2_36','resultadosmaxJ2_36','-append')
save('MapsR2dR2maxJ2','resultadosR2_46','resultadosdR2_46','resultadosmaxJ2_46','-append')
save('MapsR2dR2maxJ2','resultadosR2_56','resultadosdR2_56','resultadosmaxJ2_56','-append')
save('MapsR2dR2maxJ2','resultadosR2_66','resultadosdR2_66','resultadosmaxJ2_66','-append')
save('MapsR2dR2maxJ2','x','-append')

%% R2 plot
optimalgain=10.^x;
% Subplot 1/6  Gamma_0 vs Gamma_1
ejex=0;
ejey=1;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

contornos=-1:0.001:2;
figure(30)
contourf(X,Y,log10(resultadosR2_16),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',15)
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
caxis([0 2])
title(['\Gamma_2=10^{',num2str(gainscoef(3),3),'} \Gamma_3=10^{',num2str(gainscoef(4),2),'}'])
colormap(parula(3000))
colorbar('Ticks',log10([1,2,5,10,20,50,100]),'TickLabels',{'R_2<1%','2%','5%','10%','20%','50%','R_2>100%'})

figure(31)
subplot(2,3,1)
contourf(X,Y,log10(resultadosR2_16),contornos,'Linecolor','none')
set(gca,'YScale','log')
set(gca,'XScale','log')
hold on
caxis([0 2])
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_2=10^{',num2str(x(3),3),'} \Gamma_3=10^{',num2str(x(4),2),'}'])
colormap(parula(1000))

% Subplot 2/6  Gamma_0 vs Gamma_2
ejex=0;
ejey=2;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

subplot(2,3,2)
contourf(X,Y,log10(resultadosR2_26),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
caxis([0 2])
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_1=10^{',num2str(x(2),3),'} \Gamma_3=10^{',num2str(x(4),2),'}'])
colormap(parula(1000))

% Subplot 3/6  Gamma_0 vs Gamma_3
ejex=0;
ejey=3;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

subplot(2,3,3)
contourf(X,Y,log10(resultadosR2_36),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
caxis([0 2])
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_1=10^{',num2str(x(2),3),'} \Gamma_2=10^{',num2str(x(3),3),'}'])
colormap(parula(1000))

% Subplot 4/6   Gamma_1 vs Gamma_2
ejex=1;
ejey=2;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

subplot(2,3,4)
contourf(X,Y,log10(resultadosR2_46),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
caxis([0 2])
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_0=10^{',num2str(x(1),3),'} \Gamma_3=10^{',num2str(x(4),2),'}'])
colormap(parula(1000))

% Subplot 5/6    Gamma_1 vs Gamma_3
ejex=1;
ejey=3;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

subplot(2,3,5)
contourf(X,Y,log10(resultadosR2_56),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_0=10^{',num2str(x(1),3),'} \Gamma_2=10^{',num2str(x(3),3),'}'])
colormap(parula(1000))
caxis([0 2])

% Subplot 6/6  Gamma_2 vs Gamma_3
ejex=2;
ejey=3;
gainscoef=x;
gains0=10.^linspace(lowerbounds(ejex+1),upperbounds(ejex+1),puntos);
gains1=10.^linspace(lowerbounds(ejey+1),upperbounds(ejey+1),puntos);
[X,Y]=meshgrid(gains0,gains1);

subplot(2,3,6)
contourf(X,Y,log10(resultadosR2_66),contornos,'Linecolor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
caxis([0 2])
hold on
plot(optimalgain(1+ejex),optimalgain(1+ejey),'or','markersize',4,'markerfacecolor','r')
plot(optimalgain(1+ejex),optimalgain(1+ejey),'xr','markersize',12,'markerfacecolor','r')
xlabel(['\Gamma_',num2str(ejex)])
ylabel(['\Gamma_',num2str(ejey)])
xticks(10.^(lowerbounds(ejex+1):nticks:upperbounds(ejex+1)))
yticks(10.^(lowerbounds(ejey+1):nticks:upperbounds(ejey+1)))
title(['\Gamma_0=10^{',num2str(x(1),3),'} \Gamma_1=10^{',num2str(x(2),3),'}'])
colormap(parula(1000))