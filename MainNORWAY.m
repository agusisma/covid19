%% MATLAB codes for Norway

clear;
clc;
tic
%% load data
load Norway.txt; % load data: date | month | susceptible | active cases | cummilative recovered | cummulative death
DATA = Norway;
load CoriNorway_15MAR_2Jun.mat;
intervalC = interval;
load BnRNorway_27FEB_20MAY.mat;
intervalBR = interval;

%% Infectious time
Tinf = 9;
Std_Tinf = 1;

%%
tp  = 30;                                    % prediction time
tf  = length(DATA);                          % simulation time
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
tdp = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf+tp);
dt  = 0.01;
t   = dt:dt:tf;

for i = 2:tf
    DATA(i,7) = (DATA(i,4) + DATA(i,5) + DATA(i,6))-(DATA(i-1,4) + DATA(i-1,5) + DATA(i-1,6));
end
DATA(1:tf,7) = [DATA(1,4); DATA(2:tf,7)];

%% Data matrix
C = [1 0 0 0 0 0;     % We have data of S, I, R, D, and H. The sixth state,
     0 1 0 0 0 0;     % which is Rt, is estimated.
     0 0 1 0 0 0;
     0 0 0 1 0 0
     0 0 0 0 1 0];
%% Parameters
sigma  = 1.96; %95 Confident Interval for infectious time

%% Noise
QF = diag([10 10 10 10 5 0.2]);   % process and measurement covariance matrices
RF = diag([100 10 10 5 1]);       % are considered as tuning parameters

%% For plotting
% Adding a low pass-filter to handle short-term data fluctuations 
windowSize = 300;
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
% Plotting Rt below 1
curve11 = 0*ones(1,tf);
curve22 = 1*ones(1,tf);
x2      = [td, fliplr(td)];
fs      = 24;
fs1     = 36;

%% Simulation
for j = 1:3
% Infection time
Ti     = Tinf-Std_Tinf*sigma+(j-1)*Std_Tinf*sigma;    % infection time with standard dev. 1 day
gamma  = (1-CFR)*(1/Ti);            % recovery rate
kappa  = CFR*1/Ti;                  % death rate

%% Initialization
xhat     = [N-1; 1; 0; 0; 1; 0];   % initial condition
Pplus    = 1000*eye(6);            % since we know excatly the initial conditions
xhatEff  = 0;
% for plotting
xArray       = [];
xhatArray    = [];
xhatEffArray = [];
% extended Kalman filter
for i=1:((tf-1)/dt)
     xhatArray    = [xhatArray xhat]; 
     xhatEffArray = [xhatEffArray xhatEff];      
     % assimilating the reported data
     y = [interp1(0:1:tf-1,DATA(:,3),t);
         interp1(0:1:tf-1,DATA(:,4),t);
         interp1(0:1:tf-1,DATA(:,5),t);
         interp1(0:1:tf-1,DATA(:,6),t);
         interp1(0:1:tf-1,DATA(:,7),t)];
     y = y + sqrt(RF)*[randn randn randn randn randn]'*dt;
     % prediction
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(6)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(6)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5)+(gamma+kappa)*xhat(2)*dt-xhat(5)*dt;
     xhat(6) = xhat(6);
     xhat = xhat + sqrt(QF)*[randn randn randn randn randn randn]'*dt;
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(6)*xhat(2)*dt/N -(gamma+kappa)*xhat(6)*xhat(1)*dt/N 0 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(6)*xhat(2)*dt/N 1+(gamma+kappa)*xhat(6)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0 0;
             0 kappa*dt 0 1 0 0;
             0 (gamma+kappa)*dt 0 0 1-dt 0;
             0 0 0 0 0 1];
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);  % Kalman gain
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(6)-KF*C)*Pmin;
    xhat(6) = max(0,xhat(6));           % the reproduction number cannot be negative
    xhatEff = (xhat(1)/N)*xhat(6);      % calculating the effective repsoduction number
end

%% Plotting

xhatArray(6,:) = filter(b,a,xhatEffArray);

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatIArray  = [];
xhatI       = xhatArray(2,tf);
xhatRArray  = [];
xhatR       = xhatArray(3,tf);
xhatDArray  = [];
xhatD       = xhatArray(4,tf);
xhatHArray  = [];
xhatH       = xhatArray(5,tf);
xhatRtArray = [];
xhatRt      = xhatArray(6,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,(1/dt)*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(2,(1/dt)*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(3,(1/dt)*i);
    xhatDArray  = [xhatDArray xhatD];
    xhatD       = xhatArray(4,(1/dt)*i);
    xhatHArray  = [xhatHArray xhatH];
    xhatH       = xhatArray(5,(1/dt)*i);
    xhatRtArray = [xhatRtArray xhatRt];
    xhatRt      = xhatArray(6,(1/dt)*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatDArray  = [xhatDArray xhatD];
xhatHArray  = [xhatHArray xhatH];
xhatRtArray = [xhatRtArray xhatRt];

M(j,:) = xhatRtArray;

end

%% Relative root mean square error

RRMSEI  = 0;
RRMSER  = 0;
RRMSED  = 0;
RRMSEDR = 0;
for p = 1:length(td)
    RRMSEI  = RRMSEI +((xhatIArray(p)-DATA(p,4)')/max(1,DATA(p,4)))^2;
    RRMSER  = RRMSER +((xhatRArray(p)-DATA(p,5)')/max(1,DATA(p,5)))^2;
    RRMSED  = RRMSED +((xhatDArray(p)-DATA(p,6)')/max(1,DATA(p,6)))^2;
    RRMSEDR = RRMSEDR+((xhatHArray(p)-DATA(p,7)')/max(1,DATA(p,7)))^2;
end

RRMSEI  = (1/length(td))*RRMSEI
RRMSER  = (1/length(td))*RRMSER
RRMSED  = (1/length(td))*RRMSED
RRMSEDR = (1/length(td))*RRMSEDR

RRMSET = (RRMSEI+RRMSER+RRMSED+RRMSEDR)

for l = 1:tf
    curve2(l)      = max(M(:,l));
    xhatRtArray(l) = mean(M(:,l));
    curve1(l)      = min(M(:,l));
end

toc

tdC     = datetime(2020,3,15-1) + caldays(1:length(intervalC(:,2)));
curve3  = intervalC(:,1)';
curve4  = intervalC(:,3)';
x3      = [tdC, fliplr(tdC)];

tdBR    = datetime(2020,2,28-1) + caldays(1:length(intervalBR(:,2)));
curve7  = intervalBR(:,1)';
curve8  = intervalBR(:,3)';
x5      = [tdBR, fliplr(tdBR)];

figure(1)
subplot(2,2,1)
plot(td,DATA(:,4),'-*','linewidth',4,'MarkerSize',20,'MarkerFaceColor','r')
hold on
plot(td,xhatIArray,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
set(gca,'color','none','FontSize',fs)
xlim([min(td) max(tdC)])
title('Active Cases')
legend({'Reported','Estimated'},'Location','northwest','FontSize',fs-4)
grid on
grid minor
ax = gca;
ax.XAxis.TickValues = datetime(2020,DATA(1,2),DATA(1,1)-1) + calmonths(1:tf);
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs-4)
subplot(2,2,2)
plot(td,DATA(:,5),'-*','linewidth',4,'MarkerSize',20,'MarkerFaceColor','r')
hold on
plot(td,xhatRArray,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
set(gca,'color','none','FontSize',fs)
xlim([min(td) max(tdC)])
title('Recovered Cases')
legend({'Reported','Estimated'},'Location','northwest','FontSize',fs-4)
grid on
grid minor
ax = gca;
ax.XAxis.TickValues = datetime(2020,DATA(1,2),DATA(1,1)-1) + calmonths(1:tf);
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs-4)
subplot(2,2,3)
plot(td,DATA(:,6),'-*','linewidth',4,'MarkerSize',20,'MarkerFaceColor','r')
hold on
plot(td,xhatDArray,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
set(gca,'color','none','FontSize',fs)
xlim([min(td) max(tdC)])
title('Death Cases')
legend({'Reported','Estimated'},'Location','northwest','FontSize',fs-4)
grid on
grid minor
ax = gca;
ax.XAxis.TickValues = datetime(2020,DATA(1,2),DATA(1,1)-1) + calmonths(1:tf);
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs-4)
subplot(2,2,4)
plot(td,DATA(:,7),'-*','linewidth',4,'MarkerSize',20,'MarkerFaceColor','r')
hold on
plot(td,xhatHArray,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
set(gca,'color','none','FontSize',fs)
xlim([min(td) max(tdC)])
title('Daily New Cases')
legend({'Reported','Estimated'},'Location','northwest','FontSize',fs-4)
grid on
grid minor
ax = gca;
ax.XAxis.TickValues = datetime(2020,DATA(1,2),DATA(1,1)-1) + calmonths(1:tf);
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs-4)

figure(2)
yyaxis left
plot(td,xhatRtArray,'r','LineWidth',4)
hold on;
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r');
alpha(0.1)
hold on;
plot(tdBR,intervalBR(:,2),'-k','LineWidth',4)
hold on;
inBetween = [curve7, fliplr(curve8)];
fill(x5, inBetween, 'k');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'g','LineWidth',4)
set(gca,'color','none','FontSize',fs1)
ylabel('Effective Reproduction Number')
xlim([min(tdBR)+4, min(tdBR)+18])
ylim([0 8])
grid on
grid minor
ax = gca;
ax.XAxis.TickValues = datetime(2020,2,28-1) + calweeks(1:length(intervalBR(:,2)));
ax.XAxis.TickLabelFormat = 'd/M';
%xtickangle(45)
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs1)
yyaxis right
bar(td,DATA(:,7))
alpha(0.5)
ylabel('Daily New Cases')
legend('Hasan et al.','95% CI','Bettencourt et al.','95% CI','Rt=1','New Cases')

figure(3)
subplot(2,1,1)
yyaxis left
plot(tdC,intervalC(:,2),'-b','LineWidth',4)
hold on;
inBetween = [curve3, fliplr(curve4)];
fill(x3, inBetween, 'b');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'g','LineWidth',4)
set(gca,'color','none','FontSize',fs)
ylabel('R_t')
ylim([0, 2])
grid on
grid minor
xlim([min(tdC), max(tdC)])
ax = gca;
ax.XAxis.TickValues = datetime(2020,3,15-1) + calweeks(1:length(intervalC(:,2)));
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs)
yyaxis right
bar(td,DATA(:,7))
alpha(0.5)
ylabel('Daily New Cases')
legend('Cori et al.','95% CI','Rt=1','New Cases')

subplot(2,1,2)
yyaxis left
plot(td,xhatRtArray,'r','LineWidth',4)
hold on;
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'g','LineWidth',4)
set(gca,'color','none','FontSize',fs)
ylabel('R_t')
ylim([0, 2])
grid on
grid minor
xlim([min(tdC), max(tdC)])
ax = gca;
ax.XAxis.TickValues = datetime(2020,3,15-1) + calweeks(1:length(intervalC(:,2)));
ax.XAxis.TickLabelFormat = 'd/M';
xx = get(gca,'XTickLabel');
set(gca,'XTickLabel',xx,'fontsize',fs)
yyaxis right
bar(td,DATA(:,4))
alpha(0.5)
ylabel('Active Cases')
legend('Hasan et al.','95% CI','Rt=1','Active Cases')