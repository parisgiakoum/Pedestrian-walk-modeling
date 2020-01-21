%COPULA per produrre dei campioni che simulano la camminata
%produco delle terne di numeri che hanno le marginali piu o meno che
%vediamo noi e che sono correlate fra di loro secondo i valori di
%correlazione che possiamo misurare

%usare le copule come generatore di campioni casulae di forza, lunghezza
%del passo e tempo, anche senza avere la distribuzione congiunta

x=xlsread('Correlazione Jessica11.xlsx'); 
TempoF= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(TempoF);
meanL= mean(LunghezzaP);

mu=[meanT;meanF;meanL];
sigma =var(x);

[Fi,xi] = ecdf(TempoF);

figure()
stairs(xi,Fi,'b','LineWidth',1)
hold on
Fi_sm = ksdensity(TempoF ,xi,'function','cdf','bandwidth',0.02);
plot(xi,Fi_sm,'r-','LineWidth',1)
xlabel('X1')
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on



nu = 3;
tau = corr(x(:,1),x(:,2),'type','kendall'); %restituisce una matrice del coefficiente di correlazione a coppie tra ciascuna coppia di colonne nelle matrici di input Xe Y.
rho = copulaparam('t', tau, nu, 'type','kendall');

%Successivamente, si utilizza copularnd per generare valori casuali dalla copula
%t e trasformarli utilizzando i cdf inversi non parametrici.
n=539;
U = copularnd('t',[1 rho; rho 1],nu,n);

X1 = ksdensity(x(:,1),U(:,1),'function','icdf','width',0.02);
X2 = ksdensity(x(:,2),U(:,2),'function','icdf','width',.15);

scatterhist(X1,X2,'Direction','out')

dfittool(TempoF)
dfittool(X1)

dfittool(ForzaM)
dfittool(X2)


h1=kstest2(X1,TempoF)
h2=kstest2(X2,ForzaM)



-------------------------------------------------------------------------

[Fi,xi] = ecdf(ForzaM);

figure()
stairs(xi,Fi,'b','LineWidth',1);
hold on
Fi_sm = ksdensity(ForzaM ,xi,'function','cdf','width',0.15);
plot(xi,Fi_sm,'r-','LineWidth',1)
xlabel('X1')
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on


nu = 3;
tau = corr(x(:,2),x(:,3),'type','kendall'); %restituisce una matrice del coefficiente di correlazione a coppie tra ciascuna coppia di colonne nelle matrici di input Xe Y.
rho = copulaparam('t', tau, nu, 'type','kendall');

%Successivamente, utilizzare copularndper generare valori casuali dalla copula
%t e trasformarli utilizzando i cdf inversi non parametrici. 
% La ksdensity function consente di effettuare una stima del kernel della 
% distribuzione e di valutare il cdf inverso ai punti della copula in un'unica fase:

n = 539; 
U = copularnd('t',[1 rho; rho 1],nu,n);

X3 = ksdensity(x(:,2),U(:,1),'function','icdf','width',0.15);
X4 = ksdensity(x(:,3),U(:,2),'function','icdf','width',0.05);

scatterhist(X3,X4,'Direction','out')

dfittool(ForzaM)
dfittool(X3)

dfittool(LunghezzaP)
dfittool(X4)

h3=kstest2(X3,ForzaM)
h4=kstest2(X4,LunghezzaP)
-----------------------------------------------------------------------

[Fi,xi] = ecdf(LunghezzaP);

figure()
stairs(xi,Fi,'b','LineWidth',1);
hold on
Fi_sm = ksdensity(LunghezzaP ,xi,'function','cdf','bandwidth',0.03);
plot(xi,Fi_sm,'r-','LineWidth',1)
xlabel('X1')
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on


nu = 3;
tau = corr(x(:,3),x(:,1),'type','kendall'); %restituisce una matrice del coefficiente di correlazione a coppie tra ciascuna coppia di colonne nelle matrici di input Xe Y.
rho = copulaparam('t', tau, nu, 'type','kendall');

%Successivamente, utilizzare copularndper generare valori casuali dalla copula
%t e trasformarli utilizzando i cdf inversi non parametrici. 
% La ksdensity function consente di effettuare una stima del kernel della 
% distribuzione e di valutare il cdf inverso ai punti della copula in un'unica fase:

n = 539; 
U = copularnd('t',[1 rho; rho 1],nu,n);

X5 = ksdensity(x(:,3),U(:,1),'function','icdf','width',0.05);
X6 = ksdensity(x(:,1),U(:,2),'function','icdf','width',0.02);

scatterhist(X5,X6,'Direction','out')

dfittool(LunghezzaP)
dfittool(X5)


dfittool(TempoF)
dfittool(X6)


h5=kstest2(X5,LunghezzaP)
h6=kstest2(X6,TempoF)
