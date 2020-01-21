x1=xlsread('Correlazione Jessica112.xlsx');
tempofdp= x1(:,1);
ForzaM =x1(:,2);

meanF= mean(ForzaM);
meanT= mean(tempofdp);

mu1=[meanT,meanF];
sigma1=(var(x1));

x3=xlsread('Correlazione Jessica114.xlsx');
ForzaM= x3(:,1);
LunghezzaP =x3(:,2);

meanL= mean(ForzaM);
meanT= mean(LunghezzaP);

mu3=[meanL,meanT];
sigma3=var(x3);

LL=[mvnrnd(mu1,sigma1,1000);mvnrnd(mu3,sigma3,1000)]; %random;
gm = fitgmdist (LL, 2)


LL=[mvncdf(x1);mvncdf(x3)]; %Funzione di distribuzione cumulativa normale multivariata
gm = fitgmdist (LL, 2);

LL=[mvnpdf(x1);mvnpdf(x3)]; %Funzione di densita di probabilità multivariata
gm = fitgmdist (LL, 2);

A=mvncdf(x1);
B=mvncdf(x3);
XY=[A,B];
scatter(XY(:,1),XY(:,2),10,'.')
hold on
ezcontour(@(x,y)cdf(gm,[x y]),[0.65, 0.9],[0.48, 0.66])
ezcontour(@(x,y)pdf(gm,[x y]),[0.65, 0.9],[0.48, 0.66])