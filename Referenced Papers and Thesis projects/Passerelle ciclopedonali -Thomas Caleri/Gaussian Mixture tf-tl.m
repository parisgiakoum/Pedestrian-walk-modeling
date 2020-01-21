x1=xlsread('Correlazione Jessica112.xlsx');
tempofdp= x1(:,1);
ForzaM =x1(:,2);

meanF= mean(ForzaM);
meanT= mean(tempofdp);

mu1=[meanT,meanF];
sigma1=(cov(x1));


x2=xlsread('Correlazione Jessica113.xlsx');
tempofdp= x2(:,1);
LunghezzaP =x2(:,2);

meanL= mean(LunghezzaP);
meanT= mean(tempofdp);

mu2=[meanT,meanL];
sigma2=cov(x2);


LL=[mvnrnd(mu1,sigma1,1000);mvnrnd(mu2,sigma2,1000)]; %random;
gm = fitgmdist (LL,2);


LL=[mvncdf(x1);mvncdf(x2)]; %Funzione di distribuzione cumulativa normale multivariata
gm=fitgmdist(LL,2);

LL=[mvnpdf(x1);mvnpdf(x2)]; %Funzione di densita di probabilità multivariata
gm=fitgmdist(LL,2);

A=mvncdf(x1);
B=mvncdf(x2);
XY=[A,B];
scatter(XY(:,1),XY(:,2),10,'.');
hold on
ezcontour(@(x,y)cdf(gm,[x y]),[0.65,0.9],[0.48,0.66])
ezcontour(@(x,y)pdf(gm,[x y]),[0.65,0.9],[0.48,0.66])
ezsurf(@(x,y)pdf(gm,[x y]),[0.65,0.9],[0.48,0.66])