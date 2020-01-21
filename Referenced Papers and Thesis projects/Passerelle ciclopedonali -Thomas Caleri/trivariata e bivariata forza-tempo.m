x=xlsread('Correlazione Jessica11.xlsx');
tempofdp= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(tempofdp);
meanL= mean(LunghezzaP);


mu=[meanT,meanF,meanL];
sigma =var(x);
y = mvncdf(x,mu,sigma);
dfittool(y)
z= mvncdf(x);%Funzione di distribuzione cumulativa normale multivariata
dfittool(z)


x1=xlsread('Correlazione Jessica112.xlsx');
h=mvncdf(x1);
dfittool(h);