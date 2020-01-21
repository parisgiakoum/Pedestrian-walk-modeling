x=xlsread('Correlazione MartinaM22.xlsx');
tempofdp= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(tempofdp);
meanL= mean(LunghezzaP);


mu=[meanT,meanF,meanL];
sigma =var(x);

R= mvnrnd(mu,sigma,415);

T1= R(:,1);
T2= R(:,2);
T3= R(:,3);

h1= kstest2(tempofdp,T1);
h2= kstest2(ForzaM,T2);
h3= kstest2(LunghezzaP,T3);