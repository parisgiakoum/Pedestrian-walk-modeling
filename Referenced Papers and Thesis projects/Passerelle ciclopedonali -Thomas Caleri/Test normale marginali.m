x=xlsread('Correlazione Jessica11.xlsx');
tempofdp= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(tempofdp);
meanL= mean(LunghezzaP);


mu=[meanT,meanF,meanL];
sigma =var(x);

a= kstest(tempofdp);
b= kstest(ForzaM);
c= kstest(LunghezzaP);