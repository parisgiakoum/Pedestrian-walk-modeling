%l teorema di Sklar afferma che ogni copula è una funzione 
%di distribuzione congiunta, avente come argomenti le distribuzioni marginali. 
%Inoltre vale anche il contrario: ogni distribuzione congiunta ha una copula,
%e, se le marginali sono continue, essa è unica.




x=xlsread('Correlazione MartinaM22.xlsx');
tempofdp= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(tempofdp);
meanL= mean(LunghezzaP);


mu=[meanT,meanF,meanL];
sigma =var(x);

R = corrcoef(x);

u = copularnd('Gaussian',R,415);

T1= u(:,1);
T2= u(:,2);
T3= u(:,3);

h1= kstest2(tempofdp,T1);
h2= kstest2(ForzaM,T2);
h3= kstest2(LunghezzaP,T3);