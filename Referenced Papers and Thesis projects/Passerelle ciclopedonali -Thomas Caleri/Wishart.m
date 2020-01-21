x=xlsread('Correlazione Jessica11.xlsx');
tempofdp= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(tempofdp);
meanL= mean(LunghezzaP);


mu=[meanT,meanF,meanL];
sigma =cov(x);
df=2;

W = wishrnd(sigma,df);
Tau=inv(sigma);

W1 = iwishrnd(Tau,df);

T=0;
F=0;
L=0;
s=139;
for c=1:s
W = wishrnd(sigma,df);
T1= W(:,1);
F2= W(:,2);
L3= W(:,3);
T=cat(1,T,T1);
F=cat(1,F,F2);
L=cat(1,L,L3);
end

h1= kstest2(tempofdp,T)
h2= kstest2(ForzaM,F)
h3= kstest2(LunghezzaP,L)

% non ha tanto significato perchè è una distribuzione di probabilità sulle
% matrici semidefinite posotive
%serve in statistica per gestire la statistica delle matrici di
%correlazioni delle distribuzioni normale