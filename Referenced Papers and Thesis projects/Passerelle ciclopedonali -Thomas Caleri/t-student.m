x=xlsread('Correlazione Jessica11.xlsx');
% Matrice dei dati importati del primo campione 
TempoF= x(:,1);                           % Vettore che riporta il tempo misurato fra un passo e il successivo
ForzaM =x(:,2);                           % Vettore che riporta la Forza impressa ad ogni passo
LunghezzaP= x(:,3);                       % Vettore che riporta la lunghezza fra un passo e il successivo 

meanF= mean(ForzaM);                      % Media del vettore delle Forze
meanT= mean(TempoF);                      % Media del vettore dei Tempi
meanL= mean(LunghezzaP);                  % Media del vettore delle Lunghezze


mu=[meanT,meanF,meanL];                   % Vettore delle medie
Sigma =cov(x);                            % Matrice delle covarianze

R = corrcoef(x);                  
nu=5;             %gradi di libertà

y = mvtcdf(TempoF,R,415)
M= mvtrnd(R,nu,415); %Genera una matrice random secondo la distribuzione 
                        % t di student 
M1= M(:,1); %Marginale di M rappresentante il tempo che intercorre fra un passo e il successivo
M2= M(:,2); %Marginale di M rappresentatnte la forza impressa ad ogni passo
M3= M(:,3); %Marginale di M rappresentante la lunghezza di ogni passo

%Test di Kolmogorov Smirnoff per ogni varibile del campione con il proprio
%corrispondente generato con la t di student multivariata random.
h1= kstest2(TempoF,M1); 
h2= kstest2(ForzaM,M2);
h3= kstest2(LunghezzaP,M3);


y = mvtcdf(x,R,nu); %Funzione di distribuzione cumulativa t multivariata
z= mvtpdf(x,R,nu) %Funzione di densità di probabilità multivariata t

dfittool(z)
h = ttest(z) %test