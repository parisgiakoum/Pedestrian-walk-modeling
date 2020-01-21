%Genero una matrice random che segue la normale con media e sigma dati e
%confronto le marginali


x=xlsread('Correlazione MartinaF22.xlsx'); % Matrice dei dati importati del primo campione 
tempoF= x(:,1);                           % Vettore che riporta il tempo misurato fra un passo e il successivo
ForzaM =x(:,2);                           % Vettore che riporta la Forza impressa ad ogni passo
LunghezzaP= x(:,3);                       % Vettore che riporta la lunghezza fra un passo e il successivo 

meanF= mean(ForzaM);                      % Media del vettore delle Forze
meanT= mean(tempoF);                      % Media del vettore dei Tempi
meanL= mean(LunghezzaP);                  % Media del vettore delle Lunghezze


mu=[meanT,meanF,meanL];                   % Vettore delle medie
Sigma =var(x);                            % Matrice delle covarianze

R= mvnrnd(mu,Sigma,415);          % restituisce una matrice di vettori casuali
                                  % scelti dalla distribuzione normale multivariata      
                                  % con media Mu e covarianza Sigma

T= R(:,1);  % Marginale di R, vettore dei tempi generati secondo la normale 
F= R(:,2);  % Marginale di R, vettore delle forze generati secondo la normale 
L= R(:,3);  % Marginale di R, vettore delle lunghezze generati secondo la normale

h1= kstest2(tempoF,T);           %Test di Kolmogorov-Smirnov a due campioni
h2= kstest2(ForzaM,F);
h3= kstest2(LunghezzaP,L);

dfittool(T)         %Fitting di T
dfittool(F)         %Fitting di F
dfittool(L)         %Fitting di L
