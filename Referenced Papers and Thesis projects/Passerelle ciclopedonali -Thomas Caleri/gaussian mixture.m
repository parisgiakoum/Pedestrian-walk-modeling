x=xlsread('Correlazione MartinaM22.xlsx'); % Matrice dei dati importati del primo campione 

TempoF= x(:,1);        % Vettore che riporta il tempo misurato fra un passo e il successivo
ForzaM =x(:,2);          % Vettore che riporta la Forza impressa ad ogni passo
LunghezzaP= x(:,3);      % Vettore che riporta la lunghezza fra un passo e il successivo 

meanF= mean(ForzaM);       % Media del vettore delle Forze
meanT= mean(TempoF);     % Media del vettore dei Tempi
meanL= mean(LunghezzaP);   % Media del vettore delle Lunghezze


mu=[meanT;meanF;meanL];        % Vettore delle medie
sigma =var(x);                 % Matrice delle covarianze


GMModel = fitgmdist (x, 3) %questa è la struttura che mi dice come è fatta la mistura gaussina che fitta meglio i miei dati

Y = random(GMModel,415); %Ogni riga di Y è una variabile casuale generata dalla distribuzione della miscela gaussiana m- dimensionale gm.

T=Y(:,1);             
F=Y(:,2);
L=Y(:,3);

                          %plot delle tre marginali
dfittool(T)
dfittool(F)
dfittool(L)
                         %Test di Kolmogorov-Smirnov a due campioni
h1=kstest2(TempoF,T);
h2=kstest2(ForzaM,F);
h3=kstest2(LunghezzaP,L);

--------------------------------------------------------------------------
x=xlsread('Correlazione MartinaF22.xlsx'); % Matrice dei dati importati del primo campione 

TempoF= x(:,1);        % Vettore che riporta il tempo misurato fra un passo e il successivo
ForzaM =x(:,2);          % Vettore che riporta la Forza impressa ad ogni passo
LunghezzaP= x(:,3);      % Vettore che riporta la lunghezza fra un passo e il successivo 

meanF= mean(ForzaM);       % Media del vettore delle Forze
meanT= mean(TempoF);     % Media del vettore dei Tempi
meanL= mean(LunghezzaP);   % Media del vettore delle Lunghezze


mu=[meanT;meanF;meanL];        % Vettore delle medie
sigma =var(x);                 % Matrice delle covarianze


GMModel = fitgmdist (x, 3) %questa è la struttura che mi dice come è fatta la mistura gaussina che fitta meglio i miei dati

Y = random(GMModel,143); %Ogni riga di Y è una variabile casuale generata dalla distribuzione della miscela gaussiana m- dimensionale gm.

T=Y(:,1);             
F=Y(:,2);
L=Y(:,3);


                         %Test di Kolmogorov-Smirnov a due campioni
h1=kstest2(TempoF,T);
h2=kstest2(ForzaM,F);
h3=kstest2(LunghezzaP,L);

dfittool(T)            %plot delle tre marginali
dfittool(F)
dfittool(L)
--------------------------------------------------------------------------
x=xlsread('Correlazione MartinaM22.xlsx');
TempoF= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(TempoF);
meanL= mean(LunghezzaP);


mu=[meanT;meanF;meanL];
sigma =var(x);

%GMModel = fitgmdist (x, 3, 'RegularizationValue' , 0.1)
GMModel = fitgmdist (x, 3);

Y = random(GMModel,3) %estituisce n variabili casuali. Ogni riga di Y è una variabile casuale generata dalla distribuzione della miscela gaussiana m- dimensionale gm.

T=0;
F=0;
L=0;
s=65;
for c=1:s
Y = random(GMModel,3);
T1= Y(:,1);
F2= Y(:,2);
L3= Y(:,3);
T=cat(1,T,T1);
F=cat(1,F,F2);
L=cat(1,L,L3);
end

dfittool(T)
dfittool(F)
dfittool(L)

h1=kstest2(TempoF,T)
h2=kstest2(ForzaM,F)
h3=kstest2(LunghezzaP,L)

-------------------------------------------------------------------------


x=xlsread('Correlazione MartinaM11.xlsx');
TempoF= x(:,1);
ForzaM =x(:,2);
LunghezzaP= x(:,3);

meanF= mean(ForzaM);
meanT= mean(TempoF);
meanL= mean(LunghezzaP);


mu=[meanT;meanF;meanL];
sigma =var(x);

%GMModel = fitgmdist (x, 3, 'RegularizationValue' , 0.1)
GMModel = fitgmdist (x, 3);

Y = random(GMModel,3) %estituisce n variabili casuali. Ogni riga di Y è una variabile casuale generata dalla distribuzione della miscela gaussiana m- dimensionale gm.

T=0;
F=0;
L=0;
s=111;
for c=1:s
Y = random(GMModel,3);
T1= Y(:,1);
F2= Y(:,2);
L3= Y(:,3);
T=cat(1,T,T1);
F=cat(1,F,F2);
L=cat(1,L,L3);
end

dfittool(T)
dfittool(F)
dfittool(L)

h1=kstest2(TempoF,T)
h2=kstest2(ForzaM,F)
h3=kstest2(LunghezzaP,L)

