x=xlsread('Correlazione MartinaF22.xlsx'); % Importo i dati di input

TempoF= x(:,1); % Vettore che riporta il tempo misurato fra un passo e  
                  % il successivo

ForzaM =x(:,2);   % Vettore che riporta la forza impressa a ogni passo

LunghezzaP= x(:,3); % Vettore che riporta la distanza fra un passo e  
                    % il successivo
 
meanF= mean(ForzaM);       % Media del vettore delle Forze
meanT= mean(TempoF);       % Media del vettore dei Tempi
meanL= mean(LunghezzaP);   % Media del vettore delle Lunghezze
 
 
 
mu=[meanT;meanF;meanL];        % Vettore delle medie
sigma =var(x);                 % Matrice delle covarianze 
 
[Fi,xi] = ecdf(TempoF); % Restituisce la funzione empirica di 
                        % distribuzione cumulativa Fi valutata nei  
                        % punti xi, utilizzando i dati del vettore 
                        % TempoF.

 
% Ottengo il grafico funzione di distribuzione cumulativa empirica

figure()      
stairs(xi,Fi,'b','LineWidth',1) 
hold on

% Ottengo il grafico della curva lisciata attraverso la funzione ksdensity.

Fi_sm = ksdensity(TempoF ,xi,'function','cdf','bandwidth',0.015);

% Confronto delle due curve ottenute 

plot(xi,Fi_sm,'r-','LineWidth',1) 
xlabel('X1')                                   % Legenda grafico 
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on

--------------------------------------------------------------------------

[Fi,xi] = ecdf(ForzaM); % Restituisce la funzione empirica di 
                        % distribuzione cumulativa Fi valutata nei  
                        % punti xi, utilizzando i dati del vettore 
                        % ForzaM.

 
% Ottengo il grafico funzione di distribuzione cumulativa empirica

figure()      
stairs(xi,Fi,'b','LineWidth',1) 
hold on

% Ottengo il grafico della curva lisciata attraverso la funzione ksdensity.

Fi_sm = ksdensity(ForzaM ,xi,'function','cdf','bandwidth',0.4);

% Confronto delle due curve ottenute 

plot(xi,Fi_sm,'r-','LineWidth',1) 
xlabel('X2')                                   % Legenda grafico 
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on

--------------------------------------------------------------------------
[Fi,xi] = ecdf(LunghezzaP); % Restituisce la funzione empirica di 
                        % distribuzione cumulativa Fi valutata nei  
                        % punti xi, utilizzando i dati del vettore 
                        % LunghezzaP.

 
% Ottengo il grafico funzione di distribuzione cumulativa empirica

figure()      
stairs(xi,Fi,'b','LineWidth',1) 
hold on

% Ottengo il grafico della curva lisciata attraverso la funzione ksdensity.

Fi_sm = ksdensity(LunghezzaP ,xi,'function','cdf','bandwidth',0.03);

% Confronto delle due curve ottenute 

plot(xi,Fi_sm,'r-','LineWidth',1) 
xlabel('X3')                                   % Legenda grafico 
ylabel('Cumulative Probability')
legend('Empirical','Smoothed','Location','NW')
grid on

--------------------------------------------------------------------------