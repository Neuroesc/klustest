function [thetaR,thetaP,thetaIndx,thetaPowr,thetaLin] = getTHETAfit(input)                                                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes a spike autocorrelogram [times, spikes] and fits a decomposing or damped sine wave to it
% Then we can calculate a theta index and theta power
%
% 29/02/16 Function received from Raju Islam in Shane O'Mara's lab
% 29/02/16 modified to work with existing functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
lags = input(:,1);
count = input(:,2);
m = max(count); 

Ind = find(lags>=0);
lagsX = lags(Ind);
countX = count(Ind);
lagsX = lagsX/1000; % transformed to second for fitting equation;

%% The modified equation from Tsanov and O'Mara (2011) Theta-Modulated Head Direction Cells in the Rat Anterior Thalamus, the original equation is given at the end
coefList= {'a1', 'w1', 'tau1', 'b', 'c1', 'tau2', 'c2','tau3'};
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 6*pi, 0, 0, 0, 0, 0, 0],...
               'Upper',[m, 24*pi, 5, m, m, 0.1, m, 0.01],...
               'StartPoint',[m/2, 12*pi, 0.1, m/2, m/2, 0.05, m/2, 0.005],...
               'MaxIter', 1000000000,...
               'MaxFunEvals', 10000000);
           
fitEquation = '(a1*cos(w1*x))*exp(-abs(x)/tau1)+ b+ c1*exp(-abs(x)/tau2)- c2*exp(-abs(x)/tau3)'; % 'b' outside exponent
ft = fittype(fitEquation,'independent', 'x', 'coefficients', coefList, 'options',fo);
[fittedData,~] = fit(lagsX,countX,ft);
coef = coeffvalues(fittedData);
a1 = coef(1);
w1 = coef(2); 
tau1 = coef(3);
b = coef(4);
c1 = coef(5);
tau2 = coef(6);
c2 = coef(7);
tau3 = coef(8);
envelop = (a1*cos(w1*lagsX)).*exp(-abs(lagsX)/tau1)+ b+ c1*exp(-abs(lagsX)/tau2)- c2*exp(-abs(lagsX)/tau3); % 'b' outside exponent
thetaLin = [flipud(envelop(2:end)); envelop];

thetaIndx = a1/b; % calculate theta index

%% Calculate theta correlation and p value
[thetaR,thetaP] = corr(count,thetaLin,'type','Pearson','rows','pairwise'); % calculate theta correlation R and P

%% Theta power caluclation taken from Winter and Taube (2015)Passive Transport Disrupts Grid Signals in the Parahippocampal Cortex
% Theta power was measured by calculating a theta ratio, which was the difference between the first trough and the second peak of the autocorrelogram. 
half_lag = lags(lags >= 0);
half_theta = thetaLin(lags >= 0);
[Maxima,MaxIdx] = findpeaks(half_theta);
DataInv = 1.01*max(half_theta) - half_theta;
[~,MinIdx] = findpeaks(DataInv);
Minima = half_theta(MinIdx);
if length(MaxIdx) > 1 && ~isempty(MaxIdx) && ~isempty(MinIdx)
    thetaPowr = half_theta(MaxIdx(2)) - half_theta(MinIdx(1));
else
    thetaPowr = NaN;
end % if length(MaxIdx) > 1 && ~isempty(MaxIdx) && ~isempty(MinIdx)


%% The unmodified equation from Tsanov and O'Mara (2011) Theta-Modulated Head Direction Cells in the Rat Anterior Thalamus
% coefList= {'a1', 'w1', 'tau1', 'b', 'c1', 'tau2'};
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0, 12*pi, 0, 0, -2*m, 0],...
%                'Upper',[2*m, 24*pi, 5, 2*m, 2*m, 0.05],...
%                'StartPoint',[m, 16*pi, 0.1, m, m, 0.005],...
%                'MaxIter', 1000000000,...
%                'MaxFunEvals', 10000000);
% fitEquation= '(a1*(sin(w1*x)+1)+ b)*exp(-abs(x)/tau1)+ c1*exp(-x.^2/tau2^2)'; % 'b' inside exponent
% ft = fittype(fitEquation,'independent', 'x', 'coefficients', coefList, 'options',fo);
% [fittedData, gof] = fit(lagsX,countX,ft);
% coef= coeffvalues(fittedData);
% a1= coef(1);
% w1= coef(2); 
% tau1= coef(3);
% b= coef(4);
% c1= coef(5);
% tau2= coef(6);
% envelop= (a1*(sin(w1*lagsX)+1)+ b).*exp(-abs(lagsX)/tau1)+ c1*exp(-lagsX.^2/tau2^2);
% 
% 
% thetaIndex= a1/b;
% f1= w1/ (2*pi);




 + a2*sin(b2*x+c2).*exp(-((x-b1)/c1)^2)




