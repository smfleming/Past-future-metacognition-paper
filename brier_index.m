function [BS,UNC,CI,DI,NDI,ANDI] = brier_index(accuracy,confidence,n_conf)
%S.M. 10/03/13
%Outputs: 
%BS = Brier Score (Brier, 1950)
%Murphy's Decomposition (Murphy, 1973) 
    %BS = UNC (uncertainty) + CI (calibration) - DI(discrimination/resolution)
%Adjusted Normalized DI (Yaniv et al., 1991)
    %NDI = DI adjusted from variance, 
    %ANDI = NDI adjusted from sample size and scale size
%Inputs: 
    %accuracy = vector of accuracy
    %confidence = vector of confidence
    %n_conf = number of confidence levels available on the scale used

if isrow(accuracy)
    accuracy = accuracy';
end
if isrow(confidence)
    confidence = confidence';
end
n = length(accuracy);
values = unique(confidence);
n_used = length(values); 
N = [];
N_1 = [];
N_0 = [];
for i = 1:n_used
    a = sum(accuracy(confidence==values(i)));
    b = numel(accuracy(confidence==values(i)));
    c = b - a;
    N = [N b];
    N_1 = [N_1 a];
    N_0 = [N_0 c];
    clear a b c
end

sr = mean(accuracy);
UNC = sr*(1-sr);

CI = [];
for i = 1:n_used
    a = N(i)*((values(i)-(N_1(i)/N(i))).^2);
    CI = [CI a];
    clear a 
end
CI = sum(CI) / n;

DI = [];
for i = 1:n_used
    a = N(i)*(((N_1(i)/N(i))-sr).^2);
    DI = [DI a];
    clear a 
end
DI = sum(DI) / n;

BS = UNC + CI - DI;

NDI = DI/UNC;

ANDI = (n*NDI - n_conf + 1)/(n - n_conf + 1);