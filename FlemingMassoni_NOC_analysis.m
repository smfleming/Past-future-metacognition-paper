% Script for reproducing analyses of past/future metacogntion experiment
% Fleming, Massoni, Gajdos, Vergnaud in prep
%
%
% Steve Fleming 2012/2016
% stephen.fleming@ucl.ac.uk

clear all
close all
saveplots = 0;

if saveplots
    addpath('~/Dropbox/Utils/graphics/export_fig/')
    figDir = '~/Dropbox/Research/Metacognition/Paris/results/';
end

% Import data
DATA = importdata('data_meta_bf.txt');
subs = unique(DATA.data(:,2))';

% Initial descriptives and exclusion of outliers
for s = subs
    
    currdata = DATA.data(DATA.data(:,2) == s,:);
    acc(s) = mean(currdata(:,9));
    std_Rconf(s) = std(currdata(:,11));
    std_Pconf(s) = nanstd(currdata(:,14));
    
end

% Plot some descriptives
figure;
subplot(1,3,1);
hist(acc)
ylabel('Frequency');
xlabel('Percent correct');

subplot(1,3,2);
hist(std_Rconf)
ylabel('Frequency');
xlabel('SD(R conf)');

subplot(1,3,3);
hist(std_Pconf)
ylabel('Frequency');
xlabel('SD(P conf)');

% Define exclusions (variance of conf ~= 0)
exc = std_Rconf < 0.02 | std_Pconf < 0.02;

% Correlation of backward/forward conf, controlling for performance
ratings = [0.5 0.6 0.7 0.8 0.9 1.0];
j = 1;
allData = [];
allPconf = [];
for s = subs
    if ~exc(s)
        
        currdata = DATA.data(DATA.data(:,2) == s,:);
        
        %% EXTRACT VARIABLES
        accCode = currdata(:,9);
        RT = currdata(:,10);
        stimCode = currdata(:,7) + 1;
        respCode = currdata(:,8);
        
        temp = currdata(:,14);  % temporary variable for Pconf
        Ptrial = ~isnan(temp);  % only trials associated with P judgment
        Rtrial = isnan(temp);
        Pconf = temp(Ptrial);
        Rconf = currdata(Rtrial,11);
        
        mean_acc(j) = mean(currdata(:,9));
        mean_Racc(j) = mean(currdata(Rtrial,9));
        mean_Pacc(j) = mean(currdata(Ptrial,9));
        mean_Rconf(j) = mean(Rconf);
        mean_Pconf(j) = mean(Pconf);
        
        % Transform ratings from probabilities to discrete levels for ROC
        % code
        for r = 1:length(ratings)
            Rrate(Rconf == ratings(r)) = r;
            Prate(Pconf == ratings(r)) = r;
        end
        
        % Get under and overconfidence for R and P conf
        overRcon(j) = mean_Rconf(j) - mean_Racc(j);
        overPcon(j) = mean_Pconf(j) - mean_Pacc(j);
        
        %% HISTOGRAMS OF COUNTS
        for i = 1:length(ratings)
            histConfR(j,i) = sum(Rconf == ratings(i));
            histConfP(j,i) = sum(Pconf == ratings(i));
        end
        
        %% BRIER SCORE AND DECOMPOSITIONS (calls brier_index)
        [BS,UNC,CI,DI,NDI,ANDI] = brier_index(accCode(Ptrial),Pconf,length(ratings));
        PANDI(j) = ANDI;
        PC(j) = CI;
        [BS,UNC,CI,DI,NDI,ANDI] = brier_index(accCode(Rtrial),Rconf,length(ratings));
        RANDI(j) = ANDI;
        RC(j) = CI;
        
        %% AUROC2
        typeIcode = zeros(1,length(stimCode));
        typeIcode(stimCode == 2 & respCode == 2) = 1;
        typeIcode(stimCode == 1 & respCode == 2) = 2;
        typeIcode(stimCode == 2 & respCode == 1) = 3;
        typeIcode(stimCode == 1 & respCode == 1) = 4;
        
        [Raroc(j)] = type2roc(accCode(Rtrial)', Rrate, length(ratings));
        [Paroc(j)] = type2roc(accCode(Ptrial)', Prate, length(ratings));
        
        j=j+1;
    end
end

%% STATS
Nsubs = length(overRcon);

% Bias (omitting outlier)
[STATS.overcon.corr.r STATS.overcon.corr.p] = corrcoef(overRcon([1:5 7:length(overRcon)]), overPcon([1:5 7:length(overRcon)]));
[h p CI s] = ttest(overRcon([1:5 7:length(overRcon)]), overPcon([1:5 7:length(overRcon)]));
STATS.overcon.ttest.p = p;
STATS.overcon.ttest.T = s.tstat;
% Calibration (omitting outlier)
[STATS.calib.corr.r STATS.calib.corr.p] = corrcoef(RC([1:5 7:length(RC)]), PC([1:5 7:length(RC)]));
[h p CI s] = ttest(RC([1:5 7:length(RC)]), PC([1:5 7:length(RC)]));
STATS.calib.ttest.p = p;
STATS.calib.ttest.T = s.tstat;
% Metacognitive accuracy
[STATS.auroc2.corr.r STATS.auroc2.corr.p] = corrcoef(Raroc, Paroc);
[h p CI s] = ttest(Raroc, Paroc);
STATS.auroc2.ttest.p = p;
STATS.auroc2.ttest.T = s.tstat;
[STATS.ANDI.corr.r STATS.ANDI.corr.p] = corrcoef(RANDI, PANDI);
[h p CI s] = ttest(RANDI, PANDI);
STATS.ANDI.ttest.p = p;
STATS.ANDI.ttest.T = s.tstat;

%% FIGURES
%% Distribution of ratings
h = figure;
meanHist = nanmean(histConfR);
semHist = nanstd(histConfR)./sqrt(j-1);
set(gcf,'Position', [1500 500 520 230]);
subplot(1,2,1)
barWithError(meanHist', semHist', 0.8)
set(gca, 'XTickLabel', {'0.5','0.6','0.7','0.8','0.9','1.0'}, 'FontSize', 12)
xlabel('Retrospective confidence');
ylabel('Frequency');
box off

meanHist = nanmean(histConfP);
semHist = nanstd(histConfP)./sqrt(j-1);
subplot(1,2,2)
barWithError(meanHist', semHist', 0.8)
set(gca, 'XTickLabel', {'0.5','0.6','0.7','0.8','0.9','1.0'}, 'FontSize', 12)
xlabel('Prospective confidence');
ylabel('Frequency');
box off
if saveplots
    export_fig([figDir 'ratingDistribution.png'], '-transparent', '-painters', '-r300', h)
end

%% Correlation plots and statistics
% Correlation between overPcon and overRcon
h = figure; set(gcf,'Position', [1500 500 280 400]);
plot(overRcon, overPcon,'ko ','LineWidth',1, 'MarkerSize', 10);
p = polyfit(overRcon, overPcon, 1);
xvals = linspace(-0.2,0.35,200); pfit = polyval(p,xvals);
hold on
plot(xvals, pfit,'k','LineWidth',2);
% Compute line of best fit omitting outlier
p = polyfit(overRcon([1:5 7:length(overRcon)]), overPcon([1:5 7:length(overPcon)]), 1);
pfit = polyval(p,xvals);
plot(xvals, pfit,'k--','LineWidth',2);
set(gca,'FontSize', 14);
line([-0.2 0.35],[0 0],'LineStyle','--','Color','k')
line([0 0],[-0.3 0.6],'LineStyle','--','Color','k')
ylabel('Prospective over-confidence');
xlabel('Retrospective over-confidence');
box off
axis square
if saveplots
    export_fig([figDir 'confCorr.png'], '-transparent', '-painters', '-r300', h)
end

% Correlation between R-AUROC2 and P-AUROC2
h = figure; set(gcf,'Position', [1500 500 280 400]);
plot(Raroc, Paroc,'ko ','LineWidth',1, 'MarkerSize', 10);
p = polyfit(Raroc, Paroc, 1);
xvals = linspace(-1.5,1,200); pfit = polyval(p,xvals);
hold on
plot(xvals, pfit,'k','LineWidth',2);
set(gca,'YLim',[0.3 0.8],'XLim',[0.3 0.8], 'FontSize', 14);
line([0.3 0.8],[0.5 0.5],'LineStyle','--','Color','k')
line([0.5 0.5],[0.3 0.8],'LineStyle','--','Color','k')
ylabel('Prospective AUROC2');
xlabel('Retrospective AUROC2');
box off
axis square
if saveplots
    export_fig([figDir 'metaCorr.png'], '-transparent', '-painters', '-r300', h)
end

% Correlation between R and P calibration
h = figure; set(gcf,'Position', [1500 500 280 400]);
plot(RC, PC,'ko ','LineWidth',1,'MarkerSize', 10);
p = polyfit(RC, PC, 1);
xvals = linspace(-0.05,0.15,200); pfit = polyval(p,xvals);
hold on
plot(xvals, pfit,'k','LineWidth',2);
% Compute line of best fit omitting outlier
p = polyfit(RC([1:5 7:length(RC)]), PC([1:5 7:length(RC)]), 1);
xvals = linspace(-0.05, 0.15,200); pfit = polyval(p,xvals);
hold on
plot(xvals, pfit,'k--','LineWidth',2);
set(gca,'YLim',[-0.05 0.35],'XLim',[-0.05 0.15], 'FontSize', 14);
line([-0.05 0.15],[0 0],'LineStyle','--','Color','k')
line([0 0],[-0.05 0.35],'LineStyle','--','Color','k')
ylabel('Prospective calibration');
xlabel('Retrospective calibration');
box off
axis square
if saveplots
    export_fig([figDir 'calibCorr.png'], '-transparent', '-painters', '-r300', h)
end

% Correlation between R and P ANDI
h = figure; set(gcf,'Position', [1500 500 280 400]);
plot(RANDI, PANDI,'ko ','LineWidth',1,'MarkerSize', 10);
p = polyfit(RANDI, PANDI, 1);
xvals = linspace(-0.1,0.4,200); pfit = polyval(p,xvals);
hold on
plot(xvals, pfit,'k','LineWidth',2);
set(gca,'YLim',[-0.2 0.1],'XLim',[-0.1 0.4], 'FontSize', 14);
line([-0.1 0.4],[0 0],'LineStyle','--','Color','k')
line([0 0],[-0.2 0.1],'LineStyle','--','Color','k')
ylabel('Prospective ANDI');
xlabel('Retrospective ANDI');
box off
axis square
if saveplots
    export_fig([figDir 'ANDI_corr.png'], '-transparent', '-painters', '-r300', h)
end

%% Mean overconfidence, metacognitive accuracy and ANDI
mean_conf = mean([overPcon([1:5 7:length(overRcon)])' overRcon([1:5 7:length(overRcon)])']);
sem_conf = std([overPcon([1:5 7:length(overRcon)])' overRcon([1:5 7:length(overRcon)])'])./sqrt(Nsubs-1);
h = figure; set(gcf,'Position', [1500 500 270 250]);
bar(1, mean_conf(1), 1, 'FaceColor', [1 1 1], 'LineWidth', 2)
hold on
bar(2, mean_conf(2), 1, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar([1 2], mean_conf, sem_conf, 'k ', 'LineWidth', 2, 'LineStyle', 'None')
set(gca, 'XLim', [0 3], 'XTick', [], 'FontSize', 14)
ylabel('Overconfidence')
box off
if saveplots
    export_fig([figDir 'con_mean.png'], '-transparent', '-painters', '-r300', h)
end

mean_meta = mean([Paroc' Raroc']);
sem_meta = std([Paroc' Raroc'])./sqrt(Nsubs);
h = figure; set(gcf,'Position', [1500 500 270 250]);
bar(1, mean_meta(1), 1, 'FaceColor', [1 1 1], 'LineWidth', 2)
hold on
bar(2, mean_meta(2), 1, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar([1 2], mean_meta, sem_meta, 'k ', 'LineWidth', 2, 'LineStyle', 'None')
set(gca, 'XLim', [0 3], 'YLim', [0.3 0.7], 'XTick', [], 'FontSize', 14)
ylabel('AUROC2')
box off
if saveplots
    export_fig([figDir 'meta_mean.png'], '-transparent', '-painters', '-r300', h)
end

mean_ANDI = mean([PANDI' RANDI']);
sem_ANDI = std([PANDI' RANDI'])./sqrt(Nsubs);
h = figure; set(gcf,'Position', [1500 500 270 250]);
bar(1, mean_ANDI(1), 1, 'FaceColor', [1 1 1], 'LineWidth', 2)
hold on
bar(2, mean_ANDI(2), 1, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar([1 2], mean_ANDI, sem_ANDI, 'k ', 'LineWidth', 2, 'LineStyle', 'None')
set(gca, 'XLim', [0 3], 'XTick', [], 'FontSize', 14)
ylabel('ANDI')
box off
if saveplots
    export_fig([figDir 'ANDI_mean.png'], '-transparent', '-painters', '-r300', h)
end

mean_calib = mean([PC([1:5 7:length(RC)])' RC([1:5 7:length(RC)])']);
sem_calib = std([PC([1:5 7:length(RC)])' RC([1:5 7:length(RC)])'])./sqrt(Nsubs-1);
h = figure; set(gcf,'Position', [1500 500 270 250]);
bar(1, mean_calib(1), 1, 'FaceColor', [1 1 1], 'LineWidth', 2)
hold on
bar(2, mean_calib(2), 1, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar([1 2], mean_calib, sem_calib, 'k ', 'LineWidth', 2, 'LineStyle', 'None')
set(gca, 'XLim', [0 3], 'XTick', [], 'FontSize', 14)
ylabel('Calibration')
box off
if saveplots
    export_fig([figDir 'calib_mean.png'], '-transparent', '-painters', '-r300', h)
end

%% Plot regression coefficients for past/future accuracy
% Calculated by Sebastien on 25/04/16
% Order of betas:
% Pconf1 Rconf4:Rconf1 acc_5:acc_0 RT

% Unstandarised betas
Rbeta = [0.0379 -0.0073 0.0015 0.0100 0.1325 -0.0025 -0.0011 -0.0013 -0.0004 -0.0022 0.0395 -0.0624];
Rbeta_SE = [0.0221259 0.0152354 0.0148035 0.0140334 0.0189095 0.0034212 0.0032547 0.0032585 0.0035524 0.0032691 0.0052513 0.0075398];
Pbeta = [0.1790 0.0365 0.0541 0.1003 0.1167 0.0073 0.0002 0.0011 0.0075 0.0084 -0.0009 0.0004];
Pbeta_SE = [0.0269526 0.0161909 0.0149986 0.0157694 0.0187198 0.0040601 0.0040338 0.0044743 0.0040231 0.0039676 0.0047925 0.0025324];

h = figure;
subplot(1,2,1);
set(gcf,'Position', [1500 500 900 300]);
% P/R conf betas
errorbar(1, Rbeta(1), Rbeta_SE(1),'ro-','LineWidth',2,'MarkerSize',8);
hold on
errorbar(2:5, Rbeta(2:5), Rbeta_SE(2:5),'ro-','LineWidth',2,'MarkerSize',8);
% Acc
errorbar(6:10, Rbeta(6:10), Rbeta_SE(6:10),'ko-','LineWidth',2,'MarkerSize',8);
errorbar(11, Rbeta(11), Rbeta_SE(11),'ko-','LineWidth',2,'MarkerSize',8);
% RT
errorbar(12, Rbeta(12), Rbeta_SE(12),'bo-','LineWidth',2,'MarkerSize',8);
set(gca,'YLim',[-0.1 0.25],'XTick',[1:12],'XTickLabel',{'P-1','R-4','R-3','R-2','R-1','o-5','o-4','o-3','o-2','o-1','o-0','RT'},'FontSize',12);
line([0 12],[0 0],'LineStyle','--','Color','k')
ylabel('Regression coefficient ~ R-conf','FontSize',14);
box off

subplot(1,2,2);
% P/R conf betas
errorbar(1, Pbeta(1), Pbeta_SE(1),'ro-','LineWidth',2,'MarkerSize',8);
hold on
errorbar(2:5, Pbeta(2:5), Pbeta_SE(2:5),'ro-','LineWidth',2,'MarkerSize',8);
% Acc
errorbar(6:10, Pbeta(6:10), Pbeta_SE(6:10),'ko-','LineWidth',2,'MarkerSize',8);
errorbar(11, Pbeta(11), Pbeta_SE(11),'ko-','LineWidth',2,'MarkerSize',8);
% RT
errorbar(12, Pbeta(12), Pbeta_SE(12),'bo-','LineWidth',2,'MarkerSize',8);
set(gca,'YLim',[-0.1 0.25],'XTick',[1:12],'XTickLabel',{'P-1','R-4','R-3','R-2','R-1','o-5','o-4','o-3','o-2','o-1','o-0','RT'},'FontSize',12);
line([0 12],[0 0],'LineStyle','--','Color','k')
ylabel('Regression coefficient ~ P-conf','FontSize',14);
box off

if saveplots
    export_fig([figDir 'regression_plots.png'], '-transparent', '-painters', '-r300', h)
end
