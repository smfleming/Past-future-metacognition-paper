% Prospective confidence learning model comparison
% Fleming, Massoni, Gajdos, Vergnaud in prep
%
%
% Steve Fleming 2012/2016
% stephen.fleming@ucl.ac.uk

clear all
close all
saveplots = 0;

if saveplots
    figDir = '~/Dropbox/Research/Metacognition/Paris/results/';
end

DATA = importdata('data_meta_bf.txt');

subs = unique(DATA.data(:,2))';
plots = 1;
exampleSubs = [1 15 37];

% Initial descriptives and exclusion of outliers
for s = subs
    
    currdata = DATA.data(DATA.data(:,2) == s,:);
    acc(s) = mean(currdata(:,9));
    std_Rconf(s) = std(currdata(:,11));
    std_Pconf(s) = nanstd(currdata(:,14));
    
end

%% Load data, fit individual subjects
% Define exclusions (variance of conf = 0)
exc = std_Rconf < 0.02 | std_Pconf < 0.02;
allData = [];
count = 1;
cols = {'r','g','y'};
for s = subs
    if ~exc(s)
        currdata = DATA.data(DATA.data(:,2) == s,:);
        data.acc = currdata(:,9);
        data.Rconf = currdata(:,11);
        data.Pconf = currdata(:,14);
        data.Ptrial = ~isnan(data.Pconf);
        data.sub = currdata(:,2);
        
        model = {'intercept','obj','subj'};
        for m = 1:length(model)
            clear pArray
            switch model{m}
                case 'obj'
                    pArray(1) = 0.1;   % starting alpha
                case 'subj'
                    pArray(1) = 0.1;
                case 'intercept'
                    pArray = [];
            end
            if ~strcmp(model{m},'intercept')
                [params{m}(count,:) dev(count,m) out] = fitPconf(data, model{m}, pArray);
            else
                [temp dev(count,m) out] = fitPconf(data, model{m}, pArray);
            end
            LL(count,m) = sum(log(normpdf(data.Pconf(data.Ptrial),out.Ypred,out.stats.sfit)));
            BIC(count,m) = -2.*LL(count,m) + length(pArray).*(log(40));    % log(40) is number of trials
            
            if plots
                figure(1);
                subplot(7,7,count);
                plot(data.Pconf(data.Ptrial),'LineWidth',2);
                hold on
                plot(out.Ypred,'LineWidth',2);
                set(gca,'YLim',[0.5 1]);
                if m == 3
                    text(20, 0.6, ['alpha = ' num2str(params{m}(count))]);
                end
            end
            if count == exampleSubs(1)
                h1 = figure(length(model)+1);
                subplot(3,length(model),m);
                plot(data.Pconf(data.Ptrial),'LineWidth',2);
                hold on
                plot(out.Ypred,'r','LineWidth',2);
                set(gca,'YLim',[0.5 1],'FontSize',12);
                xlabel('Trials','FontSize',12);
                if m == 1
                    ylabel('P-conf','FontSize',14);
                    title('Intercept-only')
                elseif m == 2
                    title('Model A (outcomes)')
                elseif m == 3
                    title('Model B (R-conf)')
                end
            end
            if count == exampleSubs(2)
                figure(length(model)+1);
                subplot(3,length(model),m+length(model));
                plot(data.Pconf(data.Ptrial),'LineWidth',2);
                hold on
                plot(out.Ypred,'r','LineWidth',2);
                set(gca,'YLim',[0.5 1],'FontSize',12);
                xlabel('Trials','FontSize',12);
                if m == 1
                    ylabel('P-conf','FontSize',14);
                end
            end
            if count == exampleSubs(3)
                figure(length(model)+1);
                subplot(3,length(model),m+(length(model)*2));
                plot(data.Pconf(data.Ptrial),'LineWidth',2);
                hold on
                plot(out.Ypred,'r','LineWidth',2);
                set(gca,'YLim',[0.5 1],'FontSize',12);
                xlabel('Trials','FontSize',12);
                if m == 1
                    ylabel('P-conf','FontSize',14);
                elseif m == 3
                    legend({'Data','Model'},'Location','SouthEast')
                    legend boxoff
                end
            end
        end
        allData = [allData; currdata];
        count = count+1;
    end
end

box off
if saveplots
    print(h1,'-dpng','-r300',[figDir 'ExampleFits.png']);
end

diff_BIC = BIC(:,2) - BIC(:,3);
h2=figure; bar(sort(diff_BIC));
line([0 40],[-6 -6],'LineStyle','--','Color','r','LineWidth',1.5)
line([0 40],[6 6],'LineStyle','--','Color','r','LineWidth',1.5)
box off
ylabel('Difference in BIC')
xlabel('Subjects')
set(gca, 'FontSize', 20)
if saveplots
    print(h2,'-dpng','-r300',[figDir 'DiffBIC.png']);
end
