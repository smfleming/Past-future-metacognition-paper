function [params resErrs out] = fitPconf(data, model, pArray)

options = optimset('Display','Iter');
% LB = zeros(1,length(pArray));
% UB = ones(1,length(pArray));
% [params resErrs] = fmincon(@fitfunc,pArray,[],[],[],[],LB,UB,[],options);
[params resErrs] = fminsearch(@fitfunc,pArray,options);

    function resErrs = fitfunc(pArray)
        
        switch model              
                
            case 'obj'
                
                a = pArray(1);
                q(1) = 0.75;
                for t = 1:length(data.acc)
                    % Update rule
                    q(t+1) = q(t) + a.*(data.acc(t) - q(t));
                    % Reset q if new subject
                    if t ~= 1 & data.sub(t) ~= data.sub(t-1)
                        q(t+1) = 0.75;
                    end
                end
                
            case 'subj'
                
                y = pArray(1);
                q(1) = 0.75;
                for t = 1:length(data.acc)
                    % Update rule
                    q(t+1) = q(t) + y.*(data.Rconf(t) - q(t));
                end
                
            case 'intercept'
                q(1) = 0.75;
                for t = 1:length(data.acc)
                    q(t+1) = q(t);
                end
        end
        
        % Do regression on observed Pconf, get deviance to minimise
        [betas resErrs stats] = glmfit(q(data.Ptrial), data.Pconf(data.Ptrial), 'normal', 'link', 'identity');
        out.Ypred = glmval(betas, q(data.Ptrial), 'identity');
        out.stats = stats;
    end
end