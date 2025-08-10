function [ ] = plotCVstats(CVstats, features2test, params)
%Figures for visualization of classifier performance
n_folds = size(CVstats.training.Err,2);

if params.error
    figure
    errorbar(features2test, 100*mean(CVstats.training.Err,2), 100*std(CVstats.training.Err,[],2)./sqrt(n_folds), 'Linewidth',2)
    hold on
    switch params.CV
        case 'folds'
            errorbar(features2test, 100*mean(CVstats.testing.Err,2), 100*std(CVstats.testing.Err,[],2)./sqrt(n_folds), 'Linewidth',2)
        case 'LOO'
            plot(features2test, 100*mean(CVstats.testing.Err,2), '-o', 'Linewidth',2);
        case 'LOOTR'
            plot(features2test, 100*mean(CVstats.testing.Err,2), '-o', 'Linewidth',2);
    end
    yline(100*CVstats.random.Err, 'k--','LineWidth',1.5);
    yline(100*CVstats.random.Err_signif, 'k:','LineWidth',2);
    ylim([0,100])
    xlabel('Number of features used')
    ylabel('Percentage of error')
    legend({'Training Error', 'Testing Error', 'Chance level', 'p<0.05'},'Location', 'best')
    if isfield(params,'phase')
        % Then it's Tumbler
        switch params.CV
            case 'folds'
                title(['Classification errors (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr-%s_%dfolds_%s', params.name, params.phase, n_folds, params.suffix), {'png'}, []);
            case 'LOO'
                title('Classification errors (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr-%s_LOO_%s', params.name, params.phase, params.suffix), {'png'}, []);
            case 'LOOTR'
                title('Classification errors (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr-%s_LOOTR_%s', params.name, params.phase, params.suffix), {'png'}, []);
        end
    else
        % Then it's VEP
        switch params.CV
            case 'folds'
                title(['Classification errors (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr_%dfolds_%s', params.name, n_folds, params.suffix), {'png'}, []);
            case 'LOO'
                title('Classification errors (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr_LOO_%s', params.name, params.suffix), {'png'}, []);
            case 'LOOTR'
                title('Classification errors (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErr_LOOTR_%s', params.name, params.suffix), {'png'}, []);
        end
    end
end

if params.errorBal
    figure
    errorbar(features2test, 100*mean(CVstats.training.Err_bal,2), 100*std(CVstats.training.Err_bal,[],2)./sqrt(n_folds), 'Linewidth',2)
    hold on
    switch params.CV
        case 'folds'
            errorbar(features2test, 100*mean(CVstats.testing.Err_bal,2), 100*std(CVstats.testing.Err_bal,[],2)./sqrt(n_folds), 'Linewidth',2)
        case 'LOO'
            err_bal_loo = 100*(1-(sum(squeeze(CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(CVstats.testing.CM(:,:,1,:)),[2,3])...
                +sum(squeeze(CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
            plot(features2test, err_bal_loo, '-o', 'Linewidth',2);
        case 'LOOTR'
            err_bal_loo = 100*(1-(sum(squeeze(CVstats.testing.CM(:,:,1,1)),2)./sum(squeeze(CVstats.testing.CM(:,:,1,:)),[2,3])...
                +sum(squeeze(CVstats.testing.CM(:,:,2,2)),2)./sum(squeeze(CVstats.testing.CM(:,:,2,:)),[2,3]))/2);
            plot(features2test, err_bal_loo, '-o', 'Linewidth',2);
    end
    yline(100*CVstats.random.Err_bal, 'k--');
    yline(100*CVstats.random.Err_bal, 'k--','LineWidth',1.5);
    yline(100*CVstats.random.Err_bal_signif, 'k:','LineWidth',2);
    ylim([0,100])
    xlabel('Number of features used')
    ylabel('Percentage of error')
    legend({'Balanced Training Error','Balanced Testing Error', 'Chance level', 'p<0.05'},'Location', 'best')
    
    if isfield(params,'phase')
        % Then it's Tumbler
        switch params.CV
            case 'folds'
                title(['Classification errors (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal-%s_%dfolds_%s', params.name, params.phase, n_folds, params.suffix), {'png'}, []);
            case 'LOO'
                title('Classification errors (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal-%s_LOO_%s', params.name, params.phase, params.suffix), {'png'}, []);
            case 'LOOTR'
                title('Classification errors (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal-%s_LOOTR_%s', params.name, params.phase, params.suffix), {'png'}, []);
        end
    else
        % Then it's VEP
        switch params.CV
            case 'folds'
                title(['Classification errors (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal_%dfolds_%s', params.name, n_folds, params.suffix), {'png'}, []);
            case 'LOO'
                title('Classification errors (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal_LOO_%s', params.name, params.suffix), {'png'}, []);
            case 'LOOTR'
                title('Classification errors (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_ClassErrBal_LOOTR_%s', params.name, params.suffix), {'png'}, []);
        end
    end
end

if params.AUC && ~contains(params.CV,'LOO')
    figure
    errorbar(features2test, 100*mean(CVstats.training.AUC,2),100*std(CVstats.training.AUC,[],2)./sqrt(n_folds), 'Linewidth',2)
    hold on
    errorbar(features2test, 100*mean(CVstats.testing.AUC,2), 100*std(CVstats.testing.AUC,[],2)./sqrt(n_folds), 'Linewidth',2)
    yline(CVstats.random.AUC*100, 'k--');
    xlabel('Number of features used')
    ylabel('AUC (%)')
    ylim([0,100])
    legend({'Training','Testing'},'Location', 'best')
    title('AUC (mean across folds)')
    if isfield(params,'phase')
        % Then it's Tumbler
        saveCurrentFig(params.saveFigFolder,...
            sprintf('%s_AUC-%s_%dfolds_%s', params.name, params.phase, n_folds, params.suffix), {'png'}, []);
    else
        % Then it's VEP
        saveCurrentFig(params.saveFigFolder,...
            sprintf('%s_AUC_%dfolds_%s', params.name, n_folds, params.suffix), {'png'}, []);
    end
end

if params.CM
    n_feats = length(features2test);
    
    if n_feats <= 4
        png_dims = [];
    elseif n_feats<=6
        png_dims = [1000,800];
    elseif n_feats<=9
        png_dims = [1000,1200];
    elseif n_feats<=12
        png_dims = [1200,1200];
    else
        error('Not coded yet')
    end
    
    figure
    for f=1:n_feats
        if n_feats<=6
            subplot(2,ceil(n_feats/2),f)
        elseif n_feats<=12
            subplot(3,ceil(n_feats/3),f)
        else
            error('Not coded yet')
        end
        plotConfMat(squeeze(mean(squeeze(CVstats.training.CM(f,:,:,:)),1))', params.classes)
        title([num2str(features2test(f)),' features used'])
    end
    
    if isfield(params,'phase')
        % Then it's Tumbler
        switch params.CV
            case 'folds'
                suptitle(['Training results (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM-%s_%dfolds_%s', params.name, params.phase, n_folds, params.suffix), {'png'}, png_dims);
            case 'LOO'
                suptitle('Training results (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM-%s_LOO_%s', params.name, params.phase, params.suffix), {'png'}, png_dims);
            case 'LOOTR'
                suptitle('Training results (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM-%s_LOOTR_%s', params.name, params.phase, params.suffix), {'png'}, png_dims);
        end
    else
        % Then it's VEP
        switch params.CV
            case 'folds'
                suptitle(['Training results (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM_%dfolds_%s', params.name, n_folds, params.suffix), {'png'}, png_dims);
            case 'LOO'
                suptitle('Training results (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM_LOO_%s', params.name, params.suffix), {'png'}, png_dims);
            case 'LOOTR'
                suptitle('Training results (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TrainCM_LOOTR_%s', params.name, params.suffix), {'png'}, png_dims);
        end
    end
    
    figure
    for f=1:n_feats
        if n_feats<=6
            subplot(2,ceil(n_feats/2),f)
        elseif n_feats<=12
            subplot(3,ceil(n_feats/3),f)
        else
            error('Not coded yet')
        end
        if contains(params.CV,'LOO')
            plotConfMat(squeeze(mean(squeeze(CVstats.testing.CM(f,:,:,:)),1))', params.classes, true)
        else
            plotConfMat(squeeze(mean(squeeze(CVstats.testing.CM(f,:,:,:)),1))', params.classes)
        end
        title([num2str(features2test(f)),' features used'])
    end
    
    if isfield(params,'phase')
        % Then it's Tumbler
        switch params.CV
            case 'folds'
                suptitle(['Test results (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM-%s_%dfolds_%s', params.name, params.phase, n_folds, params.suffix), {'png'}, png_dims);
            case 'LOO'
                suptitle('Test results (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM-%s_LOO_%s', params.name, params.phase, params.suffix), {'png'}, png_dims);
            case 'LOOTR'
                suptitle('Test results (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM-%s_LOOTR_%s', params.name, params.phase, params.suffix), {'png'}, png_dims);
        end
    else
        % Then it's VEP
        switch params.CV
            case 'folds'
                suptitle(['Test results (mean across ', num2str(n_folds),' folds)'])
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM_%dfolds_%s', params.name, n_folds, params.suffix), {'png'}, png_dims);
            case 'LOO'
                suptitle('Test results (Leave-one-out)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM_LOO_%s', params.name, params.suffix), {'png'}, png_dims);
            case 'LOOTR'
                suptitle('Test results (Leave-one-out trial-based)')
                saveCurrentFig(params.saveFigFolder,...
                    sprintf('%s_TestCM_LOOTR_%s', params.name, params.suffix), {'png'}, png_dims);
        end
    end
end
end

