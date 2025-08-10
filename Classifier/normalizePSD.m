function [psd_norm] = normalizePSD(psd, trialInfo, options)
switch options.normStyle
    case 'none'
        psd_norm = psd;
    case 'acrossTrials'
        tags = groupTrials4Norm(trialInfo, options.normTrialsGroup);
        tags_id = unique(tags);
        psd_norm = nan(size(psd));
        
        for tg = 1:length(tags_id)
            nb_trials_tag = sum(tags == tags_id(tg));
            switch options.normModel
                case 'additive'
                    mu = mean(psd(:,:,tags == tags_id(tg)), 3);
                    sigma = std(psd(:,:,tags == tags_id(tg)),[],3);
                    psd_norm(:,:,tags == tags_id(tg)) = (psd(:,:,tags == tags_id(tg)) - repmat(mu,1,1,nb_trials_tag))./repmat(sigma,1,1,nb_trials_tag);
                case 'gain'
                    mu = mean(psd(:,:,tags == tags_id(tg)), 3);
                    psd_norm(:,:,tags == tags_id(tg)) = psd(:,:,tags == tags_id(tg))./repmat(mu,1,1,nb_trials_tag);
                otherwise
                    error('Not coded yet')
            end
        end
    case 'acrossChans'
        switch options.normModel
            case 'additive'
                mu = mean(psd, 1);
                sigma = std(psd,[],1);
                psd_norm = (psd - repmat(mu,nb_chans,1,1))./repmat(sigma,1,1,1);
            case 'gain'
                mu = mean(psd, 1);
                psd_norm = psd./repmat(mu,nb_chans,1,1);
            otherwise
                error('Not coded yet')
        end
    case 'acrossChans&Trials'
        tags = groupTrials4Norm(Labels, options.normTrialsGroup);
        tags_id = unique(tags);
        psd_norm = nan(size(psd));
        
        for tg = 1:length(tags_id)
            nb_trials_tag = sum(tags == tags_id(tg));
            switch options.normModel
                case 'additive'
                    mu = mean(psd(:,:,tags == tags_id(tg)), [1,3]);
                    sigma = std(psd(:,:,tags == tags_id(tg)),[],[1,3]);
                    psd_norm(:,:,tags == tags_id(tg)) = (psd(:,:,tags == tags_id(tg)) - repmat(mu,nb_chans,1,nb_trials_tag))./repmat(sigma,1,1,nb_trials_tag);
                case 'gain'
                    mu = mean(psd(:,:,tags == tags_id(tg)), [1,3]);
                    psd_norm(:,:,tags == tags_id(tg)) = psd(:,:,tags == tags_id(tg))./repmat(mu,nb_chans,1,nb_trials_tag);
                otherwise
                    error('Not coded yet')
            end
        end
    otherwise
        error('Unknown normalization style');
end
end