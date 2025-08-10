function fileName = makePSDFileName(type, sID, options)
% Create standardized file names for labels and features data (works both
% for psd features and amp features)
% Used for Tumbler and VEP tasks

switch options.event
    case 'EyesOpening'
        str_ev = 'eyeOp';
            case 'TrialStart'
        str_ev = 'TrStart';
    otherwise
        error('Unknown event')
end

switch type
    case 'labels'
        if isfield(options,'phase')
            % Then it's tumbler task
            if options.chunks == 0
                fileName = sprintf('%s_%s_labels_%s_noChk.mat', sID, str_ev, options.phase);
            else
                fileName = sprintf('%s_%s_labels_%s_%ssChk.mat', sID, str_ev, options.phase,...
                    replace(sprintf('%g',options.chunks),'.','-'));
            end
        elseif isfield(options,'epochWdW')
            % Then it's VEP
            if options.chunks == 0
                fileName = sprintf('%s_%s_%s_labels_noChk.mat', sID, str_ev, options.epochWdW);
            else
                fileName = sprintf('%s_%s_%s_labels_%ssChk.mat', sID, str_ev, options.epochWdW,...
                    replace(sprintf('%g',options.chunks),'.','-'));
            end
        else
            error('Unknown case')
        end
    case 'amp'
        if isfield(options,'phase')
            % Then it's tumbler task
            if strcmp(options.normStyle, 'none')
                if options.chunks == 0
                    fileName = sprintf('%s_%s_amp_%s_noNorm_noChunk.mat', sID, str_ev, options.phase);
                else
                    fileName = sprintf('%s_%s_amp_%s_noNorm_%sChunk.mat', sID, str_ev, options.phase,...
                        replace(num2str(options.chunks),'.','-'));
                end
            else
                error('Not checked')
                str_trGroup = 'per';
                for s = 1:numel(options.normTrialsGroup)
                    switch options.normTrialsGroup{s}
                        case 'perSubject'
                            str_trGroup = sprintf('%sSbj-', str_trGroup);
                        case 'perCondition'
                            str_trGroup = sprintf('%sCnd-', str_trGroup);
                        case 'perBlock'
                            str_trGroup = sprintf('%sBlk-', str_trGroup);
                        case 'perTrialType'
                            str_trGroup = sprintf('%sTrlTyp-', str_trGroup);
                        case 'perAnswer'
                            str_trGroup = sprintf('%sAns-', str_trGroup);
                        otherwise
                            error('Unknown trials grouping')
                    end
                end
                % remove the last '-'
                str_trGroup = str_trGroup(1:end-1);
                
                switch options.normStyle
                    case 'acrossTrials'
                        str_style = 'acrTrls';
                    case 'acrossChans'
                        str_style = 'acrChns';
                    case 'acrossChans&Trials'
                        str_style = 'acrChns-Trls';
                end
                
                switch options.normModel
                    case 'additive'
                        str_model = 'normAdd';
                    case 'gain'
                        str_model = 'normGain';
                end
                
                fileName = sprintf('%s_%s_amp_%s_%s_%s_%s.mat', sID, str_ev, options.phase,...
                    str_model, str_trGroup, str_style);
            end
        elseif isfield(options,'epochWdW')
            % Then it's VEP
            if strcmp(options.normStyle, 'none')
                if options.chunks == 0
                    fileName = sprintf('%s_%s_%s_amp_noNorm_noChunk.mat', sID, str_ev, options.epochWdW);
                else
                    fileName = sprintf('%s_%s_%s_amp_noNorm_%sChunk.mat', sID, str_ev, options.epochWdW,...
                        replace(num2str(options.chunks),'.','-'));
                end
            else
                error('Not checked')
                str_trGroup = 'per';
                for s = 1:numel(options.normTrialsGroup)
                    switch options.normTrialsGroup{s}
                        case 'perSubject'
                            str_trGroup = sprintf('%sSbj-', str_trGroup);
                        case 'perCondition'
                            str_trGroup = sprintf('%sCnd-', str_trGroup);
                        case 'perBlock'
                            str_trGroup = sprintf('%sBlk-', str_trGroup);
                        case 'perTrialType'
                            str_trGroup = sprintf('%sTrlTyp-', str_trGroup);
                        case 'perAnswer'
                            str_trGroup = sprintf('%sAns-', str_trGroup);
                        otherwise
                            error('Unknown trials grouping')
                    end
                end
                % remove the last '-'
                str_trGroup = str_trGroup(1:end-1);
                
                switch options.normStyle
                    case 'acrossTrials'
                        str_style = 'acrTrls';
                    case 'acrossChans'
                        str_style = 'acrChns';
                    case 'acrossChans&Trials'
                        str_style = 'acrChns-Trls';
                end
                
                switch options.normModel
                    case 'additive'
                        str_model = 'normAdd';
                    case 'gain'
                        str_model = 'normGain';
                end
                
                fileName = sprintf('%s_%s_amp_%s_%s_%s_%s.mat', sID, str_ev, options.phase,...
                    str_model, str_trGroup, str_style);
            end
        else
            error('Unknown case')
        end
        
    case 'psd'
        if isfield(options,'phase')
            % Then it's tumbler task
            if strcmp(options.normStyle, 'none')
                if options.chunks == 0
                    fileName = sprintf('%s_%s_psd_%s_noNorm_noChunk.mat', sID, str_ev, options.phase);
                else
                    fileName = sprintf('%s_%s_psd_%s_noNorm_%sChunk.mat', sID, str_ev, options.phase,...
                        replace(num2str(options.chunks),'.','-'));
                end
            else
                str_trGroup = 'per';
                for s = 1:numel(options.normTrialsGroup)
                    switch options.normTrialsGroup{s}
                        case 'perSubject'
                            str_trGroup = sprintf('%sSbj-', str_trGroup);
                        case 'perCondition'
                            str_trGroup = sprintf('%sCnd-', str_trGroup);
                        case 'perBlock'
                            str_trGroup = sprintf('%sBlk-', str_trGroup);
                        case 'perTrialType'
                            str_trGroup = sprintf('%sTrlTyp-', str_trGroup);
                        case 'perAnswer'
                            str_trGroup = sprintf('%sAns-', str_trGroup);
                        otherwise
                            error('Unknown trials grouping')
                    end
                end
                % remove the last '-'
                str_trGroup = str_trGroup(1:end-1);
                
                switch options.normStyle
                    case 'acrossTrials'
                        str_style = 'acrTrls';
                    case 'acrossChans'
                        str_style = 'acrChns';
                    case 'acrossChans&Trials'
                        str_style = 'acrChns-Trls';
                end
                
                switch options.normModel
                    case 'additive'
                        str_model = 'normAdd';
                    case 'gain'
                        str_model = 'normGain';
                end
                
                fileName = sprintf('%s_%s_psd_%s_%s_%s_%s.mat', sID, str_ev, options.phase,...
                    str_model, str_trGroup, str_style);
                
                %             if options.chunks == 0
                %                 fileName = sprintf('%s_%s_psd_%s_%s_%s_%s_noChunk.mat', sID, str_ev, options.phase,...
                %                     str_model, str_trGroup, str_style);
                %             else
                %                 fileName = sprintf('%s_%s_psd_%s_%s_%s_%s_%sChunk.mat', sID, str_ev, options.phase,...
                %                     str_model, str_trGroup, str_style, replace(num2str(options.chunks),'.','-'));
                %             end
            end
            
        elseif isfield(options,'epochWdW')
            % Then it's VEP
            if strcmp(options.normStyle, 'none')
                if options.chunks == 0
                    fileName = sprintf('%s_%s_%s_psd_noNorm_noChunk.mat', sID, str_ev, options.epochWdW);
                else
                    fileName = sprintf('%s_%s_%s_psd_noNorm_%sChunk.mat', sID, str_ev, options.epochWdW,...
                        replace(num2str(options.chunks),'.','-'));
                end
            else
                str_trGroup = 'per';
                for s = 1:numel(options.normTrialsGroup)
                    switch options.normTrialsGroup{s}
                        case 'perSubject'
                            str_trGroup = sprintf('%sSbj-', str_trGroup);
                        case 'perCondition'
                            str_trGroup = sprintf('%sCnd-', str_trGroup);
                        case 'perBlock'
                            str_trGroup = sprintf('%sBlk-', str_trGroup);
                        case 'perTrialType'
                            str_trGroup = sprintf('%sTrlTyp-', str_trGroup);
                        case 'perAnswer'
                            str_trGroup = sprintf('%sAns-', str_trGroup);
                        otherwise
                            error('Unknown trials grouping')
                    end
                end
                % remove the last '-'
                str_trGroup = str_trGroup(1:end-1);
                
                switch options.normStyle
                    case 'acrossTrials'
                        str_style = 'acrTrls';
                    case 'acrossChans'
                        str_style = 'acrChns';
                    case 'acrossChans&Trials'
                        str_style = 'acrChns-Trls';
                end
                
                switch options.normModel
                    case 'additive'
                        str_model = 'normAdd';
                    case 'gain'
                        str_model = 'normGain';
                end
                
                fileName = sprintf('%s_%s_%s_psd_%s_%s_%s.mat', sID, str_ev, options.epochWdW,...
                    str_model, str_trGroup, str_style);
            end
        else
            error('Unknown case')
        end
    otherwise
        error('Not coded yet');
end
end