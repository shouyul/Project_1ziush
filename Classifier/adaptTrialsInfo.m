function TrialsInfo2 = adaptTrialsInfo(TrialsInfo, TumbVis, nb_chunks)
% Adds Chunk and TumbVis info to the TrialsInfo table
if isempty(TumbVis)
    if nb_chunks == 1
        Chunk = ones(size(TrialsInfo,1),1);
        TrialsInfo2 = addvars(TrialsInfo, Chunk,...
            'NewVariableNames', {'Chunk'});
    else
        Chunk = ones(size(TrialsInfo,1),1);
        TrialsInfo2 = addvars(TrialsInfo, Chunk,...
            'NewVariableNames', {'Chunk'});
        for ch = 2:nb_chunks
            TrialsInfo2 = cat(1,TrialsInfo2,addvars(TrialsInfo, Chunk.*ch,...
                'NewVariableNames', {'Chunk'}));
        end
    end
else
    if nb_chunks == 1
        Chunk = ones(size(TrialsInfo,1),1);
        TrialsInfo2 = addvars(TrialsInfo, Chunk, TumbVis(:,1),...
            'NewVariableNames', {'Chunk', 'TumbVis'});
    else
        Chunk = ones(size(TrialsInfo,1),1);
        TrialsInfo2 = addvars(TrialsInfo, Chunk, TumbVis(:,1),...
            'NewVariableNames', {'Chunk', 'TumbVis'});
        for ch = 2:nb_chunks
            TrialsInfo2 = cat(1,TrialsInfo2,addvars(TrialsInfo, Chunk.*ch, TumbVis(:,ch),...
                'NewVariableNames', {'Chunk', 'TumbVis'}));
        end
    end
end
TrialsInfo2 = sortrows(TrialsInfo2,'Trial');
TrialsInfo2 = sortrows(TrialsInfo2,'Block');
end