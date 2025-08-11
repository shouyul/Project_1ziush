labels = {'Avec Objet', 'Sans Objet'};
days = {'2023-09-18'};
N_blocks = 10;
N_obj = 5;
N_withoutobj = 5;
N_tot = N_obj + N_withoutobj;

for d = 1:numel(days)
    fid = fopen(sprintf('trials_%s.txt',days{d}), 'wt');
    for b = 1:N_blocks
        count_obj=0;
        count_withoutobj=0;
        for tr = 1:N_tot
            if mod(tr,3) == 0
                nextinput = '%i. %s\n';
            else
                nextinput = '%i. %s\t\t';
            end
            
            tr_dayStart = tr + (b-1)*N_tot;
            if count_obj == N_obj
                count_withoutobj=count_withoutobj+1;
                fprintf(fid,nextinput,tr_dayStart,labels{2});
                continue
            end
            if count_withoutobj == N_withoutobj
                count_obj=count_obj+1;
                fprintf(fid,nextinput,tr_dayStart,labels{1});
                continue
            end
            
            a = rand();
            if a < N_obj/N_tot
                count_obj=count_obj+1;
                fprintf(fid,nextinput,tr_dayStart,labels{1});
            else
                count_withoutobj=count_withoutobj+1;
                fprintf(fid,nextinput,tr_dayStart,labels{2});
            end
        end        
        fprintf(fid,'PAUSE\n\n');
    end
    fclose(fid);
end