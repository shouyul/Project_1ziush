function [ AUC ] = AUC_manual(fpr, tpr)
h=tpr(1:end-1)+tpr(2:end);
w=diff(fpr);
AUC=sum(h.*w)/2;
end

