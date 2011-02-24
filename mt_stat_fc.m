%fold change. for 2-class problems
function [res,xvar] = mt_stat_fc(signal,labels,xvar)
   res = (mean(signal(~~labels(:,2),:),1) - mean(signal(~~labels(:,1),:),1))';
   xvar = [];
   return;

