%for continuous vectors. Spearman rank correlation
function [sval,s0] = mt_stat_cor(signal,labels,s0)
   n_arrays = size(signal,1);
   sval = corr(signal,labels,'type','Spearman');

