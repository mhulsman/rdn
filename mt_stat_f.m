%for continous vectors. F-test
function [sval,s0] = mt_stat_f(signal,labels,s0)
   n_arrays = size(signal,1);
   

   C_full = [ones(n_arrays,1) labels];
   C_simple = ones(n_arrays,1);

   p_full = C_full \ signal;
   p_simple = C_simple \ signal;


   res_full = sum((C_full * p_full - signal).^2);
   res_simple = sum((C_simple * p_simple - signal).^2);
   
   
   sval = ((res_full - res_simple) ./ (res_full ./ (n_arrays - size(C_full,2))));

   sval = sval' .* sign(p_full(2,:))';




