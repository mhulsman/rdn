function medians  = fast_median(data,dim)
   if(length(size(data)) > 3)
      error('fast_median is only tested up to 3 dimensions')
   end;

   if(nargin < 2)
      dim = 1;
   end;

   data = sort(data,dim);
   s = size(data);
   nCompare = s(dim);
   half_idx = floor(nCompare / 2);
   cs = cumprod(s);

   if(dim == 1)
      steps = (1:cs(dim):cs(end)) + half_idx;
   else
      if(length(s) == 2 || dim == 3)
        steps = half_idx * cs(dim -1) + (1:cs(dim - 1));
      else
         if(length(s) == 3)
            steps =repmat(half_idx * cs(dim -1) + (1:cs(dim - 1))',1,s(3)) + repmat((0:cs(2):(cs(end)-1)),cs(dim-1),1);
         end;
      end;
   end;
   
   medians = data(steps);

   if(half_idx * 2 == nCompare)
      if(dim == 1)
         medians = (medians + data(steps - 1)) / 2;
      else
         medians = (medians + data(steps - cs(dim - 1)))/2;
      end;
   end;

   s(dim) = 1;

   medians = reshape(medians,s);
