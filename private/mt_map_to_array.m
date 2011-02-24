function map = mt_map_to_array(probes,data)

   [narray,nprobe] = size(probes.pm);

   values = zeros(probes.nrows,probes.ncols);

   pm_pos = probes.pm_pos(1,:)' + probes.ncols*(probes.pm_pos(2,:)' - 1); 
   mm_pos = probes.mm_pos(1,:)' + probes.ncols*(probes.mm_pos(2,:)' - 1); 

   [narrays,nprobe] = size(probes.pm);

   values(pm_pos) = data;
   if(length(mm_pos) == length(data))
       values(mm_pos) = data;
   end;
   map  = values';

