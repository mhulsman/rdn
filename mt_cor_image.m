%MT_COR_IMAGE - Calculates array location based correction factors
%
%  PROBES = MT_COR_IMAGE(PROBES,FILTER_SIZE,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate image cor factors
%   FILTER_SIZE     Size of median box filter (default: 9)
%   VARARGIN        'use_optical' : Use optical background information from mt_bg_est
%
% OUTPUT
%   PROBES		    Probes structure with image cor. factors added
%
% DESCRIPTION
% Calculates image correction factors (after application of ealier 
% corrections) and adds them to the probes structure as field image_factors.   
%
% SEE ALSO
% MT_BG_EST, MT_COR_QQ, MT_COR_HYBAMP, MT_REAL_SIGNAL
%
% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_cor_image(probes,filter_size,varargin)

for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'use_optical',
            use_optical = 1;
      end;
   end;
end;

[narray,nprobe] = size(probes.pm);
pad_size = floor(filter_size /2 );

signal = mt_real_signal(probes);
mt = median(signal);

if(~isfield(probes,'image_factors'))
   probes.image_factors = zeros(size(probes.pm));
end;

if(exist('use_optical'))
   o = probes.seqbg_factors(1);
   ssignal = 2.^signal - o;
   ssignal(ssignal<1) = 1;
   ssignal = log2(ssignal);
   mt = 2.^mt - o;
   mt(mt < 1) = 1;
   data =  ssignal - repmat(log2(mt),narray,1);
else
   data = signal - repmat(mt,narray,1);
end;


nrows = probes.nrows + pad_size * 2;
ncols = probes.ncols + pad_size * 2;

values = zeros(probes.nrows + pad_size * 2,probes.ncols + pad_size * 2, 'single');

pm_pos = int32(nrows * pad_size + pad_size + probes.pm_pos(1,:)' + ncols*(probes.pm_pos(2,:)' - 1)); 

%create array image with probe indexes
values(pm_pos) = 1:nprobe;

% median indexes (w.r.t. middle point)
t = repmat(int32((-pad_size:pad_size)'),1,filter_size) + repmat(int32((-pad_size:pad_size) * nrows),filter_size,1);

number_of_outer_chunks = round((nprobe * (filter_size.^2))/10.^7);
ochunk = 1:round(length(pm_pos)/number_of_outer_chunks):length(pm_pos);
if ochunk(end)~=length(pm_pos), ochunk = [ochunk length(pm_pos)]; end

median_resids = zeros(size(probes.pm), 'single');

for w = 1:(length(ochunk) - 1)
   mlength = fprintf(1,'(Chunk %d/%d), Calculating for median area: ',w,length(ochunk) - 1);

   idxs = repmat(t(:)',length(pm_pos(ochunk(w):ochunk(w+1))),1) + repmat(pm_pos(ochunk(w):ochunk(w+1)),1,numel(t));
   % convert to probe indexes
   res = values(idxs);

   % which elements have a real probe (idx > 0), count per pixel
   median_elements = sum((res > 0),2);
   u_median_elements = unique(median_elements);


   % for each number of elements, create idxs
   med_idxs = struct();
   for i = 1:length(u_median_elements)
      %which probes have this number of elements
      med_idxs(i).probe_idxs = int32(find(median_elements == u_median_elements(i)));
      
      % get probe idxs
      im_res = res(med_idxs(i).probe_idxs,:)';
      
      % reshape the non-0 idxs to use as index array
      med_idxs(i).median_idxs = reshape(int32(im_res(im_res > 0)),u_median_elements(i),length(med_idxs(i).probe_idxs));
   end;

   %calculate median-filtered resids for each number of elements per pixel
   for i = 1:length(med_idxs)
      olength = fprintf(1,num2str(size(med_idxs(i).median_idxs,1)));
         
      nprobes = length(med_idxs(i).probe_idxs);
      
      % if number of probes is too large, split into chunks
      number_of_chunks = ceil(nprobes / 1000);
      nchunk = 1:round(nprobes/number_of_chunks):nprobes;
      if nchunk(end)~=nprobes
         nchunk = [nchunk nprobes];
      end
      % for each chunk
      for j = 1:(length(nchunk) - 1)
         if(length(nchunk) > 2)
            nlength = fprintf(1,' (probes: %d/%d)',nchunk(j),nprobes);
         else
            nlength = 0;
         end;
         % idx sub part
         di = med_idxs(i).median_idxs(:,nchunk(j):nchunk(j+1));
         % size of idx
         q = size(di);
         % calculate the part of resids in this chunk
         tmp =reshape(data(:,di),[narray,q(1),q(2)]);
         median_resids(:,ochunk(w) - 1 + med_idxs(i).probe_idxs(nchunk(j):nchunk(j+1))) =squeeze(mt_fast_median(tmp,2));
         fprintf(1,repmat('\b',1,nlength));
      end;
      fprintf(1,repmat('\b',1,olength));

   end;
   fprintf(1,repmat('\b',1,mlength));

end;
fprintf(1,'Calculating median area: completed\n')
if(exist('use_optical'))
   o = probes.seqbg_factors(1);
   nsignal = log2(o + 2.^(ssignal - median_resids));
   probes.image_factors = probes.image_factors + (signal - nsignal);
else
   probes.image_factors = probes.image_factors + median_resids; 
end;

return;


