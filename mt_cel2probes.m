%MT_CEL2PROBES - Convert one or more Affymetrix CEL files to probe intensities
%
%  PROBES = MT_CEL2PROBES(CEL,CDF,GIN,PA,SEQ,RMOUTLIERS,PROGRESS)
%
% INPUT
%		CEL				CEL structure or cell array of CEL-structures
%		CDF				CDF structure
%       GIN				GIN structure (optional)
%       PA              Probe annotation structure (optional)
%       SEQ             SEQ structure (optional)
%		RMOUTLIERS	    Flag indicating whether outliers should be removed (default: 0)
%		PROGRESS		Flag indicating whether to print progress messages (default: 1)
%
% OUTPUT
%		PROBES			Structure containing matrices PM and MM of size (#CELs x 
%									#probes), gene NAMEs and indices IND of probes into genes
% DESCRIPTION
% Converts a number of Affymetrix CEL files to perfect match (PM) and mismatch
% (MM) intensity values. The data is stored in a self-explanatory structure.
%
% SEE ALSO
% MT_READCEL, MT_READCDF, MT_READGIN, MT_READPROBE_ANNOT, MT_READSEQ

% REVISION HISTORY
%   %M1/M2/M3 Marc Hulsman, Add GIN,PA,SEQ inclusion

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_cel2probes (cel,cdf,gin,probeannot,seq,rmoutliers,progress)

	if (nargin < 7), progress   = 1; end;
	if (nargin < 6), rmoutliers = 0; end;
	if (nargin < 5), seq = []; end;
	if (nargin < 4), probeannot = []; end;
	if (nargin < 3), gin        = []; end;

	if (progress), fprintf(1,'Converting\n'); end;

	% Convert single CEL to cell array.

	if (length(cel)==1), cel = {cel}; end;

	nprobes = 0;
	for c = 1:length(cel)
		nprobes = max(cel{c}.nprobes,nprobes);
	end;

	% Pre-allocate memory.

  npm = 0;
  nmm = 0;
  for i = 1:cdf.nunits
    for j = 1:cdf.gene{i}.nblocks
      npm = npm + size(cdf.gene{i}.block{j}.pm,1);
      if(isfield(cdf.gene{i}.block{j}, 'mm'))
          nmm = nmm + size(cdf.gene{i}.block{j}.mm,1);
      end;
    end;
  end;

	%probes.array = cdf.name;
	probes.nrows = cdf.rows;
	probes.ncols = cdf.cols;

	probes.ind = zeros(npm,1,'uint32');
    probes.indmm = zeros(nmm, 1, 'uint32');
    probes.sequence = cell(npm, 1);
    probes.position = zeros(npm, 1, 'uint32');

	pm_pos = zeros(npm,2, 'uint32');
	mm_pos = zeros(nmm,2, 'uint32');

	% Collate outlier information, if requested.
	
	noutliers = 0;
	for c = 1:length(cel)
		noutliers = noutliers + cel{c}.noutliers;
	end;
	outlier = spalloc(cel{1}.rows,cel{1}.cols,noutliers);
	for c = 1:length(cel)
		outlier = outlier + cel{c}.outlier;
	end;

	if (progress), fprintf(1,'  CDF ',c); end;

	kpm = 0; kmm = 0; l = 0; oi = 0;

	for i = 1:cdf.nunits

		if (progress & (mod(i,1000)==0)), fprintf(1,'.'); end;

		for j = 1:cdf.gene{i}.nblocks

         l = l + 1; 

   %  		nprobes        = length(cdf.gene{i}.block{j}.pm);
         nprobespm        = size(cdf.gene{i}.block{j}.pm,1);
         if(isfield(cdf.gene{i}.block{j}, 'mm'))
            nprobesmm = size(cdf.gene{i}.block{j}.mm, 1);
         else
            nprobesmm = 0;
         end;
         
         probes.name{l} = cdf.gene{i}.block{j}.name;
         
         %M1
         if(~isempty(probeannot))
            pos = cdf.gene{i}.block{j}.pm;
            paidx = probeannot.probeind(sub2ind(size(probeannot.probeind),uint32(pos(:,1)),uint32(pos(:,2))));
            if(any(~paidx))
               probes.position(kpm + (1:nprobespm)) = 0;
               for q = 1:length(paidx)
                  probes.sequence{kpm + q} = '';
               end;
            else
               probes.position(kpm + (1:nprobespm)) = probeannot.position(paidx);
               for q = 1:length(paidx)
                  probes.sequence{kpm + q} = probeannot.sequence{paidx(q)};
               end;
            end;
         end;

         pm_pos(kpm +(1:nprobespm),:) = cdf.gene{i}.block{j}.pm;
         if(nprobesmm > 0)
             mm_pos(kmm +(1:nprobesmm),:) = cdf.gene{i}.block{j}.mm;
         end;
         probes.ind(kpm +(1:nprobespm)) = l;
         probes.indmm(kmm + (1:nprobesmm)) = l;
         kpm = kpm + nprobespm;
         kmm = kmm + nprobesmm;

		end;
	end;
   
   
   [names,i] = sort(lower(probes.name));
   idouble = sortrows([(1:length(names))',i'],2);
   i_rev = idouble(:,1);

   %M2
   if(~isempty(gin))
      [gin_names,g_i] = sort(lower(gin.name));
      if(~(sum(strcmp(gin_names,names)) == length(names)  && sum(strcmp(names,gin_names)) == length(names)))
         fprintf(1,'ERROR: GIN name list and CDF name list do not match\n');
      else
        probes.desc = gin.desc(g_i(i_rev));
      end;
   end;
   
   %M3
   if(~isempty(seq))
      [seq_names,s_i] = sort(lower(seq.name));
      if(~(sum(strcmp(seq_names,names)) == length(names) && sum(strcmp(names,seq_names)) == length(names)))
         fprintf(1,'ERROR: SEQ name list and CDF name list do not match\n');
      else
         probes.gene_sequence = seq.sequence(s_i(i_rev));
      end;
   end;


	fprintf (1,'\n');

	pm_pos = pm_pos(1:kpm,:); mm_pos = mm_pos(1:kmm,:);

	probes.pm_pos = pm_pos';
	probes.mm_pos = mm_pos';

	pm_pos = pm_pos(1:kpm,1) + cdf.cols*(pm_pos(1:kpm,2)-1);
	mm_pos = mm_pos(1:kmm,1) + cdf.cols*(mm_pos(1:kmm,2)-1);

	probes.pm = zeros(length(cel),size(pm_pos,1), 'single');
    probes.mm = zeros(length(cel),size(mm_pos,1), 'single');

    probes.array_names = {};

	for c = 1:length(cel)
		if (progress), fprintf(1,'  CEL %d ',c); end;
      probes.array_names{c} = cel{c}.name;
      probes.array_filenames{c} = cel{c}.filename;
      probes.pm(c,:) = cel{c}.mean(pm_pos);
 	  probes.mm(c,:) = cel{c}.mean(mm_pos);
	  fprintf(1,'\n');
	end;

    probes.ind = probes.ind(1:kpm);
    probes.indmm = probes.indmm(1:kmm);

	% Assumption: no probesets disappear completely.
	probes.ol  = find(outlier(pm_pos))';

	if (rmoutliers)
		probes = mt_probes_rmoutliers(probes);									
	end;

	if (progress), fprintf(1,'\n'); end;

return

