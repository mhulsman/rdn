%MT_CREATE_AMP_MODEL - Determines the amplification model parameters
%
%  AMPM = MT_CREATE_AMP_MODEL(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate the amp. model params
%   VARARGIN        'splines'   : Use splines instead of dinucleotide counts
%
% OUTPUT
%   AMPM		    Matrix (nprobe * nvars) with hyb. params
%
% DESCRIPTION
% Determines parameters for the amplification model. Counts the number of all 
% dinucleotides from the 3'end of the probeset
%
% SEE ALSO
% MT_COR_HYBAMP

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function D = mt_create_amp_model(probes,varargin)

degree = 2;
points = 15;

basis = 1;
for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'splines',
            use_splines = 1;
      end;
   end;
end;

nprobe = size(probes.pm,2);
narray = size(probes.pm,1);
ngene = length(probes.name);

seqpos = zeros(nprobe,1);
fprintf(1,'Finding probe positions');
for i = 1:nprobe
   if(mod(i,10000) == 0)
     fprintf(1,'.');
   end;
   gnr = probes.ind(i);
   tseq = fliplr(lower(probes.gene_sequence{gnr}));
   if(~isempty(tseq))
      pos = strfind(tseq,fliplr(lower(probes.sequence{i})));
   else
      pos = [];
   end;

   if(length(pos) == 0)
      pos = 1;
   else
      pos = pos + 12; %middle of probe
   end;

   seqpos(i) = pos(1);
end;
fprintf(1,'\n');

partseq = {};

first_probe = zeros(nprobe,1);
last_probe = zeros(nprobe,1);

fprintf(1,'Extracting partial sequences');
for i = 1:ngene
   if(mod(i,1000) == 0)
      fprintf(1,'.');
   end;
   t = find(probes.ind == i);

   [dummy,fi] = max(seqpos(t));
   first_probe(i) = t(fi);
   [dummy,li] = min(seqpos(t));
   last_probe(i) = t(li);


   tseq = fliplr(lower(probes.gene_sequence{i}));
   if(~isempty(tseq))
      partseq{i} = tseq(seqpos(last_probe(i)):seqpos(first_probe(i)));
   end;
end;
fprintf(1,'\n');
   
dna_letters = ['A','C','G','T'];
seq_letter = {};

fprintf(1,'Creating nucleotide-specific boolean sequences');
for i = 1:length(dna_letters)
   fprintf(1,'.');
   tmp = {};
   for j = 1:ngene
      tmp{j} = (lower(partseq{j}) == lower(dna_letters(i)));
   end;
   seq_letter{i} = tmp;
end;
fprintf(1,'\n');



D = zeros(nprobe,length(dna_letters).^2, 'single');
fpos = seqpos(last_probe(probes.ind));
seqpos = seqpos - fpos;
counter = 0;

%fprintf(1,'Determining nucleotide counts');
%for i = 1:length(dna_letters)
%   fprintf(1,'.');
%   counter = counter + 1;
%   tx = seq_letter{i};
%   for j = 1:ngene
%      tx{j} = cumsum(tx{j});
%   end;
%   for j = 1:nprobe
%      k = tx{probes.ind(j)};
%      if(isempty(k) | seqpos(j) == 0)
%         D(j,counter) = 0;
%      else
%         D(j,counter) = k(seqpos(j));
%      end;
%   end;
%end;
%fprintf(1,'\n');
fprintf(1,'Determining dinucleotide counts');
for i = 1:length(dna_letters)
   for k = 1:length(dna_letters)
      fprintf(1,'.');
      counter = counter + 1;
      tx1 = seq_letter{i};
      tx2 = seq_letter{k};
      tx = {};
      for j = 1:ngene
         k1 = tx1{j};
         k2 = tx2{j};
         tx{j} = cumsum(k1(1:(end-1)) & k2(2:end));
         tx{j} = [tx{j} max(tx{j})];
      end;
      for j = 1:nprobe
         k = tx{probes.ind(j)};
         if(isempty(k) | seqpos(j) == 0)
            D(j,counter) = 0;
         else
            D(j,counter) = k(seqpos(j));
         end;
      end;
   end;
end;
fprintf(1,'\n');

%fprintf(1,'Determining trinucleotide counts');
%for i = 1:length(dna_letters)
%   for w = 1:length(dna_letters)
%      for q = 1:length(dna_letters)
%         fprintf(1,'.');
%         counter = counter + 1;
%         tx1 = seq_letter{i};
%         tx2 = seq_letter{w};
%         tx3 = seq_letter{q};
%         tx = {};
%         for j = 1:ngene
%            k1 = tx1{j};
%            k2 = tx2{j};
%            k3 = tx3{j};
%            tx{j} = cumsum(k1(1:(end-2)) & k2(2:(end - 1)) & k3(3:end));
%            tx{j} = [tx{j} max(tx{j}) max(tx{j})];
%         end;
%         for j = 1:nprobe
%            k = tx{probes.ind(j)};
%            if(isempty(k) | seqpos(j) == 0)
%               D(j,counter) = 0;
%            else
%               D(j,counter) = k(seqpos(j));
%            end;
%         end;
%      end;
%   end;
%end;
%fprintf(1,'\n');

if(exist('use_splines'))
   fprintf(1,'Creating splines');
   D2 = zeros(nprobe,size(D,2) * points, 'single');

   rmax = max(D);

   for i = 1:size(D,2)
      range = [0 rmax];
      for j = range(1):range(2)
         t = find(D(:,i) == j);
         D2(t,(1 + ((i-1) * points)):(i * points)) = repmat(mt_spline_support(points,j,degree,range),length(t),1);
      end;
   end;
   D = D2;
end;
fprintf(1,'\n');
return;


