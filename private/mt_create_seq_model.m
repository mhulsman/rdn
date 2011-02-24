%MT_CREATE_SEQ_MODEL - Determines the hybridization model parameters
%
%  SEQM = MT_CREATE_SEQ_MODEL(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate the hyb. model params
%   VARARGIN        'basis',x   : Create model for single nucleotides (1), dinucleotides (2), 
%                                 trinucleotides(3). (default: 2)
%
% OUTPUT
%   SEQM		    Matrix (nprobe * nvars) with hyb. params
%
% DESCRIPTION
% Determines parameters for the hybridization model. Uses b-splines to
% reduce overfitting and number of variables. Variables are calculated
% for ((single|bi|tri)nucleotides w.r.t to probe position, as well as
% nucleotide (pair) counts for the whole probe sequence and the left and right part. 
%
% SEE ALSO
% MT_BG_EST, MT_COR_HYBAMP

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function C = mt_create_seq_model(probes,varargin)

basis = 1;
for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'basis',
            i = i + 1;
            basis = varargin{i};
      end;
   end;
end;

nprobe = size(probes.pm,2);
narray = size(probes.pm,1);

q = upper(char(probes.sequence));

nletters = size(q,2);

dna_letters = ['A','C','G','T'];
probseq = zeros(nprobe,length(dna_letters) * nletters);

nucleo_count = {};
for i = 1:3
   nucleo_count{i} = zeros(nprobe,4);
end;
for i = 1:length(dna_letters)
   t =  (q == dna_letters(i));
   nucleo_count{1}(:,i) = sum(t,2);
   nucleo_count{2}(:,i) = sum(t(:,1:floor(nletters/2)),2);
   nucleo_count{3}(:,i) = sum(t(:,ceil(nletters/2):nletters),2);
   probseq(:,(1 + (i-1) * nletters):(i*nletters)) = t;
end;

points = 5;
degree = 3;

if(basis == 1)
   C = zeros(nprobe,points * (4 + 30),'single');
   %walk over dna letters and probe pos
   for i = 1:length(dna_letters)
      for j = 1:nletters
         %find probes that have a specific letter at specific pos. 
         t = find(probseq(:,(i -1) * nletters + j));
         %calculate b-spline factors and add them to seqmodel param matrix
         C(t,(1 + (i-1) * points):(i * points)) = C(t,(1 + (i-1) * points):(i * points)) + repmat(mt_spline_support(points,j,degree,[1 nletters]),length(t),1);
      end;
   end;
   C_start = (points * 4);
else 
   if(basis == 2)
      C = zeros(nprobe,points * (16 + 30),'single');
      cset = 0;
      %walk over dinucleotides
      for i = 1:length(dna_letters)
         for j = 1:length(dna_letters)
            cset = cset + 1;
            combine = probseq(:,(1 + (i-1) * nletters):((i*nletters) - 1)) & probseq(:,(2 + (j-1) * nletters):(j*nletters));
            %walk over probe position
            for k = 1:(nletters-1)
               %find probes that have a specific letter at specific pos. 
               t = find(combine(:,k));
               %calculate b-spline factors and add them to seqmodel param matrix
               C(t,(1 + (cset-1) * points):(cset * points)) = C(t,(1 + (cset-1) * points):(cset * points)) + repmat(mt_spline_support(points,k,degree,[1 (nletters-1)]),length(t),1);
            end;      
         end;
      end;
      C_start = (points * 16);
    else
      C = zeros(nprobe,points * (64 + 30),'single');
      cset = 0;
      %walk over trinucleotides
      for i = 1:length(dna_letters)
         for j = 1:length(dna_letters)
            for g = 1:length(dna_letters)
               cset = cset + 1;
               combine = probseq(:,(1 + (i-1) * nletters):((i*nletters) - 2)) & probseq(:,(2 + (j-1) * nletters):((j*nletters) - 1)) & probseq(:,(3 + (g-1) * nletters):((g*nletters)));
               %walk over probe position
               for k = 1:(nletters-2)
               %find probes that have a specific letter at specific pos. 
                  t = find(combine(:,k));
               %calculate b-spline factors and add them to seqmodel param matrix
                  C(t,(1 + (cset-1) * points):(cset * points)) = C(t,(1 + (cset-1) * points):(cset * points)) + repmat(mt_spline_support(points,k,degree,[1 (nletters-2)]),length(t),1);
               end;     
            end;
         end;
      end;
      C_start = (points * 64);
    end;
end;


%walk over all, left, right probe sequence part
for q = 1:3
   
   %walk over nucleotide letters
   for i = 1:4
      %determine range for b-splines
      range = [min(nucleo_count{q}(:,i)) max(nucleo_count{q}(:,i))];
      %for each integer within range, calculate b-spline factors
      for j = range(1):range(2)
         t = find(nucleo_count{q}(:,i) == j);
         C(t,C_start + (i-1) * points +  (1:points)') = C(t,C_start + (i - 1) * points + (1:points)') + repmat(mt_spline_support(points,j,degree,range),length(t),1);
      end;
   end;
   C_start = C_start + (points * 4);

   counter = 0;
   %walk over nucleotide pairs
   for i = 1:length(dna_letters)
      for k = (i+1):length(dna_letters)
        counter = counter + 1;
        count = nucleo_count{q}(:,i) + nucleo_count{q}(:,k); 
        
        %determine range for b-splines
        range = [min(count) max(count)];
        %for each integer within range, calculate b-spline factors
        for j = range(1):range(2)
            t = find(count == j);
            C(t,C_start + (counter-1) * points +  (1:points)') = C(t,C_start + (counter - 1) * points + (1:points)') + repmat(mt_spline_support(points,j,degree,range),length(t),1);
        end;
      end;
   end;
   C_start = C_start + (points * 6);
end;

return;


