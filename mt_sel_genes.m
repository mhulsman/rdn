function p = select_genes(probes,genenrs)

p = probes;
fn = fieldnames(probes);
gene_dim1 = {'name','desc','chrom_loc','seqlength','gene_names','gene_description','chrom_loc_det','chrom_loc','gene_sequence',...
        'overall_factors','s_excl'};

for i = 1:length(gene_dim1)
    if(isfield(probes,gene_dim1{i}))
        f = getfield(probes,gene_dim1{i});
        f = f(genenrs);
        p = setfield(p,gene_dim1{i},f);
    end;
end;

gene_dim2 = {'batch_factors','array_factors'};
for i = 1:length(gene_dim2)
    if(isfield(probes,gene_dim2{i}))
        f = getfield(probes,gene_dim2{i});
        f = f(:,genenrs);
        p = setfield(p,gene_dim2{i},f);
    end;
end;





if(isfield(probes,'ind'))
    %renumber genes
    nprobepm = size(probes.ind,2);
    
    z = zeros(length(probes.name),1);
    translate = genenrs;
    z(translate) = 1:length(translate);
    
    %probe ids
    p.ind = z(probes.ind)';
    pids = p.ind ~= 0;
    p.ind = p.ind(pids);
   
    if(isfield(probes,'indmm'))
        nprobemm = size(probes.indmm, 2);
        p.indmm = z(probes.indmm)';
        mids = p.indmm ~= 0;
        p.indmm = p.indmm(mids);
    else
        mids = pids;
    end;
    probe_dim1 = {'position','sequence','lastpos','probe_factors',...
                  'seqbg','removed_factors'};
    for i = 1:length(probe_dim1)
        if(isfield(probes,probe_dim1{i}))
            f = getfield(probes,probe_dim1{i});
            f = f(pids);
            p = setfield(p,probe_dim1{i},f);
        end;
    end;
    probe_dim2_mm = {'mm','mm_pos'};
    for i = 1:length(probe_dim2_mm)
        if(isfield(probes, probe_dim2_mm{i}))
            f = getfield(probes,probe_dim2_mm{i});
            f = f(:,mids);
            p = setfield(p,probe_dim2_mm{i},f);
        end;
    end;

    probe_dim2 = {'pm','resids','probe_weights','seq_correction',...
                  'ampseq_correction','amp_correction','image_factors','qq_factors',...
                  'pm_pos','resids'};
    for i = 1:length(probe_dim2)
        if(isfield(probes, probe_dim2{i}))
            f = getfield(probes,probe_dim2{i});
            f = f(:,pids);
            p = setfield(p,probe_dim2{i},f);
        end;
    end;
end;


if(isfield(probes,'conc_idx'))
  tmp = z(probes.spike_idx);
  spike_keep = find(tmp > 0);
  p.spike_idx = tmp(spike_keep);
  p.conc_probeset_idx = probes.conc_probeset_idx(spike_keep);
  p.spike_name = probes.spike_name(spike_keep);
end;
