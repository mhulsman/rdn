function p = mt_sel_arrays(probes,arraynrs)

p = probes;

%2-dims
dim2 = {'pm','mm','annotations','bg_estimate','array_factors','resids','probe_weights','seq_factors',...
        'ampseq_correction','seq_correction','amp_correction','image_factors','qq_factors','annot_data'};
for i = 1:length(dim2)
    if(isfield(probes,dim2{i}))
        f = getfield(probes,dim2{i});
        f = f(arraynrs,:);
        p = setfield(p,dim2{i},f);
    end;
end;

dim1 = {'array_names','array_filenames','donor_names','donor_index','seq_norm','conc_array_idx','array_ind'};
for i = 1:length(dim1)
    if(isfield(probes,dim1{i}))
        f = getfield(probes,dim1{i});
        f = f(arraynrs);
        p =  setfield(p,dim1{i},f);
    end;
end;

