function nprobes = mt_countprobes(probes)

nprobes = zeros(length(probes.name),1);

for i = 1:length(probes.ind)
    r = probes.ind(i);
    nprobes(r) = nprobes(r) + 1;
end;
