function dat = cummin(dat)

minval = inf;
for i = 1:length(dat)
    if(dat(i) < minval)
        minval = dat(i);
    else
        dat(i) = minval;
    end;
end;
