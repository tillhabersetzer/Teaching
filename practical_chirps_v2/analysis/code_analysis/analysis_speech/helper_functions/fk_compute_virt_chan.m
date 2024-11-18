function data=fk_compute_virt_chan(cfg,data)

    beamfilts = cat(1,cfg.filter{:});
    
    nrpt = length(data.trial);
    
    fprintf(['using precomputed filters \n']);
    ft_progress('init', 'text');
    
    datatemp = data;
    
    for i = 1:nrpt
        
        ft_progress(i/nrpt, 'scanning repetition %d of %d', i,nrpt);
        
        data.trial{i} = beamfilts*datatemp.trial{i};
        
    end
    
    ft_progress('close');
    
    data.label = cellstr(num2str([1:size(beamfilts,1)]'));
        
    data=rmfield(data, 'cfg');
    data=rmfield(data, 'grad');
    
    if isfield(cfg, 'outputfilename')
        save(cfg.outputfilename, 'data');
    end
end