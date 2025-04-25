function metadata = load_chanInfo(cfg,metadata,files)
%%
ch = metadata.ch;

fprintf('%s%s/ses-%s/ieeg/%s_ses-%s_electrodes.tsv is used.\n',cfg(2).proj_dirinput,metadata.sub_label,metadata.ses_name,metadata.sub_label,metadata.ses_name)
elecName = fullfile(files(1).folder, '/',files(1).name);
cc_elecs = readtable(elecName,'FileType','text','Delimiter','\t');

% pre-allocation
elec_incl       = false(size(ch));
elec_soz        = false(size(ch));
elec_pp        = false(size(ch));
elec_pr        = false(size(ch));
elec_silicon    = false(size(ch));
elec_resected   = false(size(ch));
elec_edge       = false(size(ch));
elec_screw      = false(size(ch));
elec_wm         = false(size(ch));
elec_gm         = false(size(ch));
elec_csf        = false(size(ch));
elec_amyg       = false(size(ch));
elec_hipp       = false(size(ch));
elec_lesion     = false(size(ch));
elec_gliosis    = false(size(ch));
elec_strip      = false(size(ch));
elec_grid       = false(size(ch));
elec_seeg       = false(size(ch));
elec_depth      = false(size(ch));

elec_name       = cc_elecs.name;

% excluding electrodes other (so no grid, strip, depth)
for i=1:size(ch,1)
    
    stimname = regexp(ch{i},'[a-z_A-Z]*','match');
    stimnum = regexp(ch{i},'[0-9]*','match');
    if ~isempty(stimname) && ~isempty(stimnum)
        test1 = sprintf('%s0%d',stimname{:}, str2double(stimnum{:}));
        test2 = sprintf('%s%d',stimname{:}, str2double(stimnum{:}));
    else
        test1 = ch{i};
        test2 = ch{i};
    end
    
    idx = false(size(ch));
    
    % if channel name is present in equal writing in electrodes.tsv
    if sum(cellfun(@(x) strcmpi(x,ch{i}),elec_name)) == 1
    
        idx = cellfun(@(x) strcmpi(x,ch{i}),elec_name);
        
        % if channelname with extra 0 is present in electrodes.tsv, but not
        % present without extra 0
    elseif sum(cellfun(@(x) strcmpi(x,test1),elec_name)) >0 && sum(cellfun(@(x) strcmpi(x,test2),elec_name)) == 0
        
        if sum(cellfun(@(x) strcmpi(x,test1),elec_name)) ==1 % if electrode name is once in electrodes list, ignore case
            idx = cellfun(@(x) strcmpi(x,test1),elec_name);
        elseif sum(cellfun(@(x) strcmp(x,test2),elec_name)) == 1 % if electrode name is twice in electrodes list, do not ignore case
            idx = cellfun(@(x) strcmp(x,test2),elec_name);
        end
        
        % if channelname is present without extra 0 but not with extra 0 in
        % electrodes.tsv
    elseif sum(cellfun(@(x) strcmpi(x,test2),elec_name)) >0 && sum(cellfun(@(x) strcmpi(x,test1),elec_name)) == 0
        
        if sum(cellfun(@(x) strcmpi(x,test2),elec_name)) ==1 % if electrode name is once in electrodes list, ignore case
            idx = cellfun(@(x) strcmpi(x,test2),elec_name);
        elseif sum(cellfun(@(x) strcmp(x,test2),elec_name)) == 1 % if electrode name is twice in electrodes list, do not ignore case
            idx = cellfun(@(x) strcmp(x,test2),elec_name);
        end
        
    end
    
    % fill in the properties found in electrodes.tsv for specific channel
    if sum(idx) == 1
        
        elec_incl(i)        = ~strcmp(cc_elecs.group{idx},'other');
        elec_soz(i)         = strcmp(cc_elecs.soz{idx},'yes');
        elec_pp(i)         = strcmp(cc_elecs.pp{idx},'yes');
        elec_pr(i)         = strcmp(cc_elecs.pr{idx},'yes');
        elec_silicon(i)     = strcmp(cc_elecs.silicon{idx},'yes');
        elec_resected(i)    = strcmp(cc_elecs.resected{idx},'yes');
        elec_edge(i)        = strcmp(cc_elecs.edge{idx},'yes');

        if strcmp(cc_elecs.group{idx},'grid')
            elec_grid(i) = true;
        elseif strcmp(cc_elecs.group{idx},'strip')
            elec_strip(i) = true;
        elseif strcmp(cc_elecs.group{idx},'depth')
            elec_depth(i) = true;
            elec_seeg(i) = true;
        end

        if contains(metadata.format_info,'seeg')
            elec_screw(i)   = strcmp(cc_elecs.screw{idx},'yes');
            elec_wm(i)      = strcmp(cc_elecs.whitematter{idx},'yes');
            elec_gm(i)      = strcmp(cc_elecs.graymatter{idx},'yes');
            elec_csf(i)     = strcmp(cc_elecs.csf{idx},'yes');
            elec_amyg(i)    = strcmp(cc_elecs.amygdala{idx},'yes');
            elec_hipp(i)    = strcmp(cc_elecs.hippocampus{idx},'yes');
            elec_lesion(i)  = strcmp(cc_elecs.lesion{idx},'yes');
            elec_gliosis(i) = strcmp(cc_elecs.gliosis{idx},'yes');
            
        end        
        
    else
        elec_incl(i)        = false;
        elec_soz(i)         = false;
        elec_pp(i)         = false;
        elec_pr(i)         = false;
        elec_silicon(i)     = false;
        elec_resected(i)    = false;
        elec_edge(i)        = false;
        
        if contains(metadata.format_info,'seeg')
            elec_screw(i)   = false;
            elec_wm(i)      = false;
            elec_gm(i)      = false;
            elec_csf(i)     = false;
            elec_amyg(i)    = false;
            elec_hipp(i)    = false;
            elec_lesion(i)  = false;
            elec_gliosis(i) = false;
            
        end
    end
end

metadata.ch2use_included    = logical(elec_incl);
metadata.ch2use_silicon     = logical(elec_silicon);
metadata.ch2use_soz         = logical(elec_soz);
metadata.ch2use_pp         = logical(elec_pp);
metadata.ch2use_pr         = logical(elec_pr);
metadata.ch2use_resected    = logical(elec_resected);
metadata.ch2use_edge        = logical(elec_edge);

metadata.idx_grid   = logical(elec_grid);
metadata.idx_strip  = logical(elec_strip);
metadata.idx_depth  = logical(elec_depth);
metadata.idx_seeg   = logical(elec_seeg);

if contains(metadata.format_info,'seeg')
    
    metadata.ch2use_screw   = logical(elec_screw);
    metadata.ch2use_wm      = logical(elec_wm);
    metadata.ch2use_gm      = logical(elec_gm);
    metadata.ch2use_csf     = logical(elec_csf);
    metadata.ch2use_amyg    = logical(elec_amyg);
    metadata.ch2use_hipp    = logical(elec_hipp);
    metadata.ch2use_lesion  = logical(elec_lesion);
    metadata.ch2use_gliosis = logical(elec_gliosis);
    
end

end
