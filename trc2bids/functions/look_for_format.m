function metadata = look_for_format(metadata,format_idx,annots)

annots_format = annots(format_idx,2);

% putting format in order ECoG,strip,depth
annotsformatsplit = strsplit([annots_format{:}],{';','Format;'});
annotsformatsplit = annotsformatsplit(~cellfun(@isempty,annotsformatsplit));
ecogloc = find(cellfun('length',regexp(lower(annotsformatsplit),'ecog')) == 1);
striploc = find(cellfun('length',regexp(lower(annotsformatsplit),'strip')) == 1);
depthloc = find(cellfun('length',regexp(lower(annotsformatsplit),'depth')) == 1);
seegloc = find(cellfun('length',regexp(lower(annotsformatsplit),'seeg')) == 1);

locs_all = sort([ecogloc, striploc, depthloc,seegloc,size(annotsformatsplit,2)+1]);

% find number of letters in each channel
strSizeCh = NaN(size(metadata.ch));
for chan = 1:size(metadata.ch,1)
    strSizeCh(chan,1) = sum(isstrprop(metadata.ch{chan},'alpha'));
end

% ECoG
idx_grid = false(size(metadata.ch) );
if ~isempty(ecogloc)
    ecogformat = cell(1,size(ecogloc,2));
    for i=1:size(ecogloc,2)
        ecogformat{i} = annotsformatsplit(ecogloc(i)+1: locs_all(find(locs_all==ecogloc(i))+1)-1);
    end
    ecogformatall =  strcat([ecogformat{:}],';');
    ecogformatfin = ['ECoG;' ecogformatall{:}];
    
    % add in metadata which electrodes are in strip
    for n = 1:size(ecogformatall,2)
        
        size_gridpart = strlength(extractBefore(ecogformatall{n},'[')) ;     
        % find which electrodes have same name, and same size (to avoid FML
        % to be part of ML, because both contains ML (RESP0991))
        
        % so apparently, the number in Format is not similar to the number
        % of electrodes with this name. This is the case in RESP0588. In
        % this patient, there are unplugged electrodes named C29-32, so
        % these don't contain any data (only noise), but are still in the
        % electrode list. This does not mean that the number in Format is
        % incorrect, so we should look how to exclude this situation from
        % getting this error.
        idx_gridpart = ismember(strSizeCh, size_gridpart) == 1 & ...
            contains(metadata.ch,extractBefore(ecogformatall{n},'[')) & ...
            metadata.ch2use_included == 1;
                
        % calculate number of electrodes in format. This helps in both
        % C[4x8] and in C[4x8,4x8] (when there are two grids located next
        % to eachother, with the same naming) (see RESP0621)
        ecogformatparts = strsplit(ecogformatall{n},',');
        formatparts = NaN(size(ecogformatparts));
        for i = 1:size(ecogformatparts,2)
            if contains(ecogformatparts{i},{'[','x'}) && ~contains(ecogformatparts{i},{']'}) % like C[4x8 (as part of C[4x8,4x8])
                formatparts(i) = str2double(extractBetween(ecogformatparts{i},'[','x')) * str2double(extractAfter(ecogformatparts{i},'x'));
            elseif contains(ecogformatparts{i},{']','x'}) && ~contains(ecogformatparts{i},{'['}) % like 4x8] (as part of C[4x8,4x8])
                formatparts(i) = str2double(extractBefore(ecogformatparts{i},'x')) * str2double(extractBetween(ecogformatparts{i},'x',']'));
            elseif contains(ecogformatparts{i},{'[','x',']'}) % like C[4x8]
                formatparts(i) = str2double(extractBetween(ecogformatparts{i},'[','x')) * str2double(extractBetween(ecogformatparts{i},'x',']'));
            end
        end
        
        if sum(idx_gridpart) == sum(formatparts)
            idx_grid(idx_gridpart) = true;
        else            
                error('Error in "look_for_format.m", the number of electrodes in ecog-grid is not equal to the total of electrodes with the name %s',extractBefore(ecogformatall{n},'['))
        end
    end
    
else
    ecogformatfin = [];
end

% strip
idx_strip = false(size(metadata.ch) );
if ~isempty(striploc)
    stripformat = cell(1,size(striploc,2));
    for i=1:size(striploc,2)
        stripformat{i} = annotsformatsplit(striploc(i)+1: locs_all(find(locs_all==striploc(i))+1)-1);
    end
    stripformatall = strcat([stripformat{:}],';');
    stripformatfin = ['strip;' stripformatall{:}];
    
    % add in metadata which electrodes are in strip
    for n = 1:size(stripformatall,2)
        
        size_strippart = strlength(extractBefore(stripformatall{n},'[')) ;     
        % find which electrodes have same name, and same size (to avoid FML
        % to be part of ML, because both contains ML (RESP0991))
        idx_strippart1 = ismember(strSizeCh, size_strippart) == 1 & ...
            contains(metadata.ch,extractBefore(stripformatall{n},'[')) & ...
            metadata.ch2use_included == 1;
        
        % in RESP401, all electrodes are named sOc, except 4, which is
        % called sOC4... this lower part is added to make sure that the
        % number of electrodes is determined correctly in this patient as
        % well
        idx_strippart2 = ismember(strSizeCh, size_strippart) == 1 & ...
            contains(metadata.ch,extractBefore(stripformatall{n},'['),'IgnoreCase',true) & ...
            metadata.ch2use_included == 1;

        if sum(idx_strippart1) == str2double(extractBetween(stripformatall{n},'[','x')) * str2double(extractBetween(stripformatall{n},'x',']'))
            idx_strip(idx_strippart1) = true;
        elseif sum(idx_strippart2) == str2double(extractBetween(stripformatall{n},'[','x')) * str2double(extractBetween(stripformatall{n},'x',']'))
            idx_strip(idx_strippart2) = true;
        else
            error('Error in "look_for_format.m", the number of electrodes in strip is not equal to the total of electrodes with the name %s',extractBefore(stripformatall{n},'['))
        end
        
    end
else
    stripformatfin = [];
end

% depth
idx_depth = false(size(metadata.ch) );
if ~isempty(depthloc)
    depthformat = cell(1,size(depthloc,2));
    for i=1:size(depthloc,2)
        depthformat{i} = annotsformatsplit(depthloc(i)+1: locs_all(find(locs_all==depthloc(i))+1)-1);
    end
    depthformatall = strcat([depthformat{:}],';');
    depthformatfin = ['depth;' depthformatall{:}];
    
    % add in metadata which electrodes are in depth
    for n = 1:size(depthformatall,2)
        
        size_depthpart = strlength(extractBefore(depthformatall{n},'[')) ;     
        % find which electrodes have same name, and same size (to avoid FML
        % to be part of ML, because both contains ML (RESP0991))
        idx_depthpart = ismember(strSizeCh, size_depthpart) == 1 & ...
            contains(metadata.ch,extractBefore(depthformatall{n},'[')) & ...
            metadata.ch2use_included == 1;
                
        if sum(idx_depthpart) == str2double(extractBetween(depthformatall{n},'[','x')) * str2double(extractBetween(depthformatall{n},'x',']'))
            idx_depth(idx_depthpart) = true;
        else
            error('Error in "look_for_format.m", the number of electrodes in depth is not equal to the total of electrodes with the name %s',extractBefore(depthformatall{n},'['))
        end
    end
    
    
else
    depthformatfin = [];
end

% seeg
idx_seeg = false(size(metadata.ch) );
if ~isempty(seegloc)
    seegformat = cell(1,size(seegloc,2));
    for i=1:size(seegloc,2)
        seegformat{i} = annotsformatsplit(seegloc(i)+1: locs_all(find(locs_all==seegloc(i))+1)-1);
    end
    seegformatall = strcat([seegformat{:}],';');
    seegformatfin = ['seeg;' seegformatall{:}];
  
    % add in metadata which electrodes are in strip
    for n = 1:size(seegformatall,2)
                
        size_seegpart = strlength(extractBefore(seegformatall{n},'[')) ;     
        % find which electrodes have same name, and same size (to avoid FML
        % to be part of ML, because both contains ML (RESP0991))
        idx_seegpart = ismember(strSizeCh, size_seegpart) == 1 & ...
            contains(metadata.ch,extractBefore(seegformatall{n},'[')) & ...
            metadata.ch2use_included == 1;
        
        if sum(idx_seegpart) == str2double(extractBetween(seegformatall{n},'[','x')) * str2double(extractBetween(seegformatall{n},'x',']'))
            idx_seeg(idx_seegpart) = true;
        else
            error('Error in "look_for_format.m", the number of electrodes in seeg is not equal to the total of electrodes with the name %s',extractBefore(seegformatall{n},'['))
        end
    end

else
    seegformatfin = [];
end

metadata.format_info=[ecogformatfin, stripformatfin, depthformatfin, seegformatfin];
metadata.idx_seeg = idx_seeg;
metadata.idx_depth = idx_depth;
metadata.idx_grid = idx_grid;
metadata.idx_strip = idx_strip;
end
