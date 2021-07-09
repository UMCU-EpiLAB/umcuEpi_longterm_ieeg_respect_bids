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
        idx_gridpart = ismember(strSizeCh, size_gridpart) == 1 & ...
            contains(metadata.ch,extractBefore(ecogformatall{n},'[')) == 1;
                
        if sum(idx_gridpart) == str2double(extractBetween(ecogformatall{n},'[','x')) * str2double(extractBetween(ecogformatall{n},'x',']'))
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
        idx_strippart = ismember(strSizeCh, size_strippart) == 1 & ...
            contains(metadata.ch,extractBefore(stripformatall{n},'[')) == 1;
                
        if sum(idx_strippart) == str2double(extractBetween(stripformatall{n},'[','x')) * str2double(extractBetween(stripformatall{n},'x',']'))
            idx_strip(idx_strippart) = true;
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
            contains(metadata.ch,extractBefore(depthformatall{n},'[')) == 1;
                
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
            contains(metadata.ch,extractBefore(seegformatall{n},'[')) == 1;
        
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
