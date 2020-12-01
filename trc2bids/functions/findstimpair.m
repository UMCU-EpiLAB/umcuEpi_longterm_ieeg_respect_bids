function [stimchan,stimnum] = findstimpair(stimnums,stimchans,metadata)

ch_label = metadata.ch_label;

stimchan = cell(size(stimnums));
stimnum = zeros(size(stimnums));
for j=1:size(stimnums,2)
    test1 = sprintf('%s%d',stimchans{j},str2double(stimnums{j}));
    test2 = sprintf('%s0%d',stimchans{j},str2double(stimnums{j}));
    
    % check if one of the tests is present in ch_label, and that this
    % ch_label is one of the ecog channels
    if sum(strcmpi(test1,ch_label))>0 && ismember(find(strcmpi(test1,ch_label)),find(metadata.ch2use_included)) 
        if sum(strcmpi(test1,ch_label))==1
            stimchan{j} = ch_label{strcmpi(test1,ch_label)};
            stimnum(j) = find(strcmpi(test1,ch_label)==1);
        else
            stimchan{j} = ch_label{strcmp(test1,ch_label)};
            stimnum(j) = find(strcmp(test1,ch_label)==1);
        end
    elseif sum(strcmpi(test2,ch_label))>0 && ismember(find(strcmpi(test2,ch_label)),find(metadata.ch2use_included)) 
        if sum(strcmpi(test2,ch_label))==1
            stimchan{j} = ch_label{strcmpi(test2,ch_label)};
            stimnum(j) = find(strcmpi(test2,ch_label)==1);
        else
            stimchan{j} = ch_label{strcmp(test2,ch_label)};
            stimnum(j) = find(strcmp(test2,ch_label)==1);
        end
        
    end
end
