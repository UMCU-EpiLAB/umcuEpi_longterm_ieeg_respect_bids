function [annots_new, eventsannots ] = add_event2annotation(event,evname,eventsannots, annots_new, header)

fs = header.Rate_Min;

if size(eventsannots.type,1) == 0
    cc = 1;
else
    cc = size(eventsannots.type,1)+1;
end

if(~isempty(event))
    for i=1:numel(event)
        
        eventsannots.type{cc,1}    = evname                           ;
        
        %         if strcmp(evname,'artefact')
        eventsannots.sub_type{cc,1} = event{i}.type;
        eventsannots.s_start{cc,1} = event{i}.pos(1)/fs; % time in seconds 
        eventsannots.samp_start{cc,1} = num2str(event{i}.pos(1))          ;
        eventsannots.s_end{cc,1}   = event{i}.pos(end)/fs; % time in seconds 
        eventsannots.samp_end{cc,1} = num2str(event{i}.pos(end))          ;
        eventsannots.duration{cc,1} = eventsannots.s_end{cc} - eventsannots.s_start{cc};
        eventsannots.pulsewidth{cc,1} = 'n/a';
        eventsannots.stim_type{cc,1} = 'n/a';
        eventsannots.site_name{cc,1} = 'n/a';
        eventsannots.stim_cur{cc,1} = 'n/a';
        eventsannots.freq{cc,1} = 'n/a';
        eventsannots.pulses{cc,1} = 'n/a';
        eventsannots.isi{cc,1} = 'n/a';
       
        if isfield(event{i},'notes') && ~isempty(event{i}.notes)
            eventsannots.notes{cc,1} = event{i}.notes;
        else
            eventsannots.notes{cc,1} = 'n/a';
        end
        
        clear name 
        if size(event{i}.ch_names_on,2) == 1
            name = event{i}.ch_names_on{1}              ;
        elseif size(event{i}.ch_names_on,2) > 1
            for j=1
                name =[event{i}.ch_names_on{1}];
            end
            
            for j=2:size(event{i}.ch_names_on,2)
                name = [name ,',', event{i}.ch_names_on{j}]; %#ok<AGROW>
            end
        else % if size < 1, this is wrong!
            error('Number of channels in seizure onset is annotated incorrectly!')
        end
        eventsannots.ch_name_on{cc,1} = name;
        
        clear name
        if size(event{i}.ch_names_off,2) == 1
            name = event{i}.ch_names_off{1}              ;
        elseif size(event{i}.ch_names_off,2) > 1
            for j=1
                name =[event{i}.ch_names_off{1}];
            end
            
            for j=2:size(event{i}.ch_names_off,2)
                name = [name ,',', event{i}.ch_names_off{j}]; %#ok<AGROW>
            end
        else % if size <1, this is wrong!
            error('Number of channels in seizure offset is annotated incorrectly!')
        end
        
        eventsannots.ch_name_off{cc,1} = name;
        
        annots_new([annots_new{:,1}]==event{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==event{i}.pos(end),:)=[];
        
        cc          = cc + 1                               ;
        
        %         end
    end
end
