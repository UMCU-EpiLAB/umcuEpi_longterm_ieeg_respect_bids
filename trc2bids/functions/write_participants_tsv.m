function write_participants_tsv(cfg,header,metadata)

for i=1:size(cfg,2)
    
    if ~isempty(cfg(i).proj_diroutput)
        
        filename = fullfile(cfg(i).proj_diroutput,'participants.tsv');
        
        % find session number
        splitfile = strsplit(cfg(1).ieeg_dir{i},'\');
        sesnum = splitfile{contains(splitfile,'ses')};
        
        files = dir(cfg(i).proj_diroutput);
        pat_exist = [];
        if contains([files(:).name],'participants')
            % read existing scans-file
            participants_tsv = read_tsv(filename);
            
            if any(contains(participants_tsv.participant_id,deblank(header.name))) % look whether the name is already in the participants-table
                if ~isempty(find(contains(participants_tsv.participant_id,deblank(header.name)) & ...
                        contains(participants_tsv.session,sesnum)==1, 1)) % patient and session is already in participants-table
                    partnum = find(contains(participants_tsv.participant_id,deblank(header.name)) & ...
                        contains(participants_tsv.session,sesnum)==1, 1); %find patient number and session number
                    pat_exist = 1;
                elseif isempty(find(contains(participants_tsv.participant_id,deblank(header.name)) & ...
                        contains(participants_tsv.session,sesnum)==1, 1)) % this session is not yet in participants-table
                     partnum = size(participants_tsv,1)+1;
                end
            else % if participant is not yet in the table, the number is the last one plus one
                partnum = size(participants_tsv,1)+1;
            end
            
            participant_id = participants_tsv.participant_id;
            age = participants_tsv.age;
            session = participants_tsv.session;
            sex = participants_tsv.sex;
        else
            partnum = 1;
        end
        
        % set RESPect name and session number and sex
        participant_id{partnum,1}   = ['sub-' deblank(header.name)];
        session{partnum,1} = sesnum;
        
        if strcmpi(strip(metadata.gender),'male') || strcmpi(strip(metadata.gender),'female')
            sex{partnum,1} = strip(metadata.gender);
        elseif strcmp(strip(metadata.gender), 'unknown') && size(sex,1) < partnum
            sex{partnum,1} = strip(metadata.gender);
        end
        
        % set age of RESPect patient (comparing with current participants-table)
        if pat_exist == 1
            if any(contains(fieldnames(metadata),'age')) 
                age(partnum,1) = str2double(metadata.age);
                
            elseif age(partnum,1) == header.age && age(partnum,1) ~= 0 % if age in participants.tsv is not equal to 0  and equal to header.age
                age(partnum,1)    = header.age;
                
            elseif age(partnum,1) ~= 0 && header.age == 0 % if age is not equal to 0 (assumed to be correct)
                
            elseif age(partnum,1) == 0 && header.age ~= 0 % if age is equal to 0 and header.age is not (latter is assumed to be correct)
                age(partnum,1) = header.age;
                
            elseif age(partnum,1) ~= 0 && header.age ~= 0 && age(partnum,1) ~= header.age % if both ages are not 0 and conflicting, keep current age, because that might have been annotated in the trc file with Included, Format annotations
                warning('ages between this file and other file are in conflict!')
                
            elseif age(partnum,1) == 0 && header.age == 0
                warning('age is 0 years... assumed to be incorrect!')
            end
        else
            if any(contains(fieldnames(metadata),'age'))
                age(partnum,1) = str2double(metadata.age);
            elseif header.age == 0
                warning('age is 0 years... assumed to be incorrect! Annotate "Age;X" in a trc file!!')
                age(partnum,1) = header.age;
            else 
                age(partnum,1) = header.age;
            end
        end
        
        % extract RESPect numbers from RESPect names
        numname = zeros(size(participant_id));
        for n=1:size(participant_id,1)
            numname(n) = str2double(participant_id{n}(regexp(participant_id{n},'[0-9]')));
        end

        % convert session names to numbers (a adds 0.1, b adds 0.2)
        sesname = zeros(size(session));
        for n=1:size(session,1)
            sesname(n) = str2double(session{n}(regexp(session{n},'[0-9]')));
            sespart = session{n}(regexp(session{n},'[a-z]'));
            if strcmp(sespart,'a')
                sesname(n) = sesname(n) +0.1;
            elseif strcmp(sespart,'b')
                sesname(n) = sesname(n) + 0.2;
            end
        end
        
        % sorts table based on RESPect number and session number
        [~,I] = sortrows([numname,sesname]);
        
        participant_id_sort = participant_id(I);
        age_sort = age(I);
        sex_sort = sex(I);
        session_sort = session(I);
        
        % makes a table from name, session and age
        participants_tsv  = table(participant_id_sort, session_sort, age_sort, sex_sort, ...
            'VariableNames',{'participant_id','session', 'age', 'sex'});
        
        % save participants.tsv
        if ~isempty(participants_tsv)
            delete(filename)
            write_tsv(filename, participants_tsv);
            %fileattrib(filename,'-w -x','o') % make not-writable and not-executable for other users
            %fileattrib(filename,'+w +x','g') % make writable and executable (required for folders to open them) for group users
        end
    end
end