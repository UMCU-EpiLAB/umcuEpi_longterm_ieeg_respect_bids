%% check electrode positions 
% of two sessions
% date: November 2020

%% patient characteristics - matlab

clear
cfg(1).sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};
cfg(1).no_fieldtrip = 'yes';
cfg(1).mode = 'electrodeposition_preMRI';

% set paths
cfg = setLocalDataPath(cfg);

if strcmpi(cfg(1).ses_label,'ses-1')
   rep_ses = 'ses-2';
elseif strcmpi(cfg(1).ses_label,'ses-2')
    rep_ses = 'ses-1';
end

tb_elecs_ses_1 = readtable(fullfile(cfg(1).ieeg_directory,...
    [cfg(1).sub_labels{:} '_' cfg(1).ses_label, '_electrodes.tsv']),'FileType','text','Delimiter','\t');

tb_elecs_ses_2 = readtable(fullfile(replace(cfg(1).ieeg_directory,cfg(1).ses_label,rep_ses),...
    [cfg(1).sub_labels{:} '_', rep_ses, '_electrodes.tsv']),'FileType','text','Delimiter','\t');

cfg(1).deriv_directory = replace(cfg(1).deriv_directory,[cfg(1).ses_label ,'/'],'');

%%
close all

cfg(1).show_labels = 'yes';
cfg(1).change_size = 'no';
cfg(1).change_color = 'no';
cfg(1).view_atlas ='no';
cfg(1).elec_offset = 'yes';
% cfg(1).atlas = 'Destrieux'; % [DKT/Destrieux]
cfg(1).atlas = 'DKT'; % [DKT/Destrieux]
cfg(1).view_elec ='yes';
cfg(1).save_fig = 'yes'; % saves rendering in derivatives

check_atlas_elec_MRI_twoses(cfg(1),tb_elecs_ses_1,tb_elecs_ses_2)

