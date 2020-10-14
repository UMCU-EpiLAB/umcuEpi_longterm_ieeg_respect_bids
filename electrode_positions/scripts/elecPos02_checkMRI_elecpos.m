%% check electrode positions and move to /Fridge/CCEP and /Fridge/chronic_ECoG
% author: Dorien van Blooijs
% date: September 2019

%% patient characteristics - matlab

clear
cfg(1).sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};
cfg(1).no_fieldtrip = 'yes';
cfg(1).mode = 'electrodeposition_preMRI';

% set paths
cfg = setLocalDataPath(cfg);

tb_elecs = readtable(fullfile(cfg(1).ieeg_directory,...
    [cfg(1).sub_labels{:} '_' cfg(1).ses_label, '_electrodes.tsv']),'FileType','text','Delimiter','\t');


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

check_atlas_elec_MRI(cfg(1),tb_elecs)

