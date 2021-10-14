%% create dataset descriptor
function create_eventDesc(proj_dir,sub_label,ses_label)

edesc_json.onset.type                                       = 'object';
edesc_json.onset.description                                = 'onset of event in seconds' ;
edesc_json.duration.type                                    = 'object';
edesc_json.duration.description                             = 'duration of event in seconds' ;
edesc_json.trial_type.type                                  = 'object';
edesc_json.trial_type.description                           = 'type of event (electrical stimulation/motor task/sensing task/artefact/sleep/sleep wake transition/eyes open)' ;
edesc_json.sub_type.type                                    = 'object';
edesc_json.sub_type.description                             = 'more description of event (sleep:nrem/rem, motor:Mario/hand/jump, sens:circle, electrical stimulation:SPES/ESM/REC2stim, seizure:clinical/subclinical)' ;
edesc_json.electrodes_involved_onset.type                   = 'object';
edesc_json.electrodes_involved_onset.description            = 'electrodes involved in onset. For example: electrodes involved in seizure onset or in artefact.' ;
edesc_json.electrodes_involved_offset.type                  = 'object';
edesc_json.electrodes_involved_offset.description           = 'electrodes involved in offset. For example: electrodes involved in the end of a seizure or in an artefact.' ;
edesc_json.offset.type                                      = 'object';
edesc_json.offset.description                               = 'offset of event in seconds' ;
edesc_json.sample_start.type                                = 'object';
edesc_json.sample_start.description                         = 'onset of event in samples' ;
edesc_json.sample_end.type                                  = 'object';
edesc_json.sample_end.description                           = 'offset of event in samples' ;
edesc_json.electrical_stimulation_type.type                 = 'object';
edesc_json.electrical_stimulation_type.description          = 'type of electrical stimulation [mono-/biphasic]';
edesc_json.electrical_stimulation_site.type                 = 'object';
edesc_json.electrical_stimulation_site.description          = 'electrode names of stimulus pair';
edesc_json.electrical_stimulation_current.type              = 'object';
edesc_json.electrical_stimulation_current.description       = 'electrical stimulation current in Ampere';
edesc_json.electrical_stimulation_frequency.type            = 'object';
edesc_json.electrical_stimulation_frequency.description     = 'electrical stimulation frequency in Hertz';
edesc_json.electrical_stimulation_pulsewidth.type           = 'object';
edesc_json.electrical_stimulation_pulsewidth.description    = 'electrical stimulation pulse width in s';
edesc_json.notes.type                                       = 'object';
edesc_json.notes.description                                = 'notes about the specific event';

if ~isempty(edesc_json)
    
    filename = fullfile(proj_dir,sub_label,ses_label,'ieeg',[sub_label,'_',ses_label, '_events.json']);
    delete(filename)
    write_json(filename, edesc_json)
    
    fileattrib(filename,'-w -x','o') % make not-writable and not-executable for other users
    fileattrib(filename,'+w +x','g') % make writable and executable (required for folders to open them) for group users

end
