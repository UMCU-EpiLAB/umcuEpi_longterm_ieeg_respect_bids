%% create dataset descriptor
function create_participantsDesc(proj_dir)

partdesc_json.participant_id.Description    = 'name of the patient' ;
partdesc_json.session.Description           = 'session number, session 1 is the first period of invasive recordings, 1a/b means that the electrodes have changes during the first session';
partdesc_json.age.Description               = 'age of the patient during the session';
partdesc_json.age.Units                     = 'years';
partdesc_json.sex.Description               = 'sex of the patient (male/female/unknown)';
partdesc_json.sex.Levels.male               = 'male';
partdesc_json.sex.Levels.female             = 'female';
partdesc_json.sex.Levels.unknown            = 'unknown';

if ~isempty(partdesc_json)
    
    filename = fullfile(proj_dir,'participants.json');
    delete(filename)
    write_json(filename, partdesc_json)
    fileattrib(filename,'-w -x','o') % make not-writable and not-executable for other users
    fileattrib(filename,'+w +x','g') % make writable and executable (required for folders to open them) for group users

end
