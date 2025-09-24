function write_coordsystemJSON(cfg)

%%%% This script writes coordsystem JSON file
%%%% Make sure data matlab to JSON library is added
%%%% This can be found here: https://github.com/gllmflndn/JSONio

%%%% Dora Hermes, Jaap van der Aar, Giulio Castegnaro, Dorien van Blooijs 2019
%%

coord_json_name = fullfile(cfg(1).ieeg_directory,...
    [cfg(1).sub_label '_' cfg(1).ses_label '_coordsystem.json']);

% This line is to create a variable for the name of the file intendedFor
filename_T1w = fullfile(cfg(1).sub_label,cfg(1).ses_label, 'anat', [cfg(1).sub_label '_' cfg(1).ses_label '_rec-deface_T1w.nii']);

if exist(fullfile(cfg(1).proj_diroutput,filename_T1w),'file')
    % assign information and methodology
    loc_json.IntendedFor = filename_T1w;
    loc_json.iEEGCoordinateSystem  = 'Other';
    loc_json.iEEGCoordinateUnits  = 'mm';
    loc_json.iEEGCoordinateSystemDescription = 'The origin of the coordinate system is between the ears and the axis are in the RAS direction. The scaling is with respect to the individuals anatomical scan and no scaling or deformation have been applied to the individuals anatomical scan';
    loc_json.iEEGCoordinateProcessingDescription = 'Surface projection Hermes or Branco';
    loc_json.iEEGCoordinateProcessingReference = 'Hermes et al., 2010 JNeuroMeth , Branco et al., 2018 JNeuroMeth';

    jsonSaveDir = fileparts(coord_json_name);

    % ensure there is a /ieeg/ folder
    if ~isfolder(jsonSaveDir)
        fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
    end

    % % write JSON file
    if ~isempty(loc_json)

        delete(coord_json_name)
        write_json(coord_json_name, loc_json)
        % fileattrib(coord_json_name,'-w -x','o') % make not-writable and not-executable for other users
        % fileattrib(coord_json_name,'+w +x','g') % make writable and executable (required for folders to open them) for group users
        % % cannot be used in windows, to do make different system!

    end

end
end





