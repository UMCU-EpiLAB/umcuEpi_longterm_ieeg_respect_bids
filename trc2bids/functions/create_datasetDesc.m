%% create dataset descriptor
function create_datasetDesc(proj_dir,sub_label)


if contains(sub_label,'RESP')
    ddesc_json.Name               = 'RESPect' ;
    ddesc_json.BIDSVersion        = 'Brain Imaging Data Structure Specification v1.6.0';
    ddesc_json.License            = 'This dataset is made available under the Public Domain Dedication and License CC v1.0, whose full text can be found at https://creativecommons.org/publicdomain/zero/1.0/. We hope that all users will follow the ODC Attribution/Share-Alike Community Norms (http://www.opendatacommons.org/norms/odc-by-sa/); in particular, while not legally required, we hope that all users of the data will acknowledge by citing Demuru M, van Blooijs D, Zweiphenning W, Hermes D, Leijten F, Zijlmans M, on behalf of the RESPect group. “A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS”. Submitted to NeuroInformatics in 2020 in any publications.';
    ddesc_json.Authors            = {'van Blooijs D.', 'Demuru M.', 'Zweiphenning W', 'Leijten F', 'Zijlmans M.'};
    ddesc_json.Acknowledgements   = 'Huiskamp G.J.M.';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2021, A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS' ;
    ddesc_json.Funding            = {'Epi-Sign Project', 'Alexandre Suerman Stipendium 2015', 'EpilepsieNL #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'see HowToAcknowledge'};
    ddesc_json.DatasetDOI         = '10.18112/openneuro.ds003399.v1.0.2';
    
elseif contains(sub_label,'PRIOS')
    ddesc_json.Name               = 'PRIOS' ;
    ddesc_json.BIDSVersion        = 'Brain Imaging Data Structure Specification v1.6.0';
    ddesc_json.License            = 'This dataset is made available under the Public Domain Dedication and License CC v1.0, whose full text can be found at https://creativecommons.org/publicdomain/zero/1.0/. We hope that all users will follow the ODC Attribution/Share-Alike Community Norms (http://www.opendatacommons.org/norms/odc-by-sa/); in particular, while not legally required, we hope that all users of the data will acknowledge by citing Demuru M, van Blooijs D, Zweiphenning W, Hermes D, Leijten F, Zijlmans M, on behalf of the RESPect group. “A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS”. Submitted to NeuroInformatics in 2020 in any publications.';
    ddesc_json.Authors            = {'Blok S.', 'van Blooijs D.', 'Huiskamp G.J.M.', 'Leijten F.S.S.'};
    ddesc_json.Acknowledgements   = 'persons to acknowledge';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2021, A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS' ;
    ddesc_json.Funding            = {'EpilepsieNL #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'see HowToAcknowledge'};
    ddesc_json.DatasetDOI         = 'DOI of the dataset if online';
    
elseif contains(sub_label,'REC2Stim')
    ddesc_json.Name               = 'REC2Stim' ;
    ddesc_json.BIDSVersion        = 'Brain Imaging Data Structure Specification v1.6.0';
    ddesc_json.License            = 'This dataset is made available under the Public Domain Dedication and License CC v1.0, whose full text can be found at https://creativecommons.org/publicdomain/zero/1.0/. We hope that all users will follow the ODC Attribution/Share-Alike Community Norms (http://www.opendatacommons.org/norms/odc-by-sa/); in particular, while not legally required, we hope that all users of the data will acknowledge by citing Demuru M, van Blooijs D, Zweiphenning W, Hermes D, Leijten F, Zijlmans M, on behalf of the RESPect group. “A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS”. Submitted to NeuroInformatics in 2020 in any publications.';
    ddesc_json.Authors            = {'van Blooijs D.', 'Aarnoutse E.J.', 'Ramsey N.F.', 'Huiskamp G.J.M.', 'Leijten F.S.S.'};
    ddesc_json.Acknowledgements   = 'persons to acknowledge';
    ddesc_json.HowToAcknowledge   = 'Demuru M., van Blooijs D., Zweiphenning W. et al 2021, A practical workflow for organizing clinical intraoperative and long-term iEEG data in BIDS' ;
    ddesc_json.Funding            = {'EpilepsieNL #17-07'} ;
    ddesc_json.ReferencesAndLinks = {'see HowToAcknowledge'};
    ddesc_json.DatasetDOI         = 'DOI of the dataset if online';
    
else
    error('patient is not correctly anonymized');
end

if ~isempty(ddesc_json)
    
    filename = fullfile(proj_dir,'dataset_description.json');
    delete(filename)
    write_json(filename, ddesc_json)
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ddesc_json), json_options)
    fileattrib(filename,'-w -x','o') % make not-writable and not-executable for other users
    fileattrib(filename,'+w +x','g') % make writable and executable (required for folders to open them) for group users

end
