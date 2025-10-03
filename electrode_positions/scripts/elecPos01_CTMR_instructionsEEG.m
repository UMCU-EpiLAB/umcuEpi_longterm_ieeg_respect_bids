% code to localize subdural eeg and stereo eeg electrodes
%   Created by:
%
%     Copyright (C) 2009  D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%                   2019  D. van Blooijs, Dept of Neurology and Neurosurgery, University Medical Center Utrecht

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% DvB - made it compatible with BIDS electrodes.tsv September 2019
% SJ & IH - made it compatible for windows 2025 

% This script has a part that should be run in a linux terminal, and part
% that can be run in matlab. The parts that should be run in a linux
% terminal have "run in linux terminal" in the section title.
% see also the corresponding IFU, 5_IFU_electrode positions_procespc

% TO DO:
% - STEP 13 solve that elec pos can be plotted onto mri
% - STEP 9 and 10 change with script nicole
% - STEP 4 for a lot of pt
% - STEP 4/5 copy all files with symbolic link and not only the pial
%% patient characteristics - matlab

clear
close all
cfg(1).sub_label = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};
cfg(1).sub_labels = cfg(1).sub_label; % needed for some trc2bids functions
cfg(1).no_fieldtrip = 'yes';
cfg(1).mode = 'electrodeposition_preMRI';


% set paths
cfg = setLocalDataPath(cfg);
%%
rootPath = matlab.desktop.editor.getActiveFilename;
splitPath = regexp(rootPath,'\','split');
indexSpmPath = find(contains(splitPath,'umcuEpi_longterm_ieeg_respect_bids'));
spmPath = fullfile(splitPath{1:indexSpmPath-1},'spm12');
addpath(genpath(spmPath));
freesurferPath = fullfile(splitPath{1:indexSpmPath-1},'freesurfer');
addpath(genpath(freesurferPath));

% householding
clear rootPath splitPath spmPath indexSpmPath freesurferPath
%% STEP 3: defacing MRI - RUN IN LINUX TERMINAL- check in windows
% CHECK if this is not already done manually in sourcedata, than skip this
% step
% run this part in matlab with 'ctrl enter', this will show the text to copy in the terminal 
clc

% Find the original MRI and copy this to the folder you're using.
% Rename the original MRI to (run the line below and copy the printed line in the command window):
fprintf('\n ----- RENAME T1WEIGHTED MRI TO: ----- \n %s_%s_T1w.nii\n',...
    cfg(1).sub_label,...
    cfg(1).ses_label);

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window to deface the MRI in the linux terminal:


fprintf('\n ----- OPEN LINUX TERMINAL, RUN LINE BELOW IN LINUX TERMINAL -----\n cd %ssourcedata/%s/%s/anat/ \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n su  \n password: Dr.House \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n mri_deface %s_%s_T1w.nii %s  %s %s_%s_rec-deface_T1w.nii\n',...
    cfg(1).proj_diroutput_linux,...
    cfg(1).sub_label,...
    cfg(1).ses_label,...
    cfg(1).sub_label,...
    cfg(1).ses_label,...
    cfg(1).path_talairach,...
    cfg(1).path_face,...
    cfg(1).sub_label,...
    cfg(1).ses_label);
% this takes around 5 minutes

fprintf('\n ----- OPEN MRICRON in windows, OPEN DEFACED MRI TO CHECK DEFACING ----- \n mricron \n')

%% STEP 4 : run freesurfer to segment brain add Destrieux atlases - RUN IN LINUX TERMINAL
% run this part in matlab with 'ctrl enter', this will show the text to copy in the terminal 
clc

mri_name = [cfg(1).anat_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii'];

% Copy defaced .nii to correct folder

if exist(cfg(1).anat_directory,'dir')
    copyfile(fullfile(cfg(1).proj_diroutput,'sourcedata',cfg(1).sub_label,cfg(1).ses_label,'anat',...
        [cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii']),...
        mri_name)
    
else
    
    mkdir(cfg(1).anat_directory)
    copyfile(fullfile(cfg(1).proj_diroutput,'sourcedata',cfg(1).sub_label,cfg(1).ses_label,'anat',...
        [cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii']),...
        mri_name)
    
end

% Make a freesurfer folder
if exist(cfg(1).freesurfer_directory, 'dir')
    %     fprintf('\n%s exists already\n',cfg(1).freesurfer_directory)
else
    mkdir(cfg(1).freesurfer_directory)
end

% log into Oracle virtual box linux environment
fprintf(' \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n su  \n password: Dr.House \n')

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n cp %s%s_%s_rec-deface_T1w.nii %s\n',...
    cfg(1).anat_directory_linux,...
    cfg(1).sub_label,...
    cfg(1).ses_label, ...
    '/home/epilab/Desktop/')

%Copy the printed lines in the command window into the linux terminal:
fprintf('\n -----  RUN LINE BELOW IN LINUX TERMINAL ----- \n export SUBJECTS_DIR=%s \n', ...
     '/home/epilab/Desktop/')

% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n recon-all -autorecon-all -s %s -i %s%s_%s_rec-deface_T1w.nii -cw256\n',...
    cfg(1).sub_label,...
    '/home/epilab/Desktop/',...
    cfg(1).sub_label,...
    cfg(1).ses_label)

% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
% -P is a test 23/6 for symbolic links, see if it works
fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n cp -r -P %s%s/  %s\n',...
    '/home/epilab/Desktop/',...
    cfg(1).sub_label,...
    cfg(1).freesurfer_directory_linux)

% This takes up to 12 hours to run! In the end, you will see a subject
% folder in the freesurfer folder.
%% STEP 4 copy works but not symbolic link (see step 13 for the symbolic link part) and manually do RESP*** folder one up
%% STEP 5: convert freesurfer file to .gii - RUN IN LINUX TERMINAL
% WATCH OUT: Do for left and right hemisphere if in both hemispheres (you
% see two batches of code in the command window)
% electrodes
% run this part in matlab with 'ctrl enter', this will show the text to copy in the terminal 
clc

% Make surface folder
if exist(cfg(1).surface_directory, 'dir')
    %     fprintf('%s exists already\n',cfg(1).surface_directory)
else
    mkdir(cfg(1).surface_directory)
end

% start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
for i=1:size(cfg(1).hemisphere,2)
    fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n cp -L /home/epilab/Desktop/%s/surf/%sh.pial %ssurf \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n cd %ssurf \n ----- RUN LINE BELOW IN LINUX TERMINAL -----  \n mris_convert %sh.pial %sh.pial.gii\n',cfg(1).sub_label,cfg(1).hemisphere{i},cfg(1).freesurfer_directory_linux,cfg(1).freesurfer_directory_linux,cfg(1).hemisphere{i},cfg(1).hemisphere{i})
end


%% STEP 6: generate surface (The Hull) to project electrodes to - RUN IN LINUX TERMINAL
% run this part in matlab with 'ctrl enter', this will show the text to copy in the terminal 
clc
% only for ECoG, because this is necessary to correct for brain-shift.
% ECoG electrodes are projected to the hull.

% Freesurfer creates a file called /mri/ribbon.mgz and we want to convert
% this file to a .nii file so we can read it and use it to create the hull
% (a tight balloon of where the electrodes should be on the pre-op MRI)

% Right click in the freesurfer/mri-folder and start Linux terminal.

% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %smri ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nmri_convert ribbon.mgz t1_class.nii\n',cfg(1).freesurfer_directory_linux)


%% STEP 7: Create the hull - matlab
if ~exist(cfg(1).deriv_directory,'dir')
    mkdir(cfg(1).deriv_directory)
end

for i=1:size(cfg(1).hemisphere,2)
    settings_hull = [13,... % setting for smoothing: default 13
        0.3]; % setting for threshold: default 0.3
    k = get_mask_V3(cfg(1).sub_label,... % subject name
        [cfg(1).freesurfer_directory,'mri\t1_class.nii'],... % freesurfer class file %
        cfg(1).deriv_directory,... % where you want to safe the file
        cfg(1).hemisphere{i},... % 'l' for left 'r' for right --> only use the hemisphere where electrodes are located
        settings_hull(1),...% setting for smoothing
        settings_hull(2)); % settings for  threshold
    % the hull is saved as sub-RESPXXXX_surface1_13_03.img
end

%% STEP 8: check hull - mricron
% run this part in matlab with 'ctrl enter', this will show the text to copy in the terminal 
% type 'mricron'
% load the MRI
% put the hull as overlay on top of the mri
% check whether the hull looks like it matches the dura (should be a tight
% baloon around the grey matter)

fprintf('\n ----- OPEN MRICRON in windows, OPEN DEFACED MRI AND PUT HULL AS OVERLAY ON TOP ----- \n mricron \n \n ----- CHECK WHETHER THE HULL IS A TIGHT BALLOON AROUND THE CORTEX ----- \n')

%% STEP 9: segment electrodes from ct - nicole's function

%% STEP 9A: LOAD DATA
clc;
clear all;

addpath('\Matlab\spm12_updatesr7487')
addpath('\Matlab\fieldtrip-20201103')
addpath('\Afwijkingen_berekenen\')
ft_defaults

%load files
[ctname, ctpathname] = uigetfile('*.nii', 'Give the .nii file of the subject''s CT file'); % see werkmap
ct_elek_to_gd= [ctpathname  ctname];

[t1name, t1pathname] = uigetfile('*.nii', 'Give the .nii file of the subject''s T1 MRI file'); % see werkmap
t1_to_gd= [t1pathname  t1name];

% Clear temporary variables
clear ctname ctpathname mriname t1name  

% read coordinates from XML files
%assumes that XML files are in the same folder as MRI gado
xmlpathname = t1pathname; %[mripathname 'elektroden_brainlab\']; when xml files are in another folder
list = dir([xmlpathname '*.xmlSE']); % mripathname when xml files are in the same folder
[tmp, ind]=sort_nat({list.name});
list=list(ind);

for i=1:length(list)
    %load file
    filename = [xmlpathname list(i).name]; %mripathname when xml files are in the same folder
    delimiter = ' ';
    formatSpec = '%q%q%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');

    % Read columns of data according to format string.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

    % Close the text file.
    fclose(fileID);

    % Allocate imported array to column variable names
    electrode = dataArray{:, 1};
    version170403 = dataArray{:, 2};

    % Clear temporary variables
    clearvars filename delimiter formatSpec fileID dataArray ans;

    % CHECK REGELNUMMER WAAR specifiedTargetPoint en specifiedEntryPoint
    % staat in 'electrode'
    j=i*2-1;
    k=j+1;
    points(j,1)=str2double(electrode{219}(9:end-9)); %target
    points(j,2)=str2double(electrode{220}(9:end-9));
    points(j,3)=str2double(electrode{221}(9:end-9));
    points(k,1)=str2double(electrode{224}(9:end-9)); %entry
    points(k,2)=str2double(electrode{225}(9:end-9));
    points(k,3)=str2double(electrode{226}(9:end-9));

end

sprintf('aantal elektroden %d',i)

% Clear temporary variables
clear i list tmp xmlpathname ans electrode ind j k version170403

% x and y coordinate *-1 --> now sure why, but this is necessary
points(:,1)=points(:,1).*-1;
points(:,2)=points(:,2).*-1;

%% STEP 9B: get electrode names and length as input from user
prompt = [];
for n=1:size(points,1)/2
    prompt= [prompt {['Name electrode ' num2str(n) ':']}  {['Size electrode ' num2str(n) ':']}];
end
name='Specify electrode properties';

input1=[];input2=[];input3=[];
if numel(prompt)>12 %split if too large
    prompt1 = prompt(1:12);
    input1 = inputdlg(prompt1,name,[1 35]);
    prompt(1:12)=[];
   if numel(prompt)>12 
    prompt2 = prompt(1:12);
    input2 = inputdlg(prompt2,name,[1 35]);
    prompt(1:12)=[];
   end
   %handle remainder
   input3=inputdlg(prompt,name,[1 35]);
else
    input3=inputdlg(prompt,name,[1 35]);
end

input = [input1; input2; input3];

%create electrode table
electrode_table=cell(1,2);
for n=1:2:numel(input)
    electrode_table = [electrode_table; input(n) str2num(input{n+1})];
end
%remove empty row
electrode_table(1,:)=[];

clear input input1 input2 input3 prompt prompt1 prompt2



%% STEP 9C: work with CT
ct = ft_read_mri(ct_elek_to_gd,'dataformat','nifti_spm'); % source/ct met elektroden
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.viewresult  = 'yes';
cfg.coordsys    =  'dicom'; 

voxinds = round(ft_warp_apply(pinv(ct.transform), points)); 

% threshold CT data
threshold = 2300;
ct_coreg_threshold = ct;
ct_coreg_threshold.anatomy = ct.anatomy>threshold; % SJ ct = ct_coreg

% segement brain
mriT1 = ft_read_mri(t1_to_gd,'dataformat','nifti_spm');  % SJ t1_to_gd = t1
mriT1.coordsys    =  'dicom'; 
cfg=[];
cfg.output='tpm';
mriT1.coordsys='spm';
tpm=ft_volumesegment(cfg,mriT1);

inside = tpm.gray | tpm.white | tpm.csf;
ct_coreg_threshold.anatomy = ct_coreg_threshold.anatomy & inside;

% Clear temporary variables
clear threshold 

% combine T1 and CT
combined = mriT1;
combined.anatomy = mriT1.anatomy + double(ct_coreg_threshold.anatomy)*1000;
%plot
cfg         = [];
cfg.method = 'ortho';
cfg.interactive = 'no';
cfg.funparameter = 'mask';
ft_sourceplot(cfg, combined);

%% STEP: MAKE IMAGES: plot plane through points

% Identify electrodes
% search for closest voxels to planned electrodes -- not the tip of the
% electrode!
[ind_ct(:,1),ind_ct(:,2),ind_ct(:,3)] = ind2sub(size(ct_coreg_threshold.anatomy),find(ct_coreg_threshold.anatomy)); %indexes of ct positive voxels
closest = dsearchn(ind_ct,voxinds);
closest3 = ind_ct(closest,:);
% 
% %make pictures of each electrode
% for i= 1:2:length(closest3)
%     t = [closest3(i,2) closest3(i,1) closest3(i,3)];
%     e = [closest3(i+1,2) closest3(i+1,1) closest3(i+1,3)];
% 
%     close all;
%     figure;
% 
%     extra1 = [t(1) t(2) t(3)+50];   
%     extra2 = [t(1) t(2)+150 t(3)]; 
%     extra3 = [t(1)+150 t(2) t(3)]; 
% 
%     % plot1
%     % get two vectors
%     v1 = [t(1)-e(1) t(2)-e(2) t(3)-e(3)];
%     v2 = [t(1)-extra1(1) t(2)-extra1(2) t(3)-extra1(3)];
% 
%     %get normal to plane; cross prdocut of vectors
%     n = cross(v1,v2);
% 
%     %calculate slice
%     [SL1,x1,y1,z1] = obliqueslice(double(combined.anatomy),t,n);
%     xdist1=abs(x1(1,1)-x1(1,end))+abs(x1(1,1)-x1(end,1));
%     ydist1=abs(y1(1,1)-y1(1,end))+abs(y1(1,1)-y1(end,1));
%     zdist1=(abs(z1(1,1)-z1(1,end))+abs(z1(1,1)-z1(end,1)))*2; %z is maar tot 180, y en x tot 480
% 
%     % plot mri in plane
%     subaxis(2, 2, 1, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%     if xdist1<ydist1 && xdist1<zdist1 %axial plane
%         if z1(1,1)<z1(end,end) 
%             imshow(fliplr(SL1),[])
%         else 
%             imshow(SL1,[]);
%         end  
%     elseif ydist1<xdist1 && ydist1<zdist1 %coronal plane
%         if z1(1,1)>z1(end,end) 
%             imshow(fliplr(SL1),[])
%             camroll(90)
%         else 
%             imshow(SL1,[]);
%             camroll(-90);
%         end     
%     else %sagittal
%         imshow(SL1,[]);
%         if y1(1,1)>y1(end,end)
%             camroll(-90);
%         else
%             camroll(90);
%         end
%     end  
% 
%     %plot2
%     v2 = [t(1)-extra2(1) t(2)-extra2(2) t(3)-extra2(3)];
%     n = cross(v1,v2);
%     [SL2,x2,y2,z2] = obliqueslice(double(combined.anatomy),t,n);
%     xdist2=abs(x2(1,1)-x2(1,end))+abs(x2(1,1)-x2(end,1));
%     ydist2=abs(y2(1,1)-y2(1,end))+abs(y2(1,1)-y2(end,1));
%     zdist2=(abs(z2(1,1)-z2(1,end))+abs(z2(1,1)-z2(end,1)))*2;
% 
%     subaxis(2, 2, 2, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%     if xdist2<ydist2 && xdist2<zdist2 %axial plane
%         if z2(1,1)<z2(end,end) 
%             imshow(fliplr(SL2),[])
%         else 
%             imshow(SL2,[]);
%         end     
%     elseif ydist2<xdist2 && ydist2<zdist2 %coronal plane
%         if z2(1,1)>z2(end,end) 
%             imshow(fliplr(SL2),[])
%             camroll(90)
%         else 
%             imshow(SL2,[]);
%             camroll(-90);
%         end     
%     else %sagittal
%         imshow(SL2,[]);
%         if y2(1,1)>y2(end,end)
%             camroll(-90);
%         else
%             camroll(90);
%         end
%     end  
% 
% 
%     %plot 3
%     v2 = [t(1)-extra3(1) t(2)-extra3(2) t(3)-extra3(3)];
%     n = cross(v1,v2);
%     [SL3,x3,y3,z3] = obliqueslice(double(combined.anatomy),t,n);
%     xdist3=abs(x3(1,1)-x3(1,end))+abs(x3(1,1)-x3(end,1));
%     ydist3=abs(y3(1,1)-y3(1,end))+abs(y3(1,1)-y3(end,1));
%     zdist3=(abs(z3(1,1)-z3(1,end))+abs(z3(1,1)-z3(end,1)))*2;
% 
%     subaxis(2, 2, 3, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%     if xdist3<ydist3 && xdist3<zdist3 %axial plane
%         if z3(1,1)<z3(end,end)  
%             imshow(fliplr(SL3),[])
%         else 
%             imshow(SL3,[]);
%         end  
%     elseif ydist3<xdist3 && ydist3<zdist3 %coronal plane
%         if z3(1,1)>z3(end,end)  
%             imshow(fliplr(SL3),[])
%             camroll(90)
%         else 
%             imshow(SL3,[]);
%             camroll(-90);
%         end     
%     else %sagittal
%         imshow(SL3,[]);
%         if y3(1,1)>y3(end,end)
%             camroll(-90);
%         else
%             camroll(90);
%         end
%     end    
% 
%     subaxis(2, 2, 4, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%     imshow(zeros(50,50));
%     title(['Electrode ' num2str(round(i/2)) ' ' electrode_table{round(i/2),1}]);
%     set(gcf,'position',[50,50,1000,1200])
%     saveas(gcf,[t1pathname 'Electrode ' num2str(round(i/2)) ' .jpg'])
% 
% end

clear x1 x2 x3 y1 y2 y3 extra1 extra2 extra3 ydist1 ydist2 ydist3 xdist1 xdist2 xdist3 zdist1 zdist2 zdist3

%% STEP 9D: SEGMENT ALL ELECTRODES


% get voxels between entry and target 
%calculate voxel clusters
ConnectedComponents = bwconncomp(ct_coreg_threshold.anatomy);

figure;
linearInd=[];
for i=1:2:length(closest3)
    t = [closest3(i,2) closest3(i,1) closest3(i,3)];
    e = [closest3(i+1,2) closest3(i+1,1) closest3(i+1,3)];
    [x,y,z,~]= improfile3D(ct_coreg_threshold.anatomy,t,e); %get all voxels between target and entry voxel
    line3D = [x y z];
    scatter3(line3D(:,2),line3D(:,1),line3D(:,3),'.','g'); %x en y coordinaat moeten worden omgedraaid??
    hold on; 
    linearInd{i}= sub2ind(size(ct_coreg_threshold.anatomy),y,x,z); %switch x and y? get linear index in volume
end

for p=1:numel(ConnectedComponents.PixelIdxList)
    [x1,y1,z1]=ind2sub(size(ct_coreg_threshold.anatomy),ConnectedComponents.PixelIdxList{p});
    scatter3(x1,y1,z1,'.','r');
end

linearInd(cellfun(@isempty,linearInd))=[];

%check which voxel clusters are on each electrode line
oncomp=[];
for o=1:numel(ConnectedComponents.PixelIdxList) %for every connectedcomponent
    for i=1:length(linearInd) %for every electrode
         for j=1:numel(linearInd{i}) %for every point on line
            if any(ConnectedComponents.PixelIdxList{o}==linearInd{i}(j))
                oncomp(i,j)=o;
            end
         end
    end
end

%cluster components per electrode
ConnectedClusters = ConnectedComponents;
ConnectedClusters.NumObjects=size(oncomp,1);
ConnectedClusters.PixelIdxList=cell(1,size(oncomp,1));
for i=1:size(oncomp,1)
     curr_comp=9999;
     oncomp_el = oncomp(i,:);
     oncomp_el(oncomp_el==0)=[];
     for j=1:length(oncomp_el)
         new_comp=oncomp_el(j);
         if new_comp ~= curr_comp
         ConnectedClusters.PixelIdxList{i} = [ConnectedClusters.PixelIdxList{i}; ConnectedComponents.PixelIdxList{new_comp}];
         end
         curr_comp = new_comp;
     end
end

%% STEP 9E: SEGMENT ALL CONTACTS
% %get coordinates of electrode contacts
%begin met target en entry en reconstrueer contactpunten tussenin

allTidx = closest3(1:2:end,:); %all targets
allEidx = closest3(2:2:end,:); %all entries

%SE = strel("sphere",1);%for erosion
contact_table=cell(1,3);

%make gray binary
gray_bin = tpm.gray>0.5;

D=5.8; %voxel distance for 3.5mm electrode distance with 0.6mm cubic voxels


%for every electrode
for i=1:numel(ConnectedClusters.PixelIdxList)
    %make binary 3D with only one electrode
    binaryElectr = zeros(size(ct_coreg_threshold.anatomy));
    binaryElectr(ConnectedClusters.PixelIdxList{i}) = 1;

    %find separate contact
    ElectrodeComponents = bwconncomp(binaryElectr);
    
    % all indices of CT to one matrix, save original component
    allcomp = []; allcompidx=[];
    for j=1:length(ElectrodeComponents.PixelIdxList)
        allcomp = [allcomp; ElectrodeComponents.PixelIdxList{j}];
        allcompidx= [allcompidx; ones(size(ElectrodeComponents.PixelIdxList{j},1),1)*j];
    end

    %look component with closest target
    IdxT= sub2ind(size(ct_coreg_threshold.anatomy),allTidx(i,1),allTidx(i,2),allTidx(i,3));
    %[allcompX allcompY allcompZ] = ind2sub(size(ct_coreg_threshold.anatomy),allcomp);
    %[dist, idsort] = sort(sqrt( (allcompX-allTidx(i,1)).^2 + (allcompY-allTidx(i,2)).^2 + (allcompZ-allTidx(i,3)).^2));
    %closestContact = dsearchn(allcomp,IdxT);
    closestContact = find(allcomp==IdxT);
    [closestContactx closestContacty closestContactz]=ind2sub(size(ct_coreg_threshold.anatomy),ElectrodeComponents.PixelIdxList{allcompidx(closestContact)});

    %if closest contact is a cluster of contacts (oblique electrodes)
    if length(closestContactx)>40 
        % find 15 voxels closest to target
        % sort by Euclidean distance:
        [dist, idsort] = sort(sqrt( (closestContactx-allTidx(i,1)).^2 + (closestContacty-allTidx(i,2)).^2 + (closestContactz-allTidx(i,3)).^2));
        closestContact_new = ElectrodeComponents.PixelIdxList{allcompidx(closestContact)}(idsort(1:15));
        [closestContactx closestContacty closestContactz]=ind2sub(size(ct_coreg_threshold.anatomy),closestContact_new);
    end

    %get centroid of closest cluster target
    centrTX(i) = round(mean(closestContactx));
    centrTY(i) = round(mean(closestContacty));
    centrTZ(i) = round(mean(closestContactz));

    hold on; scatter3(centrTX(i),centrTY(i),centrTZ(i),30,'filled','k');

    %look for contact closest to planned entry 
    IdxE= sub2ind(size(ct_coreg_threshold.anatomy),allEidx(i,1),allEidx(i,2),allEidx(i,3));
    closestContact = dsearchn(allcomp,IdxE);
    [closestContactx closestContacty closestContactz]=ind2sub(size(ct_coreg_threshold.anatomy),ElectrodeComponents.PixelIdxList{allcompidx(closestContact)});

    %get centroid of closest cluster target
    centrEX(i) = round(mean(closestContactx));
    centrEY(i) = round(mean(closestContacty));
    centrEZ(i) = round(mean(closestContactz));

    hold on; scatter3(centrEX(i),centrEY(i),centrEZ(i),30,'filled','k');

    %get contacts by adding a contact every 3.5mm from the target to the
    %number of contacts in electrode_table

    %grey white for target
    idxC = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i), centrTY(i), centrTZ(i));
    %voxels around
    idxC1 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i)+1, centrTY(i), centrTZ(i));
    idxC2 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i)-1, centrTY(i), centrTZ(i));
    idxC3 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i), centrTY(i)+1, centrTZ(i));
    idxC4 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i), centrTY(i)-1, centrTZ(i));
    idxC5 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i), centrTY(i), centrTZ(i)+1);
    idxC6 = sub2ind(size(ct_coreg_threshold.anatomy),centrTX(i), centrTY(i), centrTZ(i)-1);
    if gray_bin(idxC)
        gray(1)=1;
    elseif sum(gray_bin([idxC idxC1 idxC2 idxC3 idxC4 idxC5 idxC6]))>2
        gray(1)=1;
    elseif sum(gray_bin([idxC idxC1 idxC2 idxC3 idxC4 idxC5 idxC6]))>0
        gray(1)=0.5;
    else
        gray(1)=0;
    end

    contact_table = [contact_table; [electrode_table{i} num2str(1)], {[centrTX(i) centrTY(i) centrTZ(i)]}, gray(1)]; %first contact of each electrode is target point

    lastcontactX = centrTX(i); lastcontactY = centrTY(i); lastcontactZ = centrTZ(i);
    for k=2:electrode_table{i,2}
        x0=lastcontactX; x1=centrEX(i);
        y0=lastcontactY; y1=centrEY(i);
        z0=lastcontactZ; z1=centrEZ(i);

        xyd = sqrt((x1-x0).^2+(y1-y0).^2+(z1-z0).^2);
        x = x0 + D .* (x1-x0)./xyd;
        y = y0 + D .* (y1-y0)./xyd;
        z = z0 + D .* (z1-z0)./xyd;
        hold on; scatter3(x,y,z,30,'filled','b')
        lastcontactX = x; lastcontactY = y; lastcontactZ = z;

        idxC = sub2ind(size(ct_coreg_threshold.anatomy),round(x), round(y), round(z));
        idxC1 = sub2ind(size(ct_coreg_threshold.anatomy),round(x+1), round(y), round(z));
        idxC2 = sub2ind(size(ct_coreg_threshold.anatomy),round(x-1), round(y), round(z));
        idxC3 = sub2ind(size(ct_coreg_threshold.anatomy),round(x), round(y+1), round(z));
        idxC4 = sub2ind(size(ct_coreg_threshold.anatomy),round(x), round(y-1), round(z));
        idxC5 = sub2ind(size(ct_coreg_threshold.anatomy),round(x), round(y), round(z+1));
        idxC6 = sub2ind(size(ct_coreg_threshold.anatomy),round(x), round(y), round(z-1));
        if gray_bin(idxC)
            gray(k)=1;
        elseif sum(gray_bin([idxC idxC1 idxC2 idxC3 idxC4 idxC5 idxC6]))>2
            gray(k)=1;
        elseif sum(gray_bin([idxC idxC1 idxC2 idxC3 idxC4 idxC5 idxC6]))>0
            gray(k)=0.5;
        else
            gray(k)=0;
        end

        contact_table = [contact_table; [electrode_table{i} num2str(k)], {[x y z]}, gray(k)];
    end
end

%remove first empty line of contact_table
contact_table(1,:)=[];

clearvars idxC idxC1 idxC2 idxC3 idxC4 idxC5 idxC6 x x0 x1 xyd y y0 y1 z z0 z1 z2 z3 v1 v2 n oncomp oncomp_el SL1 SL2 SL3 lastcontactX ...
    lastcontactY lastcontactZ i j k o p t e dist D curr_comp closestContactz closestContacty ...
    closestContactx closestContact closestContact_new IdxE IdxT ind_ct line3D 
%% STEP 9F
% to do transform contact_table into coordinates that can be used into the
% ct scan registered differently. Maybe use the other ct from the beginnen
% and transform with ct as basis instead of t1
%% STEP 10: select electrodes from ct - matlab
% the order in which you click electrodes does not matter. Just make sure
% you click all electrodes implanted!
clc
fprintf(' ----- OPEN THE CT-SCAN AND CLICK ON ALL ELECTRODES. YOU CAN CHECK WHETHER YOU HAVE ALL ELECTRODES BY CLICKING ON VIEW RESULT ----- \n')
fprintf(' ----- SAVE WHEN FINISHED IN %s \n',cfg(1).deriv_directory)

ctmr
% view result
% save image: saves as nifti hdr and img files
% this is saved as electrodes1.hdr and electrodes1.img

%% STEP 11: sort unprojected electrodes - matlab

fprintf('------ OPEN THE ELECTRODES.TSV \n-----')
% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg(1).elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% Make ieeg folder
if exist(cfg(1).ieeg_directory, 'dir')
    fprintf('%s exists already\n',cfg(1).ieeg_directory)
else
    mkdir(cfg(1).ieeg_directory)
end

fprintf('------ OPEN THE CLICKED ELECTRODES YOU SAVED IN THE PREVIOUS STEP \n-----')

% sort unprojected electrodes
cfg(1).saveFile = sprintf('%s%s_%s_electrodes_temp.mat',cfg(1).deriv_directory,cfg(1).sub_label,cfg(1).ses_label);
sortElectrodes(tb_elecs,cfg(1)); % [electrode labels, folder to save]
fprintf('Matched electrodes are saved in %s\n',cfg(1).saveFile)
% loads img file with electrodes from previous step
% saves in electrodes_temp.mat;
%% STEP 12: plot electrodes 2 surface - matlab - run for sEEG and ECoG
% corrects for the brain shift - ONLY for ECoG
% 1xN STRIP: do not project any electrodes that are already close to the
%               surfaces, such as subtemporal and interhemispheric
% DEPTH: do not project

% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg(1).elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% log_elec_incl = ~strcmp(tb_elecs.group,'other');
% tb_elecs = tb_elecs(log_elec_incl,:);

[filename, pathname] = uigetfile('*.mat','Select electrodes_temp.mat',cfg(1).deriv_directory);
load(fullfile(pathname,filename));

% only select letters from channelname
letters = regexp(tb_elecs.name,'[a-z_A-Z]');
channame = cell(size(tb_elecs,1),1);
for chan = 1:size(tb_elecs,1)
    channame{chan,1} = tb_elecs.name{chan}(letters{chan});
end

% load json-file with formats of specific electrode strips/grids
files = dir(cfg(1).elec_input);
jsonfile = find(contains({files(:).name},'_ieeg.json')==1,1);
if ~isempty(jsonfile)
    json = loadjson(fullfile(cfg(1).elec_input, files(jsonfile).name));
    chansplit = strsplit(json.iEEGElectrodeGroups,{';','['});
    changroups = chansplit(diff(contains(chansplit,']'))==1);
else
    changroups = [];
end

% electrodes on strip subtemporal or interhemispheric should be skipped
% because they are already located close to brain surface.
fprintf('Which electrodes should be skipped since they are already close to surfaces? \n')
fprintf('(such as subtemporal, interhemispheric, depth)? \n')
if ~isempty(changroups)
    fprintf('choose from %s %s %s %s %s %s %s %s %s', changroups{:});
end
skip_elec = input(': ','s');
skip_elec = strsplit(skip_elec,{', ',',',' '});
skip_elec = skip_elec(~cellfun(@isempty,skip_elec));

% determine format of each specific electrode strip/grid and correct for
% brainshift
if contains(json.iEEGElectrodeGroups,'seeg','IgnoreCase',1)
    disp('This session is seeg, so no correction for brain shift is needed')
    
    elecmatrix_shift = elecmatrix;
    
elseif contains(json.iEEGElectrodeGroups,'ecog','IgnoreCase',1) || contains(json.iEEGElectrodeGroups,'strip','IgnoreCase',1)
    
    format = strsplit(json.iEEGElectrodeGroups,';');
    elecmatrix_shift = NaN(size(elecmatrix));
    
    for i=1:size(format,2)
        if contains(format{i},'[')
            % find channame and format of this group of electrodes
            formatelec = strsplit(format{i},{'[','x',']'});
            formatelec = formatelec(~cellfun(@isempty,formatelec));
            
            % find which electrodes belong to this specific group
            num_elecs = find(strcmp(formatelec{1},channame)==1);
            num_elecs = num_elecs(~isnan(elecmatrix(num_elecs,1)) );
            
            % check whether all electrodes in group are included
            if str2double(formatelec{2}) * str2double(formatelec{3}) == numel(num_elecs)
            else
                disp('ERROR: mismatch between found electrodes and expected number of electrodes in format!')
            end
            
            % skip electrodes mentioned above
            if ~contains(format{i},skip_elec)
                % format of specific grid/strip
                if any(contains(formatelec,'1')) % it is a 1xN strip/depth electrode
                    settings = [0,2];
                elseif any(contains(formatelec,'2')) % it is a 2xN strip
                    settings = [4,1];
                else % it is a grid (larger than 2xN)
                    settings = [5,1];
                end
                % correct location of electrodes for brain shift
                [out_els,out_els_ind] = electrodes2surf(...
                    cfg(1).sub_label,... % subject name
                    settings(1),... % 5 for grid, 4 for 2xN strip, 0 for 1xN strip
                    settings(2),... % 1 for grid or 2xN strip, 2 for 1xN strip
                    num_elecs,... % matrix indices (rows) in elecmatrix, e.g. 1:64 is for grid C1-C64
                    [cfg(1).deriv_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_electrodes_temp.mat'],... % file that contains elecmatrix
                    [cfg(1).deriv_directory, cfg(1).sub_label,'_',cfg(1).hemisphere{1},'_surface',num2str(k),'_',num2str(settings_hull(1)),'_0',num2str(settings_hull(2)*10),'.img'],... % hull we just created
                    [cfg(1).deriv_directory, cfg(1).sub_label,'_',cfg(1).hemisphere{1},'_surface',num2str(k),'_',num2str(settings_hull(1)),'_0',num2str(settings_hull(2)*10),'.img'],... % hull we just created
                    cfg(1).deriv_directory);
%          [cfg(1).anat_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii'],... % T1 file (mr.img for same image space with electrode positions)

                % saves automatically a matrix with projected electrode positions and an image
                % with projected electrodes
                % saves as electrodes_onsurface_filenumber_inputnr2
                elecmatrix_shift(num_elecs,:) = out_els;
                
            elseif contains(format{i},skip_elec)
                elecmatrix_shift(num_elecs,:) = elecmatrix(num_elecs,:);
            end
            
        end
    end
end

tb_elecs.x = elecmatrix_shift(:,1);
tb_elecs.y = elecmatrix_shift(:,2);
tb_elecs.z = elecmatrix_shift(:,3);


%% STEP 13: save electrode positions, corrected for brain shift to electrodes.tsv - matlab

saveFile = sprintf('%s%s_%s_electrodesPositioned.tsv',cfg(1).deriv_directory,cfg(1).sub_label,cfg(1).ses_label);
writetable(tb_elecs, saveFile, 'Delimiter', 'tab', 'FileType', 'text');
fprintf('Electrode positions, corrected for brainshift, are saved in %s\n',saveFile)

%% START HERE WHEN electrodes were overwritten

fprintf('------ OPEN THE ELECTRODES.TSV \n-----')
% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electrocesPositioned.tsv file',cfg(1).deriv_directory);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

elecmatrix_shift(:,1) = tb_elecs.x;
elecmatrix_shift(:,2) = tb_elecs.y;
elecmatrix_shift(:,3) = tb_elecs.z;

%% STEP 14: Write electrode positions as numbers in a nifti - matlab
% This is not necessary, only if you want to do some extra checks or so.
% e.g. it can be nice to visualize the projected electrodes in MRIcron.
% to visualize in MRIcron: open the defaced MRI and add the surface_all.img
% as overlay % does not work, cannot open the .img as overlay

[output,els,els_ind,outputStruct] = position2reslicedImage(elecmatrix_shift,[cfg(1).anat_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii']);

for filenummer=1:100
    save([cfg(1).deriv_directory cfg(1).sub_label '_' cfg(1).ses_label,'_electrodes_surface_loc_all' int2str(filenummer) '.mat'],'elecmatrix_shift');
    outputStruct.fname=[cfg(1).deriv_directory,cfg(1).sub_label,'_',cfg(1).ses_label,'_electrodes_surface_all' int2str(filenummer) '.img' ];
    if ~exist(outputStruct.fname,'file')>0
        fprintf('----- SAVING %s ------ \n', outputStruct.fname);
        % save the data
        spm_write_vol(outputStruct,output);
        break
    end
end


%% STEP 15: Convert the .gii coordinates to the MRI native space - matlab

for i=1:size(cfg(1).hemisphere,2)
    % load the Freesurfer gifti (freesurfer coordinates)
    g = gifti(fullfile(cfg(1).freesurfer_directory,'surf',[cfg(1).hemisphere{i},'h.pial.gii']));
    
    % convert from freesurfer space to original space
    % the transformation matrix is in the /freesurfer/sub/mri/orig.mgz file:
    mri_orig = fullfile(cfg(1).freesurfer_directory,'mri','orig.mgz');
    orig = MRIread_windows(mri_orig); % MRIread is a function from vistasoft (freesurfer)
    Torig = orig.tkrvox2ras;
    Norig = orig.vox2ras;
    freeSurfer2T1 = Norig*inv(Torig); %#ok<MINV>
    
    % convert freesurfer vertices to original T1 space
    vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
    vert_mat = freeSurfer2T1*vert_mat;
    vert_mat(4,:) = [];
    vert_mat = vert_mat';
    g.vertices = vert_mat; clear vert_mat
    
    % save correct coordinates back as a gifti
    gifti_name = fullfile(cfg(1).surface_directory, ...
        [cfg(1).sub_label,'_',cfg(1).ses_label,'_T1w_pial.' cfg(1).hemisphere{i} '.surf.gii']);
    
    save(g,gifti_name,'Base64Binary')
    
    fprintf('gifti %s converted to original space \n',cfg(1).hemisphere{i})
end

%% STEP 16: add labels of atlases to tsv-file - matlab

[tb_elecs_atlases, cfg(1).destrieux_labels, cfg(1).DKT_labels] = lookupAtlases(cfg(1),tb_elecs);

disp('Atlases added')

%% STEP 17: CHECK ATLAS WITH ELECTRODE POSITIONS - matlab
close all

cfg(1).show_labels = 'yes';
cfg(1).change_size = 'no';
cfg(1).change_color = 'no';
cfg(1).view_atlas ='yes';
cfg(1).elec_offset = 'yes';
cfg(1).atlas = 'Destrieux'; % [DKT/Destrieux]
% cfg(1).atlas = 'DKT'; % [DKT/Destrieux]
cfg(1).view_elec ='yes';
cfg(1).save_fig = 'no';

check_atlas_elec_MRI(cfg(1),tb_elecs_atlases)

%% REPLACE NAN WITH N/A

tb_elecs_atlases = bids_tsv_nan2na(tb_elecs_atlases);

%% STEP 18: save electrodes.tsv, make electrodes descriptor, and add hemisphere to existing ieeg_json files

addpath(cfg(1).fieldtrip_folder)
addpath(cfg(1).fieldtrip_private)
ft_defaults

writetable(tb_elecs_atlases, ...
    fullfile(cfg(1).ieeg_directory, ...
    [cfg(1).sub_label '_' cfg(1).ses_label '_electrodes.tsv' ]),...
    'Filetype','text','Delimiter','\t');

fprintf('Saved %s\n',fullfile(cfg(1).ieeg_directory, ...
    [cfg(1).sub_label '_' cfg(1).ses_label '_electrodes.tsv' ]))

% 5. write T1w.json accompanying the .nii file
mri_name = [cfg(1).anat_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii'];
create_json_mri(replace(mri_name,'nii','json'))

% 6. create electrodes descriptor

create_elecDesc(cfg(1).proj_diroutput,cfg(1))

% 7. write coordsystem
    
write_coordsystemJSON(cfg(1))
 
% 8. add hemisphere to ieeg_json files

D = dir(cfg(1).ieeg_directory);
ieeg_json_filenums = contains({D(:).name},'_ieeg.json');

for i=1:size(D,1)
    if ieeg_json_filenums(i) == 1
        ieeg_json = read_json([D(i).folder '/' D(i).name] );
        if size(cfg(1).hemisphere,2) == 1
            if strcmpi(cfg(1).hemisphere,'r')
                ieeg_json.iEEGPlacementScheme = 'right';
            elseif strcmpi(cfg(1).hemisphere,'l')
                ieeg_json.iEEGPlacementScheme = 'left';
            end
        elseif size(cfg(1).hemisphere,2) == 2
            ieeg_json.iEEGPlacementScheme = 'left,right';
        end
        write_json([D(i).folder '/' D(i).name], ieeg_json)
    end
end

%% STEP 19: save everything to second proj_diroutput directory, not used at the moment

if~isempty(cfg(2).proj_diroutput)
    cfg(1).freesurfer_directory = sprintf('%sderivatives/freesurfer/%s/%s/',cfg(2).proj_diroutput,cfg(1).sub_label,cfg(1).ses_label);
    cfg(1).anat_directory = sprintf('%s%s/%s/anat/',cfg(2).proj_diroutput,cfg(1).sub_label,cfg(1).ses_label);
    cfg(1).ieeg_directory = sprintf('%s%s/%s/ieeg/',cfg(2).proj_diroutput,cfg(1).sub_label,cfg(1).ses_label);
    cfg(1).surface_directory = sprintf('%sderivatives/surfaces/%s/%s/',cfg(2).proj_diroutput,cfg(1).sub_label,cfg(1).ses_label);
    cfg(1).elec_input = sprintf('%s%s/%s/ieeg/',cfg(2).proj_diroutput,cfg(1).sub_label,cfg(1).ses_label);
    
    % 5B. save electrodes.tsv
    writetable(tb_elecs_atlases, ...
        fullfile(cfg(1).ieeg_directory, ...
        [cfg(1).sub_label '_' cfg(1).ses_label '_electrodes.tsv' ]),...
        'Filetype','text','Delimiter','\t');
    
    disp('Saved electrodes.tsv')
    
    % 5. write T1w.json accompanying the .nii file
    if ~exist(cfg(1).anat_directory,'dir')
        mkdir(fullfile(cfg(1).anat_directory))
    end
    mri_name = [cfg(1).anat_directory, cfg(1).sub_label,'_',cfg(1).ses_label,'_rec-deface_T1w.nii'];
    create_json_mri(replace(mri_name,'nii','json'))
    
    % 6. create electrodes descriptor
    
    create_elecDesc(cfg(2).proj_diroutput,cfg)
    
    % 7. write coordsystem
    
    write_coordsystemJSON(cfg(1))
    
    % 8. add hemisphere to ieeg_json files
    
    D = dir(cfg(1).ieeg_directory);
    ieeg_json_filenums = contains({D(:).name},'_ieeg.json');
    
    for i=1:size(D,1)
        if ieeg_json_filenums(i) == 1
            ieeg_json = read_json([D(i).folder '/' D(i).name] );
            if size(cfg(1).hemisphere,2) == 1
                if strcmpi(cfg(1).hemisphere,'r')
                    ieeg_json.iEEGPlacementScheme = 'right';
                elseif strcmpi(cfg(1).hemisphere,'l')
                    ieeg_json.iEEGPlacementScheme = 'left';
                end
            elseif size(cfg(1).hemisphere,2) == 2
                ieeg_json.iEEGPlacementScheme = 'left,right';
            end
            write_json([D(i).folder '/' D(i).name], ieeg_json)
        end
    end
    
end


%% FUNCTIONS
close
function json = read_json(filename)
ft_info('reading %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    json = loadjson(filename);
else
    fid = fopen(filename, 'r');
    str = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
    json = jsondecode(str);
end
end

function write_json(filename, json)
json = remove_empty(json);
ft_info('writing %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    savejson('', json, filename);
else
    str = jsonencode(json);
    fid = fopen(filename, 'w');
    fwrite(fid, str);
    fclose(fid);
end
end

function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);
end





