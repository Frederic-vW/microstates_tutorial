% microstate analysis (including mstsa) on sleep data

% delete workspace
clear all
close all
clc
% make sure we get the correct relative paths, starting in the folder 
% of this script
cd(fileparts(matlab.desktop.editor.getActiveFilename))

% results directory
SavePath = '.\results';
if ~exist(SavePath, 'dir')
    mkdir(SavePath)
end

freq_lo = 1; % 1, 2
freq_hi = 30; % 20, 30
n_coeffs = 2000;
apply_bandpass = true;

% fitting paramters (decide whether GFP peaks only or not)
FitPars = struct();
FitPars.nClasses = 4;
FitPars.lambda= 5;
FitPars.b = 12; % 20
FitPars.PeakFit = true; % interpolation y/n?
FitPars.BControl = true;
FitPars.Rectify = false;
FitPars.Normalize = false;

% clustering parameters
ClustPars = struct();
ClustPars.MinClasses = 4;
ClustPars.MaxClasses = 4;
ClustPars.GFPPeaks = true;
ClustPars.IgnorePolarity = true;
ClustPars.MaxMaps = inf;
ClustPars.Restarts = 20;
ClustPars.UseAAHC = false;
ClustPars.Normalize = false;
disp('[+] Paths, fitting and clustering parameters set up')

% % Read the data
eeglabpath = fileparts(which('eeglab.m'));
DipFitPath = fullfile(eeglabpath,'plugins','dipfit');
%elpFile = fullfile(DipFitPath,'standard_BESA','standard-10-5-cap385.elp');
eeglab

% collect files
sleep_stages = {'W', 'N3'};
nGroups = length(sleep_stages);
data_dir = '..\data';

AllSubjects = [];
for Group = 1:nGroups
    GroupIndex{Group} = [];
    % all sleep files are in one directory
    query = ['*_', sleep_stages{Group}, '.set'];
    files_ = dir(fullfile(data_dir, query));
    FileNamesGroup = {files_.name};
    for f = 1:numel(FileNamesGroup)
        tmpEEG = pop_loadset('filename',FileNamesGroup{f},'filepath',data_dir);
        [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off');
        tmpEEG = pop_reref(tmpEEG, []); % apply average reference
        if apply_bandpass
            tmpEEG = pop_eegfilt(tmpEEG, freq_lo, freq_hi, n_coeffs);
        end
        % set the group (appears in the statistics output)
        tmpEEG.group = sprintf('Group_%i',Group);
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % store
        GroupIndex{Group} = [GroupIndex{Group} CURRENTSET]; % track the group
        AllSubjects = [AllSubjects CURRENTSET];
    end
end
eeglab redraw
disp('[+] Datasets loaded, re-referenced (avg.) and filtered')
   
% % subject-wise clustering
for i = 1:numel(AllSubjects) 
    EEG = eeg_retrieve(ALLEEG,AllSubjects(i)); % subject EEG
    fprintf(1,'Clustering dataset %s (%i/%i)\n',EEG.setname,i,numel(AllSubjects ));
    EEG = pop_FindMSTemplates(EEG, ClustPars); % mod. kmeans clustering
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, AllSubjects (i));
end
eeglab redraw
disp('[+] Subject-level clustering complete')

% % group-level microstate maps

% normative maps to orient us later
%templatepath = fullfile(fileparts(which('eegplugin_Microstates.m')),'Templates');

%EEG = pop_loadset('filename',...
%    'Normative microstate template maps Neuroimage 2002.set',...
%    'filepath', templatepath);
%[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');

% And we have a look at it
%NormativeTemplateIndex = CURRENTSET;
%pop_ShowIndMSMaps(ALLEEG(NormativeTemplateIndex), 4); 
%drawnow;

% group-level maps
for Group = 1:nGroups
    EEG = pop_CombMSTemplates(ALLEEG, GroupIndex{Group}, 0, 0, ...
        sprintf('GrandMean Group %i',Group));
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off');
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    GrandMeanIndex(Group) = CURRENTSET;
end

% across-group microstate maps (grand-grand mean)
if nGroups > 1
    EEG = pop_CombMSTemplates(ALLEEG, GrandMeanIndex, 1, 0, 'GrandGrandMean');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off');
    GrandGrandMeanIndex = CURRENTSET;
else
    GrandGrandMeanIndex = GrandMeanIndex(1);
end

% sort the grand-grand-mean according to literature templates
%[ALLEEG,EEG] = pop_SortMSTemplates(ALLEEG, GrandGrandMeanIndex, 1, NormativeTemplateIndex);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, GrandGrandMeanIndex);

% interactive re-ordering option
pop_ShowIndMSMaps(EEG, 4, GrandGrandMeanIndex, ALLEEG); 
eeglab redraw
disp('[+] Group-level clustering complete')

% sort group and subject maps wrt grand-grand-mean
% (1) sort group maps for group 1 and 2 wrt grand-grand-mean
if nGroups > 1
    ALLEEG = pop_SortMSTemplates(ALLEEG, GrandMeanIndex, 1, GrandGrandMeanIndex);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
end

% (2) sort subjects maps wrt (sorted) group maps
for Group = 1:nGroups
    ALLEEG = pop_SortMSTemplates(ALLEEG, GroupIndex{Group}, 0, GrandMeanIndex(Group));
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % store
end
eeglab redraw
disp('[+] Microstate maps sorted wrt normative maps')

% save re-sorted maps and EEGs
for f = 1:numel(ALLEEG)
    EEG = eeg_retrieve(ALLEEG,f);
    %pop_saveset(EEG, 'filename', EEG.setname, 'filepath', SavePath);
    pop_saveset(EEG, 'filename', EEG.filename, 'filepath', SavePath);
end
disp('[+] Datasets re-sorted wrt group maps')

% Visualize first data set
%EEG = eeg_retrieve(ALLEEG,1); 
%pop_ShowIndMSDyn([],EEG,0,FitPars);
%pop_ShowIndMSMaps(EEG,FitPars.nClasses);

% % Compute microstate dynamics (using the grand grand mean template)
pop_QuantMSTemplates(ALLEEG, AllSubjects, 1, FitPars, ...
    GrandGrandMeanIndex, ...
    fullfile(SavePath,'ResultsFromGrandGrandMeanTemplate.xlsx'));
disp('[+] Microstate dynamics quantified.')
