app = processC3D_exported;
app.path = 'C:\StudentSynology\BBAM_SS2025\Sandro\2026_04_28';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

% app.emgLabelCSVPath = 'C:\insidebone\10_Drive\Willi\TU_Metastasen\Data_v2\emg_labels.csv';

for i = 11 : numel(fileList)
    if ~contains(fileList(i).name, '_reduced')
        % if ~(contains(fileList(i).name, 'cmj') || contains(fileList(i).name, 'stair') || contains(fileList(i).name, 'static'))

            app.automaticMovementtypeCheckBox.Value = 0;
            app.detectforceplatesautomaticallyCheckBox.Value = 0;
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
            app.exporteachstepseparatelyCheckBox.Value = 0;
            app.croptovalidfootstrikesbufferCheckBox.Value = 0;

            % EMG data processing
            app.filterEMGCheckBox.Value = true;
            app.normalizeallmagnitudesto1CheckBox.Value = 0;
            app.EMGfactorSpinner.Value = 1;
            app.exportrawEMGCheckBox.Value = 1;
            app.renameEMGlabelsCheckBox.Value = 0;

            app.path = fileList(i).folder;
            app.c3dFileName = fileList(i).name;
            disp(['processing ' app.c3dFileName]);



            app.processC3Dfile();

            app.CreatefilesButtonPushed(0);
        % end
    end
end
