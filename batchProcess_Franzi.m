app = processC3D_exported;
app.path = 'C:\StudentSynology\Franziska_Fussball\lotta_aufnahmen_3\Gait';
app.path = 'C:\StudentSynology\Franziska_Fussball\lotta_aufnahmen_3\instep_new';
app.path = 'C:\StudentSynology\Franziska_Fussball\lotta_aufnahmen_3\passkick';
% app.path = 'C:\StudentSynology\Franziska_Fussball\lotta_aufnahmen_3\EMG_MVC_trials';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

% app.emgLabelCSVPath = 'C:\insidebone\10_Drive\Willi\TU_Metastasen\Data_v2\emg_labels.csv';

for i = 1 : numel(fileList)
    if ~contains(fileList(i).name, '_reduced')
        % if ~(contains(fileList(i).name, 'cmj') || contains(fileList(i).name, 'stair') || contains(fileList(i).name, 'static'))

            app.path = fileList(i).folder;
            app.c3dFileName = fileList(i).name;
            % app.checkMovementType();
            disp(['processing ' app.c3dFileName]);

            % EMG data processing
            app.filterEMGCheckBox.Value = true;
            app.normalizeallmagnitudesto1CheckBox.Value = 0;
            app.EMGfactorSpinner.Value = 1;
            app.exportrawEMGCheckBox.Value = 1;
            app.renameEMGlabelsCheckBox.Value = 0;

            app.detectforceplatesautomaticallyCheckBox.Value = 0;
            app.exporteachstepseparatelyCheckBox.Value = 0;
            app.croptovalidfootstrikesbufferCheckBox.Value = 0;

            app.processC3Dfile();

            app.CreatefilesButtonPushed(0);
        % end
    end
end
