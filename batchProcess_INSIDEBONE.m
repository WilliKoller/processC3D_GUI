app = processC3D_exported;
app.path = 'C:\insidebone\10_Drive\Willi\TU_Metastasen\Data_v2';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

app.emgLabelCSVPath = 'C:\insidebone\10_Drive\Willi\TU_Metastasen\Data_v2\emg_labels.csv';

for i = 35 : numel(fileList)
    if ~contains(fileList(i).name, '_reduced')
        if ~(contains(fileList(i).name, 'cmj') || contains(fileList(i).name, 'stair') || contains(fileList(i).name, 'static'))

            app.path = fileList(i).folder;
            app.c3dFileName = fileList(i).name;
            app.checkMovementType();
            disp(['processing ' app.c3dFileName]);

            % EMG data processing
            app.filterEMGCheckBox.Value = true;
            app.normalizeallmagnitudesto1CheckBox.Value = 0;
            app.EMGfactorSpinner.Value = 1;
            app.exportrawEMGCheckBox.Value = 0;
            app.renameEMGlabelsCheckBox.Value = 1;

            app.exporteachstepseparatelyCheckBox.Value = true;
            app.croptovalidfootstrikesbufferCheckBox.Value = true;
            app.BufferStartsSpinner.Value = 0.1;
            app.BufferEndsSpinner.Value = 0.1;

            app.processC3Dfile();

            app.CreatefilesButtonPushed(0);
        end
    end
end
