app = processC3D_exported;
app.path = 'C:\insidebone\10_Drive\Willi\Study_GaitModifications_MechanicalStimuli\GaitData_12Hz_eachStepSeperately\TD31';
app.emgLabelCSVPath = 'C:\insidebone\10_Drive\Willi\Study_GaitModifications_MechanicalStimuli\emg_labels.csv';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

for i = 1 : numel(fileList)
    if ~contains(fileList(i).folder, 'TD01') && ~strcmp(fileList(i).name, 'c3dfile.c3d') && ~strcmp(fileList(i).name, 'tmp.c3d') %&& isfolder(fullfile(fileList(i).folder, strrep(fileList(i).name, '.c3d', '')))
        %
        % pause(0.5);
        app.path = fileList(i).folder;
        copyfile(fullfile(app.path, fileList(i).name), fullfile(app.path, 'tmp.c3d'));

        % pause(0.5);
        app.c3dFileName = 'tmp.c3d';
        disp(['processing' fileList(i).name]);

        % app.ClimbingCheckBox.Value = 0;
        app.AMTITandemTreadmillCheckBox.Value = 1;
        app.AMTITandemTreadmillCheckBoxValueChanged(0);
        % app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 1;
        % app.detectforceplatesautomaticallyCheckBox.Value = 1;
        app.filterorderSpinner.Value = 4;
        app.cutofffrequencySpinner.Value = 12;
        app.renameEMGlabelsCheckBox.Value = 1;
        % app.normalizeEMGCheckBox.Value = 1;

        % app.filtermarkersCheckBox.Value = 0;
        % app.filterGRFsCheckBox.Value = 1;
        % app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 1;
        % app.filterEMGandexporttostoCheckBox.Value = 0;

        app.processC3Dfile();
        pause(0.2);

        stepCountLeft = numel(app.leftListBox.Items);
        stepCountRight = numel(app.rightListBox.Items);

        for stepLeft = 1 : stepCountLeft
            app.processC3Dfile();
            pause(0.2);
            app.leftListBox.Value = app.leftListBox.Items(stepLeft);
            app.rightListBox.Value = {};
            app.filterGRFs();
            app.updateCycles();
            app.CreatefilesButtonPushed(0);
            pause(0.2);
            copyfile(fullfile(app.path, 'tmp'), fullfile(app.path, [strrep(fileList(i).name, '.c3d', ''), '_left_' num2str(stepLeft)]))
            % copyfile(fullfile(app.path, 'tmp', 'EMG_filtered.sto'), fullfile(app.path, [strrep(fileList(i).name, '.c3d', ''), '_left_' num2str(stepLeft)], 'EMG_filtered_normalizedTo1.sto'))
        end

        for stepRight = 1 : numel(app.rightListBox.Items)
            app.processC3Dfile();
            pause(0.2);
            app.rightListBox.Value = app.rightListBox.Items(stepRight);
            app.leftListBox.Value = {};
            app.filterGRFs();
            app.updateCycles();
            app.CreatefilesButtonPushed(0);
            pause(0.2);
            copyfile(fullfile(app.path, 'tmp'), fullfile(app.path, [strrep(fileList(i).name, '.c3d', ''), '_right_' num2str(stepRight)]))
            % copyfile(fullfile(app.path, 'tmp', 'EMG_filtered.sto'), fullfile(app.path, [strrep(fileList(i).name, '.c3d', ''), '_right_' num2str(stepRight)], 'EMG_filtered_normalizedTo1.sto'))
        end

        delete(fullfile(app.path, 'tmp', '*'));
        rmdir(fullfile(app.path, 'tmp'));
        delete(fullfile(app.path, 'tmp.c3d'));
    end
end
