app = processC3D_exported;
app.path = 'pathToC3DFiles';
% app.emgLabelCSVPath = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\emg_labels_d2_v2.csv';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

for i = 1 : numel(fileList)
    app.path = fileList(i).folder;
    app.c3dFileName = fileList(i).name;
    if ~strcmp(app.c3dFileName, 'c3dfile.c3d')
        disp(['processing ' app.path(end-20 : end) ' - ' app.c3dFileName]);

        % app.ClimbingCheckBox.Value = 1;
        app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
        app.rotatearoundyaxisDropDown.Value = '270°';
        app.detectforceplatesautomaticallyCheckBox.Value = 1;
        app.Ignorec3deventsCheckBox.Value = 0;
        app.doPostZeroLevellingCheckBox.Value = 1;
        app.validfootstrikeeventsLabel.Visible = "on";
        app.forceplatecontactsLabel.Visible = "off";

        app.filtermarkersCheckBox.Value = 0;
        app.filterEMGandexporttostoCheckBox.Value = 0;
        app.filterGRFsCheckBox.Value = 1;

        app.processC3Dfile();

        app.CreatefilesButtonPushed(0);
    end
end
