app = processC3D_exported;
app.path = 'C:\Users\willi\ucloud\Students\Judith_climbing\Pilotdata';
app.emgLabelCSVPath = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\emg_labels_d2_v2.csv';
fileList = dir([app.path '\*.c3d']);

for i = 1: numel(fileList)
    app.c3dFileName = fileList(i).name;
    disp(['processing' app.c3dFileName]);
    
    app.ClimbingCheckBox.Value = 1;
    app.invertforcesensorsCheckBox.Value = 1;
    app.mirrorCheckBox.Value = 0;
    app.rotatearoundyaxisDropDown.Value = '180°';

    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
    app.detectforceplatesautomaticallyCheckBox.Value = 0;
    app.filtermarkersCheckBox.Value = 0;
    app.filterGRFsCheckBox.Value = 0;
    app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;

    app.filterEMGandexporttostoCheckBox.Value = 1;
    app.normalizeEMGCheckBox.Value = 0;

    % CHANGE
    app.renameEMGlabelsCheckBox.Value = 0;

    app.processC3Dfile();

    app.CreatefilesButtonPushed(0);

    app.processC3Dfile();

    app.mirrorCheckBox.Value = 1;

    app.CreatefilesButtonPushed(0);
end
