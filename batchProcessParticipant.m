app = processC3D_exported;
app.path = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\Juliana\09_05_2023';
app.emgLabelCSVPath = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\emg_labels_d2_v2.csv';
fileList = dir([app.path '\*.c3d']);

for i = 21%1: numel(fileList)
    app.c3dFileName = fileList(i).name;
    disp(['processing' app.c3dFileName]);
    
    app.ClimbingCheckBox.Value = 1;
    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
    app.detectforceplatesautomaticallyCheckBox.Value = 0;

    app.filtermarkersCheckBox.Value = 0;
    app.filterGRFsCheckBox.Value = 0;
    app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;

    app.processC3Dfile();

    app.CreatefilesButtonPushed(0);
end
