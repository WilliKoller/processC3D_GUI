classdef processC3D_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PrepareC3DfileUIFigure          matlab.ui.Figure
        Image                           matlab.ui.control.Image
        EMGSettingsPanel                matlab.ui.container.Panel
        EMGfactorSpinner                matlab.ui.control.Spinner
        EMGfactorSpinnerLabel           matlab.ui.control.Label
        exportrawEMGCheckBox            matlab.ui.control.CheckBox
        normalizeallmagnitudesto1CheckBox  matlab.ui.control.CheckBox
        renameEMGlabelsCheckBox         matlab.ui.control.CheckBox
        createCSVandopenButton          matlab.ui.control.Button
        SelectexistingCSVButton         matlab.ui.control.Button
        filterEMGCheckBox               matlab.ui.control.CheckBox
        TimingSettingsPanel             matlab.ui.container.Panel
        BufferEndsSpinner               matlab.ui.control.Spinner
        BufferEndsSpinnerLabel          matlab.ui.control.Label
        BufferStartsSpinner             matlab.ui.control.Spinner
        BufferStartsLabel               matlab.ui.control.Label
        croptovalidfootstrikesbufferCheckBox  matlab.ui.control.CheckBox
        exporteachstepseparatelyCheckBox  matlab.ui.control.CheckBox
        SetbeforeopeningC3DPanel        matlab.ui.container.Panel
        automaticMovementtypeCheckBox   matlab.ui.control.CheckBox
        MovementTypeDropDown            matlab.ui.control.DropDown
        MovementTypeDropDownLabel       matlab.ui.control.Label
        Ignorec3deventsCheckBox         matlab.ui.control.CheckBox
        AMTITandemTreadmillCheckBox     matlab.ui.control.CheckBox
        deleteoriginalfileafterrenamingCheckBox  matlab.ui.control.CheckBox
        removeblanksetcfromfilenameCheckBox  matlab.ui.control.CheckBox
        ClimbingPanel                   matlab.ui.container.Panel
        mirrorCheckBox                  matlab.ui.control.CheckBox
        invertforcesensorsCheckBox      matlab.ui.control.CheckBox
        ClimbingCheckBox                matlab.ui.control.CheckBox
        cutofffrequencySpinner_Marker   matlab.ui.control.Spinner
        cutofffrequencySpinner_2Label   matlab.ui.control.Label
        filterorderSpinner_Marker       matlab.ui.control.Spinner
        filterorderSpinnerLabel_Marker  matlab.ui.control.Label
        filtermarkersCheckBox           matlab.ui.control.CheckBox
        forceplatecontactsLabel         matlab.ui.control.Label
        v20Label                        matlab.ui.control.Label
        detectFPbyMarker                matlab.ui.control.EditField
        detectWalkingByMarker           matlab.ui.control.EditField
        cutofffrequencySpinner          matlab.ui.control.Spinner
        cutofffrequencySpinnerLabel     matlab.ui.control.Label
        filterorderSpinner              matlab.ui.control.Spinner
        filterorderSpinnerLabel         matlab.ui.control.Label
        LogTextArea                     matlab.ui.control.TextArea
        LogTextAreaLabel                matlab.ui.control.Label
        openselectedfilewithdefaultprogramButton  matlab.ui.control.Button
        rightListBox                    matlab.ui.control.ListBox
        rightListBoxLabel               matlab.ui.control.Label
        validfootstrikeeventsLabel      matlab.ui.control.Label
        leftListBox                     matlab.ui.control.ListBox
        leftListBoxLabel                matlab.ui.control.Label
        CreatefilesButton               matlab.ui.control.Button
        rotatearoundyaxisDropDown       matlab.ui.control.DropDown
        rotatearoundyaxisDropDownLabel  matlab.ui.control.Label
        detectwalkingdirectionautomaticallybymarkerCheckBox  matlab.ui.control.CheckBox
        filterGRFsCheckBox              matlab.ui.control.CheckBox
        detectforceplatesautomaticallyCheckBox  matlab.ui.control.CheckBox
        SelectaC3DfileLabel             matlab.ui.control.Label
        PrepareC3DfileforOpenSimsimulationsLabel  matlab.ui.control.Label
        Selectc3dfileButton             matlab.ui.control.Button
        UIAxesGRFs                      matlab.ui.control.UIAxes
        UIAxesSteps                     matlab.ui.control.UIAxes
    end

    %-------------------------------------------------------------------------%
    % Copyright (c) 2022 Koller W.                                            %
    %    Author:   Willi Koller,  2022                                        %
    %    email:    willi.koller@univie.ac.at                                  %
    % ----------------------------------------------------------------------- %

    % Process C3D (c) by Willi Koller, University of Vienna
    %
    % Process C3D is licensed under a
    % Creative Commons Attribution-NonCommercial 4.0 International License.
    %
    % You should have received a copy of the license along with this
    % work. If not, see <http://creativecommons.org/licenses/by-nc/4.0/>.
    % This package uses the btk toolbox
    % http://biomechanical-toolkit.github.io/docs/Wrapping/Matlab/_tutorial.html


    properties (Access = public)
        c3dFileName % Description
        path % Description
        rightFootStrikeList % Description
        leftFootStrikeList % Description
        c3d % Description
        grf_forces % Description
        markers % Description
        rotationAngle % Description
        firstFrame_ % Description
        leftStepForcePlateAssignment % Description
        leftStepForcePlateAssignmentAll % Description
        rightStepForcePlateAssignment % Description
        rightStepForcePlateAssignmentAll % Description
        acq % Description
        leftFootOffList % Description
        rightFootOffList % Description
        frequency_ % Description
        forceOffset % Description
        forceFrequency % Description
        grf_forces_adjusted % Description
        fileLoaded % Description
        EMG % Description
        EMG_lowPass % Description
        emgLabelCSVPath % Description
    end

    methods (Access = public)

        function writeLog(app, msg)
            app.LogTextArea.Value = [app.LogTextArea.Value(1:end-1); msg; ' '];
            scroll(app.LogTextArea, 'bottom');
            % drawnow; % Force text area to update
        end

        function updateCyclePlot(app)
            cla(app.UIAxesSteps);

            if ~isempty(app.c3d) && app.c3d.getNumForces() > 0
                selectedLeft = app.leftListBox.Value;
                for i = 1 : size(app.leftFootStrikeList, 2) - 1
                    footOffFrame = app.leftFootOffList( find( app.leftFootOffList > str2double(app.leftFootStrikeList{i}), 1 ) );
                    if ismember(app.leftFootStrikeList{i}, selectedLeft)
                        if ~isempty(footOffFrame)
                            fill(app.UIAxesSteps, [str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i}) footOffFrame str2double(app.leftFootStrikeList{i + 1})], [0 1 0.8 0], 'r', 'FaceAlpha', 1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        else
                            fill(app.UIAxesSteps, [str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i + 1})], [0 1 0], 'r', 'FaceAlpha', 1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        end
                    else
                        if ~isempty(footOffFrame)
                            fill(app.UIAxesSteps, [str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i}) footOffFrame str2double(app.leftFootStrikeList{i + 1})], [0 1 0.8 0], 'k', 'FaceAlpha', 0.1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        else
                            fill(app.UIAxesSteps, [str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i}) str2double(app.leftFootStrikeList{i + 1})], [0 1 0], 'k', 'FaceAlpha', 0.1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        end
                    end
                end
                selectedRight = app.rightListBox.Value;
                for i = 1 : size(app.rightFootStrikeList, 2) - 1
                    footOffFrame = app.rightFootOffList( find( app.rightFootOffList > str2double(app.rightFootStrikeList{i}), 1 ) );
                    if ismember(app.rightFootStrikeList{i}, selectedRight)
                        if ~isempty(footOffFrame)
                            fill(app.UIAxesSteps, [str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i}) footOffFrame str2double(app.rightFootStrikeList{i + 1})], [0 -1 -0.8 0], 'g', 'FaceAlpha', 1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        else
                            fill(app.UIAxesSteps, [str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i + 1})], [0 -1 0], 'g', 'FaceAlpha', 1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        end
                    else
                        if ~isempty(footOffFrame)
                            fill(app.UIAxesSteps, [str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i}) footOffFrame str2double(app.rightFootStrikeList{i + 1})], [0 -1 -0.8 0], 'k', 'FaceAlpha', 0.1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        else
                            fill(app.UIAxesSteps, [str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i}) str2double(app.rightFootStrikeList{i + 1})], [0 -1 0], 'k', 'FaceAlpha', 0.1, 'EdgeColor', "none", 'ButtonDownFcn', {@fillCallback, app});
                        end
                    end
                end
                writeFPtoCyclePlot(app);
            end

            % Nested function - indexes app property "states" using input argument "s"
            function fillCallback(src,event, app)
                % disp('clicked')
                src.Vertices(1, 1);
                src.Vertices(4, 1);
                if src.Vertices(2, 2) > 0 % left
                    for i = 1 : size(app.leftListBox.Items, 2)
                        if str2double(app.leftListBox.Items(i)) == src.Vertices(1, 1)
                            if sum(ismember(str2double(app.leftListBox.Value), src.Vertices(1, 1))) == 1  % it was selected, deselect now
                                app.leftListBox.Value(ismember(str2double(app.leftListBox.Value), src.Vertices(1, 1))) = [];
                            else % was not selected, select it now
                                app.leftListBox.Value(end+1) = app.leftListBox.Items(i);
                            end
                        end
                    end
                else % right
                    for i = 1 : size(app.rightListBox.Items, 2)
                        if str2double(app.rightListBox.Items(i)) == src.Vertices(1, 1)
                            if sum(ismember(str2double(app.rightListBox.Value), src.Vertices(1, 1))) == 1  % it was selected, deselect now
                                app.rightListBox.Value(ismember(str2double(app.rightListBox.Value), src.Vertices(1, 1))) = [];
                            else % was not selected, select it now
                                app.rightListBox.Value(end+1) = app.rightListBox.Items(i);
                            end
                        end
                    end
                end

                updateCycles(app, event);
            end
        end

        function writeFPtoCyclePlot(app)
            for i = 1 : size(app.leftStepForcePlateAssignment, 1)
                text(app.UIAxesSteps, app.leftStepForcePlateAssignment(i, 1), 0.5, ['FP' num2str(app.leftStepForcePlateAssignment(i, 2))]);
            end
            for i = 1 : size(app.rightStepForcePlateAssignment, 1)
                text(app.UIAxesSteps, app.rightStepForcePlateAssignment(i, 1), -0.5, ['FP' num2str(app.rightStepForcePlateAssignment(i, 2))]);
            end
        end

        function filterGRFs(app)

            threshold = 20;
            timeThresholdFPs = 10;
            if app.AMTITandemTreadmillCheckBox.Value == 1
                timeThresholdFPs = 5;
                threshold = 5;
            end
            if app.filterGRFsCheckBox.Value && app.c3d.getNumForces() > 0
                app.grf_forces_adjusted = [];
                Fs = app.forceFrequency;  % Sampling Frequency
                N  = app.filterorderSpinner.Value;  % Filter Order
                Fc = app.cutofffrequencySpinner.Value;  % Cutoff Frequency
                [filter_b, filter_a] = butter(N/2,Fc/(Fs/2), 'low'); % divide order by 2 because filtfilt doubles it again
                forcesFields = fieldnames(app.grf_forces);

                forcesFields = sort(forcesFields);
                fpNumbers = [];
                for i = 1 : numel(forcesFields)
                    if char(forcesFields(i)) ~= "time"
                        numbers = regexp(forcesFields{i}, '\d+\.?\d*', 'match'); % Extract numbers as cell array
                        fpNumbers(end+1) = str2double(numbers);
                    end
                end
                fpNumbersUnique = unique(fpNumbers);

                app.grf_forces_adjusted.time = app.grf_forces.time;

                for i = 1 : numel(fpNumbersUnique)
                    fieldsIdxOfThisPlate = fpNumbers==fpNumbersUnique(i);
                    fieldsOfThisPlate = forcesFields(fieldsIdxOfThisPlate);
                    for f = 1 : numel(fieldsOfThisPlate)
                        if strcmp(fieldsOfThisPlate{f}(end-1:end), 'vy')
                            fieldForThreshold = fieldsOfThisPlate{f};
                            continue;
                        end
                    end
                    % disp(fieldForThreshold);

                    Fz = app.grf_forces.(fieldForThreshold);
                    Fz(Fz<threshold) = nan;

                    % logic analog to orginial AMTITreadmillGaitCycleEvents.m from Vicon support website
                    % Find the locations of the NaN values. This will give a value for -1 for the first frame at which the data is not a NaN and a 1 for the last frame for NaN's.
                    nanValues = diff([true; isnan(Fz); true]);
                    footStrikes = find(nanValues < 0);
                    footOffs = find(nanValues > 0)-1;

                    validChanges = find(footOffs - footStrikes > app.frequency_ * 0.5);
                    forceFieldFilterPadding = floor(0.05*app.forceFrequency);

                    footStrikes = footStrikes(validChanges);
                    footOffs = footOffs(validChanges);
                    hasValidSteps = 0;
                    try
                        if Fz(1) > threshold
                            footStrikes = footStrikes(2:end);
                        end

                        if Fz(end) > threshold
                            footOffs = footOffs(1:end-1);
                        end

                        if footOffs(1) < footStrikes(1)
                            footOffs = footOffs(2:end);
                        end
                        hasValidSteps = 1;
                    catch e
                        if strcmp(e.message, 'Index exceeds array bounds.')
                            hasValidSteps = 0;
                        end
                    end

                    if strcmp(app.MovementTypeDropDown.Value, 'Gait') && hasValidSteps

                        footStrikes_Optical = floor(footStrikes / app.forceFrequency * app.frequency_) + 1 + app.firstFrame_;
                        footOffs_Optical = floor(footOffs / app.forceFrequency * app.frequency_) + app.firstFrame_;

                        footStrikes_Optical = footStrikes_Optical;
                        footOffs_Optical = footOffs_Optical;
                        frameThreshold = 5;

                        % check the selected steps on this plate
                        footStrikesOnThisForcePlate = [];
                        for step = 1 : size(app.leftStepForcePlateAssignment, 1)
                            if app.leftStepForcePlateAssignment(step, 2) == fpNumbersUnique(i)
                                footStrikesOnThisForcePlate(end+1) = app.leftStepForcePlateAssignment(step, 1);
                            end
                        end
                        for step = 1 : size(app.rightStepForcePlateAssignment, 1)
                            if app.rightStepForcePlateAssignment(step, 2) == fpNumbersUnique(i)
                                footStrikesOnThisForcePlate(end+1) = app.rightStepForcePlateAssignment(step, 1);
                            end
                        end

                        hasValidStepThatIsAlsoSelected = 0;
                        for f = 1 : numel(fieldsOfThisPlate)
                            newData = zeros(size(app.grf_forces.(fieldsOfThisPlate{f})));
                            % filteredData = filtfilt(filter_b, filter_a, app.grf_forces.(fieldsOfThisPlate{f}));
                            % filter each step seperately if it is selected in the GUI and concatenate results
                            for step = 1 : numel(footStrikes_Optical)
                                if sum(and(footStrikesOnThisForcePlate >= footStrikes_Optical(step) - frameThreshold, footStrikesOnThisForcePlate <= footStrikes_Optical(step) + frameThreshold)) > 0
                                    idxFootContact = footStrikes(step) : footOffs(step);
                                    hasValidStepThatIsAlsoSelected = 1;

                                    if contains(fieldsOfThisPlate{f}(end-4:end), '_p')
                                        % cut several frames at beginning and end to ensure that the point of application does not go back to zero - would result in problems for filtered data
                                        dataToFilter = app.grf_forces.(fieldsOfThisPlate{f})(idxFootContact(timeThresholdFPs+1 : end-timeThresholdFPs));
                                        % filter the data
                                        filteredStep = filtfilt(filter_b, filter_a, dataToFilter);
                                        % interpolate to the original length
                                        lengthOfStep = size(idxFootContact, 2);
                                        filteredStep_originalLength = interp1(1:size(filteredStep, 1), filteredStep, linspace(1, size(filteredStep, 1), lengthOfStep));
                                        % set the filtered data to the output variable
                                        newData(idxFootContact) = filteredStep_originalLength;

                                    elseif contains(fieldsOfThisPlate{f}(end-4:end), '_v') || contains(fieldsOfThisPlate{f}(end-4:end), '_m')
                                        dataToFilter = app.grf_forces.(fieldsOfThisPlate{f})(idxFootContact);
                                        dataToFilter = [zeros(forceFieldFilterPadding, 1); dataToFilter; zeros(forceFieldFilterPadding, 1)]; % add 0.05s of zeros at the start and end to make filtering better
                                        % filter the data
                                        filteredStep = filtfilt(filter_b, filter_a, dataToFilter);
                                        newData(idxFootContact) = filteredStep(1+forceFieldFilterPadding : end-forceFieldFilterPadding);
                                    else
                                        disp('not sure what field this is');
                                    end
                                end
                            end

                            if hasValidStepThatIsAlsoSelected
                                app.grf_forces_adjusted.(fieldsOfThisPlate{f}) = newData;
                            else
                                % no selected steps on this force plate, filter where Fz is over threshold
                                idxFootContact = ~isnan(Fz);
                                newData = zeros(size(app.grf_forces.(fieldsOfThisPlate{f})));
                                filteredData = filtfilt(filter_b, filter_a, app.grf_forces.(fieldsOfThisPlate{f})(idxFootContact));
                                newData(idxFootContact) = filteredData;
                                app.grf_forces_adjusted.(fieldsOfThisPlate{f}) = newData;
                            end
                        end
                    else
                        % no valid steps on this force plate, filter where Fz is over threshold
                        idxFootContact = ~isnan(Fz);

                        for f = 1 : numel(fieldsOfThisPlate)
                            newData = zeros(size(app.grf_forces.(fieldsOfThisPlate{f})));
                            try
                                filteredData = filtfilt(filter_b, filter_a, app.grf_forces.(fieldsOfThisPlate{f})(idxFootContact));
                                newData(idxFootContact) = filteredData;
                            end
                            app.grf_forces_adjusted.(fieldsOfThisPlate{f}) = newData;
                        end
                    end
                end
            else
                app.grf_forces_adjusted = app.grf_forces;
            end
        end

        function resetGUI(app)
            cla(app.UIAxesSteps);
            cla(app.UIAxesGRFs);
        end


        function readEMG(app)
            app.EMG = struct;
            if ~isempty(app.acq)
                [analogs, analogsInfo] = btkGetAnalogs(app.acq);
                ratioEmgToFrames = btkGetAnalogSampleNumberPerFrame(app.acq);
                app.EMG.frequency = btkGetAnalogFrequency(app.acq);
                analogChannels = fieldnames(analogsInfo.units);
                for i = 1 : numel(analogChannels)
                    if contains(analogsInfo.description.(analogChannels{i}), 'EMG')
                        app.EMG.(analogChannels{i}) = analogs.(analogChannels{i});
                    end
                end

                emgFields = fieldnames(app.EMG);
                emgFields = emgFields(~contains(emgFields, 'frequency'));
                if size(emgFields, 1) ~= 0
                    try
                        app.EMG.time(1) = app.markers.time(1);
                    catch
                        app.EMG.time(1) = 0;
                    end
                    for j = 2 : size(app.EMG.(emgFields{1}), 1)
                        app.EMG.time(j) = app.EMG.time(j-1) + 1/app.EMG.frequency;
                    end
                    app.EMG.time = app.EMG.time';
                else
                    app.filterEMGCheckBox.Value = 0;
                    app.renameEMGlabelsCheckBox.Value = 0;
                end
            end
        end

        function filterEMG(app)
            if ~isempty(app.EMG)
                app.EMG_lowPass = struct;

                Fs = app.EMG.frequency;  % Sampling Frequency
                band    = (2 / Fs) * [20, 400];
                [bandpass_b, bandpass_a] = butter(1, band, 'bandpass');
                [butter4_b, butter4_a] = butter(2, 10/(Fs/2), 'low'); % divide order by 2 because filtfilt doubles it again

                channels = fieldnames(app.EMG);
                channels = channels(~contains(channels, 'frequency'));
                channels = channels(~contains(channels, 'time'));
                for i = 1 : numel(channels)
                    tmp_after_bandpass = filter(bandpass_b, bandpass_a, app.EMG.(channels{i}), [], 1);
                    EMG_demeaned            = tmp_after_bandpass - mean(tmp_after_bandpass); %demeaned
                    EMG_rectified           = sqrt(EMG_demeaned.^2); %full wave rectified
                    app.EMG_lowPass.(channels{i}) = filtfilt(butter4_b, butter4_a, EMG_rectified);

                    app.EMG_lowPass.(channels{i}) = movmean(app.EMG_lowPass.(channels{i}), Fs*0.05);
                    app.EMG_lowPass.(channels{i})(app.EMG_lowPass.(channels{i})<0) = 0; %values under 0 (due to the low pass filter) are set to 0)
                    if app.normalizeallmagnitudesto1CheckBox.Value == 1
                        app.EMG_lowPass.(channels{i}) = app.EMG_lowPass.(channels{i}) / max(app.EMG_lowPass.(channels{i}));
                    else
                        app.EMG_lowPass.(channels{i}) = app.EMG_lowPass.(channels{i});
                    end
                end
            end
        end

        function processC3Dfile(app)

            if app.c3dFileName ~= 0
                d = uiprogressdlg(app.PrepareC3DfileUIFigure,'Title','Please Wait', 'Message','Opening the selected C3D file and performing some calculations ... ');
                app.SelectaC3DfileLabel.FontColor = 'black';
                app.SelectaC3DfileLabel.Text = fullfile(app.path, app.c3dFileName);

                app.openselectedfilewithdefaultprogramButton.Visible = 1;

                % get events and put them in the list for selection
                app.acq = btkReadAcquisition(fullfile(app.path, app.c3dFileName));
                c3devents = btkGetEvents(app.acq);
                app.frequency_ = btkGetPointFrequency(app.acq);
                app.firstFrame_ = btkGetFirstFrame(app.acq);
                d.Value = d.Value + 0.1;

                if app.removeblanksetcfromfilenameCheckBox.Value == 1
                    validFileName = strrep(strrep(app.c3dFileName, '.c3d', ''), ' ', '_');
                    validFileName = matlab.lang.makeValidName(validFileName);

                    enfFiles = dir(fullfile(app.path, strrep(app.c3dFileName, '.c3d', '.*enf')));


                    if ~strcmp(validFileName, strrep(app.c3dFileName, '.c3d', ''))
                        copyfile(fullfile(app.path, app.c3dFileName), fullfile(app.path, [validFileName '.c3d']));
                        if app.deleteoriginalfileafterrenamingCheckBox.Value == 1
                            delete(fullfile(app.path, app.c3dFileName));
                        end
                        app.c3dFileName = [validFileName '.c3d'];

                        if numel(enfFiles) == 1
                            copyfile(fullfile(enfFiles(1).folder, enfFiles(1).name), fullfile(app.path, [validFileName '.enf']));
                            if app.deleteoriginalfileafterrenamingCheckBox.Value == 1
                                delete(fullfile(enfFiles(1).folder, enfFiles(1).name));
                            end
                        end
                        
                        app.processC3Dfile();
                        return;
                    end
                end

                if app.Ignorec3deventsCheckBox.Value == 0
                    app.validfootstrikeeventsLabel.Visible = "on";
                    app.forceplatecontactsLabel.Visible = "off";
                    if isfield(c3devents, 'Left_Foot_Strike')
                        app.leftFootStrikeList = sprintfc('%.0f', c3devents.Left_Foot_Strike * app.frequency_);
                        app.leftListBox.Items = app.leftFootStrikeList(1 : end-1);
                        if strcmp(app.MovementTypeDropDown.Value, 'Stand-Sit-Stand') || strcmp(app.MovementTypeDropDown.Value, 'Counter Movement Jump')
                            app.leftListBox.Value = app.leftFootStrikeList([1 3]);
                        else
                            app.leftListBox.Value = app.leftFootStrikeList(1 : end-1);
                        end
                    else
                        app.leftListBox.Items = {};
                        app.leftFootStrikeList = {};
                    end
                    if isfield(c3devents, 'Left_Foot_Off')
                        app.leftFootOffList = floor(c3devents.Left_Foot_Off * app.frequency_);
                    else
                        app.leftFootOffList = [];
                    end
                    if isfield(c3devents, 'Right_Foot_Strike')
                        app.rightFootStrikeList = sprintfc('%.0f', c3devents.Right_Foot_Strike * app.frequency_);
                        app.rightListBox.Items = app.rightFootStrikeList(1 : end-1);
                        if strcmp(app.MovementTypeDropDown.Value, 'Stand-Sit-Stand') || strcmp(app.MovementTypeDropDown.Value, 'Counter Movement Jump')
                            app.rightListBox.Value = app.rightFootStrikeList([1 3]);
                        else
                            app.rightListBox.Value = app.rightFootStrikeList(1 : end-1);
                        end
                    else
                        app.rightListBox.Items = {};
                        app.rightFootStrikeList = {};
                    end
                    if isfield(c3devents, 'Right_Foot_Off')
                        app.rightFootOffList = floor(c3devents.Right_Foot_Off * app.frequency_);
                    else
                        app.rightFootOffList = [];
                    end
                end

                % drawnow;
                % figure(app.PrepareC3DfileUIFigure);

                offsetFrames = 10;
                if isfield(c3devents, 'Right_Foot_Strike') && isfield(c3devents, 'Left_Foot_Strike')
                    xminlimit = min(c3devents.Right_Foot_Strike(1), c3devents.Left_Foot_Strike(1)) * app.frequency_ - offsetFrames;
                    xmaxlimit = max(c3devents.Right_Foot_Strike(end), c3devents.Left_Foot_Strike(end)) * app.frequency_ + offsetFrames;
                elseif isfield(c3devents, 'Right_Foot_Strike')
                    xminlimit = c3devents.Right_Foot_Strike(1) * app.frequency_ - offsetFrames;
                    xmaxlimit = c3devents.Right_Foot_Strike(end) * app.frequency_ - offsetFrames;
                elseif isfield(c3devents, 'Left_Foot_Strike')
                    xminlimit = c3devents.Left_Foot_Strike(1) * app.frequency_ - offsetFrames;
                    xmaxlimit = c3devents.Left_Foot_Strike(end) * app.frequency_ - offsetFrames;
                else
                    xminlimit = 0;
                    xmaxlimit = 5000;
                end
                if xminlimit == xmaxlimit
                    xmaxlimit = xminlimit + 5000;
                end
                xlim(app.UIAxesSteps, [xminlimit xmaxlimit]);
                ylim(app.UIAxesSteps, [-1.1 1.1]);

                c3dpath = fullfile(app.path, app.c3dFileName);
                d.Value = d.Value + 0.1;

                try
                    app.c3d = osimC3D(c3dpath, 1);
                catch e
                    if contains(e.message, 'Number of labels does not match number of columns of dependent data.')

                        markersToRemove = [
                            "glut_med","glut_min","semimem","semiten","bifemlh","bifemsh","sar", ...
                            "add_long","add_brev","add_mag","tfl","pect","grac","glut_max", ...
                            "iliacus","psoas","quad_fem","gem","peri", ...
                            "rect_fem","vas_med","vas_int","vas_lat", ...
                            "med_gas","lat_gas","soleus", ...
                            "tib_post","flex_dig","flex_hal","tib_ant", ...
                            "per_brev","per_long","per_tert", ...
                            "ext_dig","ext_hal", '_force', '_length', '_l_norm', '_l', '_r', '_r_norm', ...
                            'usermo', 'pelvis', 'femur', 'tibia', 'foot', 'thorax', 'upperarm', 'forearm', 'hand', 'head', ...
                            'medial_knee', 'patfem', 'walker_knee', 'ankle_', 'hip_'
                            ];

                        reducedPath = strrep(c3dpath, '.c3d', '_reduced.c3d');

                        maxAttempts = 10;
                        for attempt = 1:maxAttempts
                            if attempt == 1
                                tmpAcq = btkReadAcquisition(c3dpath);
                            else
                                tmpAcq = btkReadAcquisition(reducedPath);
                            end

                            md = btkGetMetaData(tmpAcq,'POINT','LABELS');
                            pointLabels = strtrim(string(md.info.values));
                            for k = 1:numel(pointLabels)
                                label = pointLabels{k};
                                if contains(lower(label), lower(markersToRemove)) || startsWith(label, 'C_')
                                    btkRemovePoint(tmpAcq, label);
                                end
                            end
                            btkWriteAcquisition(tmpAcq, reducedPath);
                            btkCloseAcquisition(tmpAcq);

                            try
                                app.c3d = osimC3D(reducedPath, 1);
                                break; % success → exit loop
                            catch e
                                if attempt == maxAttempts || ...
                                        ~contains(e.message, 'Number of labels does not match number of columns of dependent data.')
                                    rethrow(e); % unexpected error or out of retries
                                end
                            end
                        end
                        app.writeLog('Some markers/model outputs were removed to open file ');
                    else
                        rethrow(e);
                    end
                end
                try
                    app.forceFrequency = app.c3d.getRate_force();
                catch ME
                    app.forceFrequency = 1000;
                end
                app.c3d.rotateData('x', -90);
                try
                    app.c3d.convertMillimeters2Meters();
                catch
                    app.writeLog('probably no markers in the file!');
                end
                app.markers = osimTableToStruct(app.c3d.getTable_markers);

                d.Value = d.Value + 0.1;
                d.Message = 'writing forces to temporary file...';

                if app.c3d.getNumForces() ~= 0
                    tmp_mot_file = fullfile(getenv('TEMP'), 'temp.mot');
                    app.c3d.writeMOT(tmp_mot_file);
                    %replace nan and -nan(ind)
                    fid  = fopen(tmp_mot_file,'r');
                    f=fread(fid,'*char')';
                    fclose(fid);
                    f = strrep(f,'-nan(ind)','0');
                    f = strrep(f,'nan','0');
                    fid  = fopen(tmp_mot_file,'w');
                    fprintf(fid,'%s',f);
                    fclose(fid);
                    app.grf_forces = load_sto_file(tmp_mot_file);
                    delete(tmp_mot_file);
                else
                    app.grf_forces = struct;
                end

                if app.Ignorec3deventsCheckBox.Value == 1
                    grfForceNames = fieldnames(app.grf_forces);
                    cnt = 1;
                    if numel(grfForceNames) > 0
                        for i = 1 : numel(grfForceNames)
                            if contains(grfForceNames{i}, '_vy')
                                disp(grfForceNames{i}(end-9 : end-3));
                                leftList{cnt} = grfForceNames{i}(end-9 : end-3);
                                rightList{cnt} = grfForceNames{i}(end-9 : end-3);
                                cnt = cnt + 1;
                            end
                        end
                    else
                        leftList = {};
                        rightList = {};
                    end

                    if ~isequal(leftList, app.leftListBox) || ~isequal(rightList, app.rightListBox)
                        app.rightListBox.Items = rightList;
                        app.leftListBox.Items = leftList;
                    end
                end

                app.forceOffset = app.firstFrame_ / app.frequency_ * app.forceFrequency;

                if isfield(c3devents, 'Right_Foot_Strike') && isfield(c3devents, 'Left_Foot_Strike')
                    xminlimit = min(c3devents.Right_Foot_Strike(1), c3devents.Left_Foot_Strike(1)) * app.forceFrequency - offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                    xmaxlimit = max(c3devents.Right_Foot_Strike(end), c3devents.Left_Foot_Strike(end)) * app.forceFrequency + offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                elseif isfield(c3devents, 'Right_Foot_Strike')
                    xminlimit = c3devents.Right_Foot_Strike(1) * app.forceFrequency - offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                    xmaxlimit = c3devents.Right_Foot_Strike(end) * app.forceFrequency + offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                elseif isfield(c3devents, 'Left_Foot_Strike')
                    xminlimit = c3devents.Left_Foot_Strike(1) * app.forceFrequency - offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                    xmaxlimit = c3devents.Left_Foot_Strike(end) * app.forceFrequency + offsetFrames / app.frequency_ * app.forceFrequency - app.forceOffset;
                else
                    xminlimit = 0;
                    xmaxlimit = 10000;
                end
                if (xminlimit - xmaxlimit <= 0)
                    xmaxlimit = xminlimit + 1;
                end

                xlim(app.UIAxesGRFs, [xminlimit xmaxlimit]);

                app.fileLoaded = 1;

                d.Value = d.Value + 0.1;
                d.Message = 'detect walking direction file...';
                detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged(app, 0);
                d.Value = d.Value + 0.1;
                d.Message = 'detect force plates for each step...';
                if app.c3d.getNumForces() ~= 0
                    detectforceplatesautomaticallyCheckBoxValueChanged(app, 0);
                end
                d.Value = d.Value + 0.1;
                d.Message = 'updating plots...';
                updateGRF_Plot(app, 0);
                updateCyclePlot(app);
                d.Value = 1;
                d.Message = 'Finished';
                writeLog(app, [app.c3dFileName ' sucessfully read']);
                app.CreatefilesButton.Enable = "on";
                close(d);

                if app.filterEMGCheckBox.Value
                    readEMG(app);
                    filterEMG(app);
                end
            else
                resetGUI(app);
                writeLog(app, 'Select a file! - data which may be displayed is wrong - Select a file first!');
                app.SelectaC3DfileLabel.Text = 'Select a file! - data which may be displayed is wrong - Select a file first!';
                app.SelectaC3DfileLabel.FontColor = 'red';
                app.CreatefilesButton.Enable = "off";
                app.fileLoaded = 0;
                app.openselectedfilewithdefaultprogramButton.Visible = 0;
            end
        end

        function mirrorMarkerData(app, filename)
            disp('mirroring marker data...')
            d.Message = 'mirroring marker data...';

            text = fileread(filename);
            lines = strsplit(text, '\n');
            data = [];
            for i = 1 : numel(lines)
                data{i} = strsplit(lines{i}, '\t', 'CollapseDelimiters', false);
            end

            for col = 3 : numel(data{7}) % skip time and Frame
                traj = [];
                for row = 7 : numel(data)-1 % skip 6 header rows
                    traj(row-6) = str2double(data{row}{col});
                end
                if sum(isnan(traj)) == 0
                    if contains(data{5}{col}, 'Z')
                        traj = traj * - 1;
                    end
                    % filteredTraj = filtfilt(filter_b, filter_a, traj);
                    for row = 7 : numel(data)-1 % skip 6 header rows
                        data{row}{col} = traj(row-6);
                        % traj(row-6) = cell2mat();
                    end
                end
                if ~isempty(data{4}{col}) && strcmp(data{4}{col}(1), 'L')
                    data{4}{col}(1) = 'R';
                elseif ~isempty(data{4}{col}) && strcmp(data{4}{col}(1), 'R')
                    data{4}{col}(1) = 'L';
                end
            end


            data{1}{4} = strrep(data{1}{4}, '\', '\\');

            finalStr = '';
            for row = 1 : numel(data)-1
                rowStr = '';
                for col = 1 : numel(data{row})
                    try rowStr = [rowStr, num2str(data{row}{col}), '\t'];
                    catch rowStr = [rowStr, data{row}{col}, '\t'];
                    end
                end
                if ~isempty(rowStr)
                    rowStr(end-2 : end) = [];
                    rowStr = [rowStr, '\n'];
                else
                    rowStr = '\n';
                end

                finalStr = [finalStr, rowStr];
            end


            fid = fopen(filename,'wt');
            fprintf(fid, finalStr);
            fclose(fid);
        end

        function checkMovementType(app)
            if app.AMTITandemTreadmillCheckBox.Value == 0
                if contains(app.c3dFileName, 'cmj', 'IgnoreCase', true)
                    app.MovementTypeDropDown.Value = 'Counter Movement Jump';
                    app.exporteachstepseparatelyCheckBox.Value = 1;
                    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
                    app.rotatearoundyaxisDropDown.Value = '180°';
                elseif contains(app.c3dFileName, 'sts', 'IgnoreCase', true)
                    app.MovementTypeDropDown.Value = 'Stand-Sit-Stand';
                    app.exporteachstepseparatelyCheckBox.Value = 1;
                    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
                    app.rotatearoundyaxisDropDown.Value = '180°';
                elseif contains(app.c3dFileName, 'sidestep', 'IgnoreCase', true)
                    app.MovementTypeDropDown.Value = 'SideStep';
                    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
                    app.rotatearoundyaxisDropDown.Value = '180°';
                    app.exporteachstepseparatelyCheckBox.Value = 0;
                else
                    app.MovementTypeDropDown.Value = 'Gait';
                    app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 1;
                    app.exporteachstepseparatelyCheckBox.Value = 1;
                end
            end
        end
    end


    % Callbacks that handle component events
    methods (Access = public)

        % Code that executes after component creation
        function startupFcn(app)
            %             movegui(app.PrepareC3DfileUIFigure, "south");
            % add OpenSim jar file to java path
            addingSuccessfull = 0;
            systemPaths = strsplit(getenv('PATH'), ';');

            for i = 1 : numel(systemPaths)
                if contains(systemPaths{i}, 'OpenSim')
                    systemPaths{i} = strrep(systemPaths{i}, '"', '');
                    pathToAdd = systemPaths{i};
                    javaaddpath(pathToAdd);
                    addingSuccessfull = 1;
                    break;
                end
            end

            if ~addingSuccessfull
                errordlg('Check if you have installed OpenSim properly and added it to System PATH variable. Then restart application! Typically the installation path is C:\OpenSim 4.2\bin');
                %                 app.delete;
            end

            hold(app.UIAxesSteps, 'on');
            hold(app.UIAxesGRFs, 'on');
        end

        % Button pushed function: Selectc3dfileButton
        function Selectc3dfileButtonPushed(app, event)

            if ~isempty(app.path) && ischar(app.path)
                [app.c3dFileName, app.path] = uigetfile('*.c3d', 'MultiSelect', 'off', 'Select C3D file of trial', app.path);
            else
                [app.c3dFileName, app.path] = uigetfile('*.c3d', 'MultiSelect', 'off', 'Select C3D file of trial');
            end

            if app.automaticMovementtypeCheckBox.Value
                app.checkMovementType();
            end

            % drawnow;
            % figure(app.PrepareC3DfileUIFigure);

            app.processC3Dfile();
        end

        % Callback function: detectWalkingByMarker, detectWalkingByMarker, 
        % ...and 1 other component
        function detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged(app, event)
            value = app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value;
            app.rotatearoundyaxisDropDown.Enable = mod(value + 1, 2);
            if app.fileLoaded
                if value % detect walking direction and
                    if isfield(app.markers, app.detectWalkingByMarker.Value)
                        app.detectWalkingByMarker.BackgroundColor = 'white';
                        app.detectWalkingByMarker.FontColor = 'black';
                        markerToTrack = app.detectWalkingByMarker.Value;

                        firstValidLocation = nan(1, 3);
                        for i = 1 : size(app.markers.(markerToTrack), 1)
                            if ~isnan(app.markers.(markerToTrack)(i, :))
                                firstValidLocation = app.markers.(markerToTrack)(i, :);
                                break;
                            end
                        end
                        lastValidLocation = nan(1, 3);
                        for i = size(app.markers.(markerToTrack), 1) : -1 : 1
                            if ~isnan(app.markers.(markerToTrack)(i, :))
                                lastValidLocation = app.markers.(markerToTrack)(i, :);
                                break;
                            end
                        end
                        distanceTravelled = lastValidLocation - firstValidLocation;
                        [~, maxInd] = max(abs(distanceTravelled));
                        if maxInd == 3
                            if distanceTravelled(maxInd) < 0
                                app.rotatearoundyaxisDropDown.Value = '270°';
                                app.rotationAngle = 270;
                            else
                                app.rotatearoundyaxisDropDown.Value = '90°';
                                app.rotationAngle = 90;
                            end
                        elseif maxInd == 1
                            if distanceTravelled(maxInd) < 0
                                app.rotatearoundyaxisDropDown.Value = '180°';
                                app.rotationAngle = 180;
                            else
                                app.rotatearoundyaxisDropDown.Value = '0°';
                                app.rotationAngle = 0;
                            end
                        end

                    else
                        app.detectWalkingByMarker.Value = 'SET THIS!';
                        app.detectWalkingByMarker.BackgroundColor = 'red';
                        app.detectWalkingByMarker.FontColor = 'white';
                        app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
                    end
                end
            end
        end

        % Button pushed function: CreatefilesButton
        function CreatefilesButtonPushed(app, event)

            if app.fileLoaded
                d = uiprogressdlg(app.PrepareC3DfileUIFigure,'Title','Please Wait', 'Message','Creating files with selected options ... ');

                [folder, baseFileNameNoExt, ~] = fileparts(fullfile(app.path, app.c3dFileName));

                output_folder = fullfile(folder, baseFileNameNoExt);

                if app.mirrorCheckBox.Value
                    output_folder = [output_folder '_mirrored'];
                end



                d.Value = d.Value + 0.1;
                d.Message = 'rotating data ...';
                if ~app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value
                    app.rotationAngle = str2double(strrep(app.rotatearoundyaxisDropDown.Value, '°', ''));
                end
                c3d_tmp = app.c3d;
                c3d_tmp.rotateData('y', app.rotationAngle);

                outputSets = struct;
                if app.exporteachstepseparatelyCheckBox.Value == 1

                    for i = 1 : size(app.rightListBox.Value, 2)
                        outputFolderName = [baseFileNameNoExt '_R' num2str(i)];
                        if ~isempty(app.rightStepForcePlateAssignment)
                            idxAssignment = find(app.rightStepForcePlateAssignment(:, 1) == str2double(app.rightListBox.Value(i)));
                            if ~isempty(idxAssignment)
                                if strcmp(app.MovementTypeDropDown.Value, 'Stand-Sit-Stand')
                                    if i == 1
                                        outputFolderName = [baseFileNameNoExt '_SitDown'];
                                    else
                                        outputFolderName = [baseFileNameNoExt '_StandUp'];
                                    end
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(idxAssignment, :);
                                elseif strcmp(app.MovementTypeDropDown.Value, 'Counter Movement Jump')
                                    if i == 1
                                        outputFolderName = [baseFileNameNoExt '_TakeOff'];
                                    else
                                        outputFolderName = [baseFileNameNoExt '_Landing'];
                                    end
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(idxAssignment, :);
                                else
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = [];
                                end
                            else
                                outputSets.(outputFolderName).rightStepForcePlateAssignment = [];
                                outputSets.(outputFolderName).leftStepForcePlateAssignment = [];
                            end
                        else
                            outputSets.(outputFolderName).rightStepForcePlateAssignment = [];
                            outputSets.(outputFolderName).leftStepForcePlateAssignment = [];
                        end
                        outputSets.(outputFolderName).firstIC = str2double(app.rightListBox.Value(i));
                        outputSets.(outputFolderName).lastIC = outputSets.(outputFolderName).firstIC; % only this one step
                        outputSets.(outputFolderName).rightICs = str2double(app.rightListBox.Value(i));
                        outputSets.(outputFolderName).leftICs = [];
                    end

                    for i = 1 : size(app.leftListBox.Value, 2)
                        outputFolderName = [baseFileNameNoExt '_L' num2str(i)];
                        if ~isempty(app.leftStepForcePlateAssignment)
                            idxAssignment = find(app.leftStepForcePlateAssignment(:, 1) == str2double(app.leftListBox.Value(i)));
                            if ~isempty(idxAssignment)
                                if strcmp(app.MovementTypeDropDown.Value, 'Stand-Sit-Stand')
                                    if i == 1
                                        outputFolderName = [baseFileNameNoExt '_SitDown'];
                                    else
                                        outputFolderName = [baseFileNameNoExt '_StandUp'];
                                    end
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(idxAssignment, :);
                                elseif strcmp(app.MovementTypeDropDown.Value, 'Counter Movement Jump')
                                    if i == 1
                                        outputFolderName = [baseFileNameNoExt '_TakeOff'];
                                    else
                                        outputFolderName = [baseFileNameNoExt '_Landing'];
                                    end
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(idxAssignment, :);
                                else
                                    outputSets.(outputFolderName).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(idxAssignment, :);
                                    outputSets.(outputFolderName).rightStepForcePlateAssignment = [];
                                end
                            else
                                outputSets.(outputFolderName).rightStepForcePlateAssignment = [];
                                outputSets.(outputFolderName).leftStepForcePlateAssignment = [];
                            end
                        else
                            outputSets.(outputFolderName).rightStepForcePlateAssignment = [];
                            outputSets.(outputFolderName).leftStepForcePlateAssignment = [];
                        end
                        outputSets.(outputFolderName).firstIC = str2double(app.leftListBox.Value(i));
                        outputSets.(outputFolderName).lastIC = outputSets.(outputFolderName).firstIC; % only this one step
                        outputSets.(outputFolderName).leftICs = str2double(app.leftListBox.Value(i));
                        outputSets.(outputFolderName).rightICs = [];
                    end
                else
                    outputSets.(baseFileNameNoExt).rightStepForcePlateAssignment = app.rightStepForcePlateAssignment;
                    outputSets.(baseFileNameNoExt).leftStepForcePlateAssignment = app.leftStepForcePlateAssignment;
                    outputSets.(baseFileNameNoExt).firstIC = min([str2double(app.leftListBox.Value) str2double(app.rightListBox.Value)]);
                    outputSets.(baseFileNameNoExt).lastIC = max([str2double(app.leftListBox.Value) str2double(app.rightListBox.Value)]);
                    outputSets.(baseFileNameNoExt).leftICs = str2double(app.leftListBox.Value);
                    outputSets.(baseFileNameNoExt).rightICs = str2double(app.rightListBox.Value);
                end

                cropCheckBoxOriginal = app.croptovalidfootstrikesbufferCheckBox.Value;
                outputFolders = fieldnames(outputSets);
                if isempty(outputFolders)
                    outputSets.(baseFileNameNoExt).rightStepForcePlateAssignment = [];
                    outputSets.(baseFileNameNoExt).leftStepForcePlateAssignment = [];
                    outputFolders = fieldnames(outputSets);
                    app.croptovalidfootstrikesbufferCheckBox.Value = 0;
                end
                for o = 1 : numel(outputFolders)
                    outputFolderName = outputFolders{o};
                    output_folder = fullfile(folder, outputFolderName);

                    % create folder with the same name as c3d file
                    if ~exist(output_folder, 'dir')
                        mkdir(output_folder)
                    end

                    try
                        c3d_tmp.writeTRC(fullfile(output_folder, 'marker_experimental.trc'));
                    catch
                        app.writeLog('Markers not exportet - probably none in file!');
                    end
                    if app.croptovalidfootstrikesbufferCheckBox.Value == 1 || app.filtermarkersCheckBox.Value == 1
                        markerData = load_marker_trc(fullfile(output_folder, 'marker_experimental.trc'));
                    end

                    exportStartTime = [];
                    exportEndTime = [];

                    if app.croptovalidfootstrikesbufferCheckBox.Value == 1
                        % identify timing of this output set
                        endFrame = [];
                        firstIC = outputSets.(outputFolderName).firstIC;
                        lastIC = outputSets.(outputFolderName).lastIC;

                        for i = 1 : numel(app.leftFootStrikeList)-1
                            if str2double(app.leftFootStrikeList{i}) == lastIC
                                endFrame = str2double(app.leftFootStrikeList{i+1});
                            end
                        end
                        for i = 1 : numel(app.rightFootStrikeList)-1
                            if str2double(app.rightFootStrikeList{i}) == lastIC
                                endFrame = str2double(app.rightFootStrikeList{i+1});
                            end
                        end
                        firstFrame = app.firstFrame_;
                        frequency = app.frequency_;
                        markerDataTime = cell2mat(markerData.Time);

                        timeOfIC = (firstIC - (firstFrame-1)) / frequency;
                        timeOfEnd = (endFrame - (firstFrame-1)) / frequency;
                        timeBufferStart = app.BufferStartsSpinner.Value;
                        timeBufferEnd = app.BufferEndsSpinner.Value;
                        exportStartTime = max([markerDataTime(1), timeOfIC - timeBufferStart]);
                        exportEndTime = min([markerDataTime(end) timeOfEnd + timeBufferEnd]);

                        % crop marker data
                        keepIdx = and(markerDataTime >= exportStartTime, markerDataTime <= exportEndTime);
                        markerNames = fieldnames(markerData);
                        for i = 1 : numel(markerNames)
                            markerData.(char(markerNames(i))) = markerData.(char(markerNames(i)))(keepIdx);
                        end
                        markerDataTime = markerDataTime(keepIdx);
                        writeTRCFile(fullfile(output_folder, 'marker_experimental.trc'), markerData, markerDataTime, 1);
                    end

                    d.Value = d.Value + 0.1 / numel(outputFolders);
                    d.Message = 'writing marker_experimental.trc file with marker locations...';
                    % c3d_tmp.writeTRC(fullfile(output_folder, 'marker_experimental.trc'));

                    if app.filtermarkersCheckBox.Value
                        disp('filtering marker data...')
                        d.Message = 'filtering marker data...';
                        % markerData = load_marker_trc(fullfile(output_folder, 'marker_experimental.trc'));
                        markerNames = fieldnames(markerData);

                        Fs = app.frequency_;  % Sampling Frequency
                        N  = app.filterorderSpinner_Marker.Value;  % Filter Order
                        Fc = app.cutofffrequencySpinner_Marker.Value;  % Cutoff Frequency
                        [filter_b, filter_a] = butter(N/2,Fc/(Fs/2), 'low'); % divide order by 2 because filtfilt doubles it again

                        for i = 1 : numel(markerNames)
                            if strcmp(markerNames{i}(end-1:end), '_X') || strcmp(markerNames{i}(end-1:end), '_Y') || strcmp(markerNames{i}(end-1:end), '_Z')
                                try
                                    markerData.(char(markerNames(i))) = filtfilt(filter_b, filter_a, cell2mat(markerData.(char(markerNames(i)))));
                                catch e                                    
                                    if contains(e.message, 'Unable to concatenate')
                                        tmp = markerData.(char(markerNames(i)));
                                        numData = nan(size(tmp));
                                        for u = 1 : size(tmp, 1)
                                            try
                                                numData(u) = cell2mat(tmp(u));
                                            end
                                        end
                                        numData = fillmissing(numData, 'spline');
                                        markerData.(char(markerNames(i))) = filtfilt(filter_b, filter_a, numData);
                                    else
                                        markerData = rmfield(markerData, (char(markerNames(i))));
                                    end
                                end
                            else
                                try
                                    markerData.(char(markerNames(i))) = cell2mat(markerData.(char(markerNames(i))));
                                catch
                                    markerData = rmfield(markerData, (char(markerNames(i))));
                                end
                            end
                        end

                        markerDataTime = markerData.Time;
                        status = movefile(fullfile(output_folder, 'marker_experimental.trc'), fullfile(output_folder, 'marker_experimental_unfiltered.trc'));
                        if status
                            writeTRCFile(fullfile(output_folder, 'marker_experimental.trc'), markerData, markerDataTime, 1);
                        else
                            writeTRCFile(fullfile(output_folder, 'marker_experimental_filtered.trc'), markerData, markerDataTime, 1);
                        end
                    end

                    if app.mirrorCheckBox.Value
                        mirrorMarkerData(app, fullfile(output_folder, 'marker_experimental.trc'));
                    end

                    nForces = app.c3d.getNumForces();

                    d.Value = d.Value + 0.1 / numel(outputFolders);
                    d.Message = 'processing forces ...';
                    if (nForces ~= 0)
                        grfFileName = fullfile(output_folder, 'grf.mot');
                        app.c3d.writeMOT(grfFileName);

                        %replace nan and -nan(ind)
                        fid  = fopen(grfFileName,'r');
                        f=fread(fid,'*char')';
                        fclose(fid);
                        f = strrep(f,'-nan(ind)','0');
                        f = strrep(f,'nan','0');
                        fid  = fopen(grfFileName,'w');
                        fprintf(fid,'%s',f);
                        fclose(fid);

                        app.grf_forces = load_sto_file(grfFileName);

                        filterGRFs(app);

                        % reorder fieldnames otherwise OpenSim cannot process it
                        temp_forceStruct.time = app.grf_forces_adjusted.time;
                        app.grf_forces_adjusted = orderfields(app.grf_forces_adjusted);
                        allFields = fieldnames(app.grf_forces_adjusted);
                        for i = 1 : (numel(allFields) - 1) / 9
                            fpNr = allFields{i * 6, 1}(14);
                            temp_forceStruct.(['ground_force_' num2str(fpNr) '_vx']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_vx']);
                            temp_forceStruct.(['ground_force_' num2str(fpNr) '_vy']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_vy']);

                            if app.mirrorCheckBox.Value
                                temp_forceStruct.(['ground_force_' num2str(fpNr) '_vz']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_vz']) * -1;
                            else
                                temp_forceStruct.(['ground_force_' num2str(fpNr) '_vz']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_vz']);
                            end
                            temp_forceStruct.(['ground_force_' num2str(fpNr) '_px']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_px']);
                            temp_forceStruct.(['ground_force_' num2str(fpNr) '_py']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_py']);
                            if app.mirrorCheckBox.Value
                                temp_forceStruct.(['ground_force_' num2str(fpNr) '_pz']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_pz']) * -1;
                            else
                                temp_forceStruct.(['ground_force_' num2str(fpNr) '_pz']) = app.grf_forces_adjusted.(['ground_force_' num2str(fpNr) '_pz']);
                            end
                            temp_forceStruct.(['ground_moment_' num2str(fpNr) '_mx']) = app.grf_forces_adjusted.(['ground_moment_' num2str(fpNr) '_mx']);
                            temp_forceStruct.(['ground_moment_' num2str(fpNr) '_my']) = app.grf_forces_adjusted.(['ground_moment_' num2str(fpNr) '_my']);

                            if app.mirrorCheckBox.Value
                                temp_forceStruct.(['ground_moment_' num2str(fpNr) '_mz']) = app.grf_forces_adjusted.(['ground_moment_' num2str(fpNr) '_mz']) * -1;
                            else
                                temp_forceStruct.(['ground_moment_' num2str(fpNr) '_mz']) = app.grf_forces_adjusted.(['ground_moment_' num2str(fpNr) '_mz']);
                            end
                        end
                        app.grf_forces_adjusted = temp_forceStruct;


                        if app.croptovalidfootstrikesbufferCheckBox.Value == 1
                            fields = fieldnames(app.grf_forces_adjusted);
                            keepIdx = and(app.grf_forces_adjusted.time >= exportStartTime, app.grf_forces_adjusted.time <= exportEndTime);
                            for i = 1 : numel(fields)
                                app.grf_forces_adjusted.(fields{i}) = app.grf_forces_adjusted.(fields{i})(keepIdx);
                            end
                        end

                        d.Value = d.Value + 0.1 / numel(outputFolders);
                        d.Message = 'writing forces to grf.mot file...';

                        write_sto_file(app.grf_forces_adjusted, grfFileName);
                    end

                    grforces = xml_read('GRF_file_all.xml');
                    grforces_generated = xml_read('GRF_file_empty.xml');

                    % remove duplicate assignments to avoid applying GRF
                    % multiple times
                    uniqueRightStepForcePlateAssignment = [];
                    uniqueLeftStepForcePlateAssignment = [];
                    if ~isempty(outputSets.(outputFolders{o}).rightStepForcePlateAssignment)
                        [~, uniqueIdx] = unique(outputSets.(outputFolders{o}).rightStepForcePlateAssignment(:,2), 'stable');
                        uniqueRightStepForcePlateAssignment = outputSets.(outputFolders{o}).rightStepForcePlateAssignment(uniqueIdx, :);
                    end
                    if ~isempty(outputSets.(outputFolders{o}).leftStepForcePlateAssignment)
                        [~, uniqueIdx] = unique(outputSets.(outputFolders{o}).leftStepForcePlateAssignment(:,2), 'stable');
                        uniqueLeftStepForcePlateAssignment = outputSets.(outputFolders{o}).leftStepForcePlateAssignment(uniqueIdx, :);
                    end

                    counter = 1;
                    if app.Ignorec3deventsCheckBox.Value == 0
                        if app.AMTITandemTreadmillCheckBox.Value || ~app.mirrorCheckBox.Value
                            for i = 1 : size(uniqueRightStepForcePlateAssignment, 1)
                                grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(uniqueRightStepForcePlateAssignment(i, 2));
                                counter = counter + 1;
                            end
                            for i = 1 : size(uniqueLeftStepForcePlateAssignment, 1)
                                grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(9 + uniqueLeftStepForcePlateAssignment(i, 2));
                                counter = counter + 1;
                            end
                        else % this is when data should be mirrored
                            for i = 1 : size(uniqueRightStepForcePlateAssignment, 1)
                                grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(9 + uniqueRightStepForcePlateAssignment(i, 2));
                                counter = counter + 1;
                            end
                            for i = 1 : size(uniqueLeftStepForcePlateAssignment, 1)
                                grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(uniqueLeftStepForcePlateAssignment(i, 2));
                                counter = counter + 1;
                            end
                        end
                        Pref = struct;
                        Pref.StructItem = false;

                        d.Value = d.Value + 0.1 / numel(outputFolders);
                        d.Message = 'writing GRF.XML ...';
                        xml_write(fullfile(output_folder, 'GRF.xml'), grforces_generated, 'OpenSimDocument', Pref);
                    else
                        for i = 1 : size(app.rightListBox.Value, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(str2double(app.rightListBox.Value{i}(end)));
                            counter = counter + 1;
                        end

                        for i = 1 : size(app.leftListBox.Value, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(str2double(app.leftListBox.Value{i}(end)) + 9);
                            counter = counter + 1;
                        end

                        Pref = struct;
                        Pref.StructItem = false;

                        d.Value = d.Value + 0.1 / numel(outputFolders);
                        d.Message = 'writing GRF.XML ...';
                        xml_write(fullfile(output_folder, 'GRF.xml'), grforces_generated, 'OpenSimDocument', Pref);
                    end

                    d.Value = d.Value + 0.1 / numel(outputFolders);
                    cycle = struct;
                    if app.Ignorec3deventsCheckBox.Value == 0
                        d.Message = 'preparing information about cycles ...';
                        cycle = struct;
                        if ~isempty(outputSets.(outputFolderName).leftICs)
                            cycle.left = struct;
                            cycle.left.start = outputSets.(outputFolderName).leftICs - app.firstFrame_;
                            cycle.left.startFrame = outputSets.(outputFolderName).leftICs - (app.firstFrame_-1);
                            cycle.left.startTime = cycle.left.startFrame / app.frequency_;

                            footStrikes = str2double(app.leftFootStrikeList) - (app.firstFrame_-1);
                            contraLateralFootStrikes = str2double(app.rightFootStrikeList) - (app.firstFrame_ - 1);
                            for i = 1 : size(cycle.left.startFrame, 2)
                                cycle.left.end(i) = footStrikes( find( footStrikes > cycle.left.start(i)+2, 1 ) ) - 1;
                                cycle.left.endFrame(i) = footStrikes( find( footStrikes > cycle.left.startFrame(i), 1 ) );
                                cycle.left.endTime(i) = cycle.left.endFrame(i) / app.frequency_;
                                if ~isempty(app.leftFootOffList)
                                    try
                                        cycle.left.footOff(i) = app.leftFootOffList( find( app.leftFootOffList > cycle.left.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                        cycle.left.footOffFrame(i) = app.leftFootOffList( find( app.leftFootOffList  > cycle.left.startFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                        cycle.left.footOffTime(i) = cycle.left.footOffFrame(i) / app.frequency_;
                                    catch e
                                        if contains(e.message, 'Unable to perform assignment because the left and right sides have a different number of elements')
                                            writeLog(app, ['No Foot Off Event for left step starting at frame ' num2str(cycle.left.startFrame(i) + app.firstFrame_)]);
                                        end
                                        cycle.left.footOff(i) = cycle.left.end(i);
                                        cycle.left.footOffFrame(i) = cycle.left.endFrame(i);
                                        cycle.left.footOffTime(i) = cycle.left.footOffFrame(i) / app.frequency_;
                                    end
                                    if strcmp(app.MovementTypeDropDown.Value, 'SideStep')
                                        cycle.left.touchDown(i) = cycle.left.footOff(i);
                                        cycle.left.touchDownFrame(i) = cycle.left.footOffFrame(i);
                                        cycle.left.touchDownTime(i) = cycle.left.touchDownFrame(i) / app.frequency_;
                                        cycle.left.leaveGround(i) = app.leftFootOffList( find( app.leftFootOffList  > cycle.left.footOff(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                        cycle.left.leaveGroundFrame(i) = app.leftFootOffList( find( app.leftFootOffList  > cycle.left.footOffFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                        cycle.left.leaveGroundTime(i) = cycle.left.leaveGroundFrame(i) / app.frequency_;
                                    end
                                end
                                if ~isempty(app.rightFootOffList)
                                    cycle.left.contraLateralFootOff(i) = app.rightFootOffList( find( app.rightFootOffList > cycle.left.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                    cycle.left.contraLateralFootOffFrame(i) = app.rightFootOffList( find( app.rightFootOffList  > cycle.left.startFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                    cycle.left.contraLateralFootOffTime(i) = cycle.left.contraLateralFootOffFrame(i) / app.frequency_;
                                end
                                if ~isempty(contraLateralFootStrikes)
                                    cycle.left.contraLateralFootStrike(i) = contraLateralFootStrikes( find( contraLateralFootStrikes > cycle.left.start(i) + 2, 1 ) ) - 1;
                                    cycle.left.contraLateralFootStrikeFrame(i) = contraLateralFootStrikes( find( contraLateralFootStrikes > cycle.left.startFrame(i), 1 ) );
                                    cycle.left.contraLateralFootStrikeTime(i) = cycle.left.contraLateralFootStrikeFrame(i) / app.frequency_;
                                end
                                if isfield(cycle.left, 'contraLateralFootStrike') && length(cycle.left.contraLateralFootStrike) >= i && ...
                                        isfield(cycle.left, 'contraLateralFootOff') && length(cycle.left.contraLateralFootOff) >= i
                                    % set a flag whether contralateral events are within start and end of the step, otherwise they are probably not correct
                                    if cycle.left.contraLateralFootOff(i) > cycle.left.start(i) && cycle.left.contraLateralFootOff(i) < cycle.left.end(i) && ...
                                            cycle.left.contraLateralFootStrike(i) > cycle.left.start(i) && cycle.left.contraLateralFootStrike(i) < cycle.left.end(i) && ...
                                            cycle.left.contraLateralFootOff(i) < cycle.left.contraLateralFootStrike(i)
                                        cycle.left.contraLateralEventsValid(i) = 1;
                                    else
                                        cycle.left.contraLateralEventsValid(i) = 0;
                                    end
                                else
                                    cycle.left.contraLateralEventsValid(i) = 0;
                                end
                            end
                        end
                        if ~isempty(outputSets.(outputFolderName).rightICs)
                            cycle.right = struct;
                            cycle.right.start = outputSets.(outputFolderName).rightICs - app.firstFrame_;
                            cycle.right.startFrame = outputSets.(outputFolderName).rightICs - (app.firstFrame_-1);
                            cycle.right.startTime = cycle.right.startFrame / app.frequency_;

                            footStrikes = str2double(app.rightFootStrikeList) - (app.firstFrame_ - 1);
                            contraLateralFootStrikes = str2double(app.leftFootStrikeList) - (app.firstFrame_ - 1);
                            for i = 1 : size(cycle.right.startFrame, 2)
                                cycle.right.end(i) = footStrikes( find( footStrikes > cycle.right.start(i)+2, 1 ) ) - 1;
                                cycle.right.endFrame(i) = footStrikes( find( footStrikes > cycle.right.startFrame(i), 1 ) );
                                cycle.right.endTime(i) = cycle.right.endFrame(i) / app.frequency_;
                                if ~isempty(app.rightFootOffList)
                                    try
                                        cycle.right.footOff(i) = app.rightFootOffList( find( app.rightFootOffList > cycle.right.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                        cycle.right.footOffFrame(i) = app.rightFootOffList( find( app.rightFootOffList  > cycle.right.startFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                        cycle.right.footOffTime(i) = cycle.right.footOffFrame(i) / app.frequency_;
                                    catch e
                                        if contains(e.message, 'Unable to perform assignment because the left and right sides have a different number of elements')
                                            writeLog(app, ['No Foot Off Event for right step starting at frame ' num2str(cycle.right.startFrame(i) + app.firstFrame_)]);
                                        end
                                        cycle.right.footOff(i) = cycle.right.end(i);
                                        cycle.right.footOffFrame(i) = cycle.right.endFrame(i);
                                        cycle.right.footOffTime(i) = cycle.right.footOffFrame(i) / app.frequency_;
                                    end
                                    if strcmp(app.MovementTypeDropDown.Value, 'SideStep')
                                        cycle.right.touchDown(i) = cycle.right.footOff(i);
                                        cycle.right.touchDownFrame(i) = cycle.right.footOffFrame(i);
                                        cycle.right.touchDownTime(i) = cycle.right.touchDownFrame(i) / app.frequency_;
                                        cycle.right.leaveGround(i) = app.rightFootOffList( find( app.rightFootOffList  > cycle.right.footOff(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                        cycle.right.leaveGroundFrame(i) = app.rightFootOffList( find( app.rightFootOffList  > cycle.right.footOffFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                        cycle.right.leaveGroundTime(i) = cycle.right.leaveGroundFrame(i) / app.frequency_;
                                    end
                                end
                                if ~isempty(app.leftFootOffList)
                                    cycle.right.contraLateralFootOff(i) = app.leftFootOffList( find( app.leftFootOffList > cycle.right.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                                    cycle.right.contraLateralFootOffFrame(i) = app.leftFootOffList( find( app.leftFootOffList  > cycle.right.startFrame(i) + app.firstFrame_, 1 ) ) - (app.firstFrame_-2);
                                    cycle.right.contraLateralFootOffTime(i) = cycle.right.contraLateralFootOffFrame(i) / app.frequency_;
                                end
                                if ~isempty(contraLateralFootStrikes)
                                    cycle.right.contraLateralFootStrike(i) = contraLateralFootStrikes( find( contraLateralFootStrikes > cycle.right.start(i) + 2, 1 ) ) - 1;
                                    cycle.right.contraLateralFootStrikeFrame(i) = contraLateralFootStrikes( find( contraLateralFootStrikes > cycle.right.startFrame(i), 1 ) );
                                    cycle.right.contraLateralFootStrikeTime(i) = cycle.right.contraLateralFootStrikeFrame(i) / app.frequency_;
                                end
                                if isfield(cycle.right, 'contraLateralFootStrike') && length(cycle.right.contraLateralFootStrike) >= i && ...
                                        isfield(cycle.right, 'contraLateralFootOff') && length(cycle.right.contraLateralFootOff) >= i
                                    % set a flag whether contralateral events are within start and end of the step, otherwise they are probably not correct
                                    if cycle.right.contraLateralFootOff(i) > cycle.right.start(i) && cycle.right.contraLateralFootOff(i) < cycle.right.end(i) && ...
                                            cycle.right.contraLateralFootStrike(i) > cycle.right.start(i) && cycle.right.contraLateralFootStrike(i) < cycle.right.end(i) && ...
                                            cycle.right.contraLateralFootOff(i) < cycle.right.contraLateralFootStrike(i)
                                        cycle.right.contraLateralEventsValid(i) = 1;
                                    else
                                        cycle.right.contraLateralEventsValid(i) = 0;
                                    end
                                else
                                    cycle.right.contraLateralEventsValid(i) = 0;
                                end
                            end
                        end

                        if app.mirrorCheckBox.Value
                            cycle_tmp = cycle;
                            cycle = struct;
                            if isfield(cycle_tmp, 'right')
                                cycle.left = cycle_tmp.right;
                            end
                            if isfield(cycle_tmp, 'left')
                                cycle.right = cycle_tmp.left;
                            end
                        end
                    end

                    firstFrame = app.firstFrame_;
                    frequency = app.frequency_;
                    forceFrequency = app.forceFrequency;

                    minStartTime = 10000;
                    maxEndTime = 0;

                    if isfield(cycle, 'left')
                        minStartTime = min([minStartTime, cycle.left.startTime]);
                        maxEndTime = max([maxEndTime cycle.left.endTime]);
                    end
                    if isfield(cycle, 'right')
                        minStartTime = min([minStartTime, cycle.right.startTime]);
                        maxEndTime = max([maxEndTime cycle.right.endTime]);
                    end

                    % duration = (btkGetPointFrameNumber(app.acq) - 1) / frequency;
                    duration = maxEndTime - minStartTime;

                    c3dExportSettings = struct;
                    c3dExportSettings.markersFiltered = app.filtermarkersCheckBox.Value;
                    if c3dExportSettings.markersFiltered
                        c3dExportSettings.markersFilterOrder = app.filterorderSpinner_Marker.Value;
                        c3dExportSettings.markersFilterCutOff = app.cutofffrequencySpinner_Marker.Value;
                    end
                    c3dExportSettings.GRFsFiltered = app.filterGRFsCheckBox.Value;
                    if c3dExportSettings.GRFsFiltered
                        c3dExportSettings.GRFFilterOrder = app.filterorderSpinner.Value;
                        c3dExportSettings.GRFFilterCutOff = app.cutofffrequencySpinner.Value;
                    end
                    c3dExportSettings.EMGFiltered = app.filterEMGCheckBox.Value;
                    if c3dExportSettings.EMGFiltered
                        % one could include emg filtering in the gui and add the settings here
                        c3dExportSettings.EMGFilter = 'Bandpass [20-400]; De-Meaned; Full-wave rectified; 4th order butterworth low-pass with 10Hz cutoff; Moving Average with Windows of 50ms; Negative Values set to zero';
                        c3dExportSettings.EMG_normalizedTo1 = app.normalizeallmagnitudesto1CheckBox.Value;
                        c3dExportSettings.EMG_multipliedBy = app.EMGfactorSpinner.Value;
                        c3dExportSettings.EMG_renamed = app.renameEMGlabelsCheckBox.Value;
                        if c3dExportSettings.EMG_renamed
                            c3dExportSettings.EMG_renamingCSV = app.emgLabelCSVPath;
                        end
                    end

                    c3dExportSettings.dataRotatedAroundY = app.rotatearoundyaxisDropDown.Value;
                    c3dExportSettings.movementType = app.MovementTypeDropDown.Value;

                    if app.AMTITandemTreadmillCheckBox.Value
                        c3dExportSettings.isAMTITreadmill = 1;
                    end

                    c3dExportSettings.exportEachStepSeperately = app.exporteachstepseparatelyCheckBox.Value;
                    c3dExportSettings.croppedToFramesOfInterest = app.croptovalidfootstrikesbufferCheckBox.Value;
                    if c3dExportSettings.croppedToFramesOfInterest
                        c3dExportSettings.timeBufferStart = app.BufferStartsSpinner.Value;
                        c3dExportSettings.timeBufferEnd = app.BufferEndsSpinner.Value;
                    end

                    d.Value = d.Value + 0.1 / numel(outputFolders);
                    d.Message = 'write settings.mat and copy c3d file ...';
                    save(fullfile(output_folder, 'settings.mat'), 'c3dExportSettings', 'cycle', 'firstFrame', 'forceFrequency', 'frequency', 'duration', '-mat');
                    % copyfile(fullfile(app.path, app.c3dFileName), fullfile(output_folder, 'c3dfile.c3d'));

                    d.Value = d.Value + 0.1 / numel(outputFolders);

                    % export filtered EMG
                    invalidFieldNames = 0;
                    if app.filterEMGCheckBox.Value
                        if app.renameEMGlabelsCheckBox.Value
                            if ~isempty(app.emgLabelCSVPath)
                                exportEMG = struct;
                                channels = fieldnames(app.EMG_lowPass);
                                channels = sort(channels);
                                channels = channels(~strcmp(channels, 'time'));
                                emgLabels = readtable(app.emgLabelCSVPath,'HeaderLines',0);
                                for i = 1 : size(channels, 1)
                                    originalLabel = table2cell(emgLabels(i, 1));
                                    newLabel = table2cell(emgLabels(i, 2));
                                    if ~strcmp(newLabel{1}, '')
                                        if isvarname(newLabel{1})
                                            exportEMG.(newLabel{1}) = app.EMG_lowPass.(originalLabel{1});
                                        else
                                            errordlg([newLabel{1} ' is an invalid name in MATLAB! Change this and create files again! Name must not include special characters and must not start with a number!']);
                                            invalidFieldNames = 1;
                                            break;
                                        end
                                    end
                                end
                                exportEMG.time = app.EMG.time;
                                if invalidFieldNames == 1
                                    exportEMG = app.EMG_lowPass;
                                end
                            else
                                errordlg('EMG labels were not renamed! No CSV was selected! Exporting was continued with default names');
                                exportEMG = app.EMG_lowPass;
                                exportEMG.time = app.EMG.time;
                            end
                        else
                            exportEMG = app.EMG_lowPass;
                            exportEMG.time = app.EMG.time;
                        end

                        fields = fieldnames(exportEMG);
                        fields = fields(~contains(fields, 'time'));
                        for f = 1 : numel(fields)
                            exportEMG.(fields{f}) = exportEMG.(fields{f}) * app.EMGfactorSpinner.Value;
                        end

                        if app.croptovalidfootstrikesbufferCheckBox.Value == 1
                            fields = fieldnames(exportEMG);
                            keepIdx = and(exportEMG.time >= exportStartTime, exportEMG.time <= exportEndTime);
                            for i = 1 : numel(fields)
                                exportEMG.(fields{i}) = exportEMG.(fields{i})(keepIdx);
                            end
                        end

                        write_sto_file(exportEMG, fullfile(output_folder, 'EMG_filtered.sto'))
                    end

                    % export raw emg
                    invalidFieldNames = 0;
                    if app.exportrawEMGCheckBox.Value
                        if app.renameEMGlabelsCheckBox.Value
                            if ~isempty(app.emgLabelCSVPath)
                                exportEMG = struct;
                                channels = fieldnames(app.EMG_lowPass); % export the same fields as for the filtered file
                                channels = sort(channels);
                                channels = channels(~strcmp(channels, 'time'));
                                emgLabels = readtable(app.emgLabelCSVPath,'HeaderLines',0);
                                for i = 1 : size(channels, 1)
                                    originalLabel = table2cell(emgLabels(i, 1));
                                    newLabel = table2cell(emgLabels(i, 2));
                                    if ~strcmp(newLabel{1}, '')
                                        if isvarname(newLabel{1})
                                            exportEMG.(newLabel{1}) = app.EMG.(originalLabel{1});
                                        else
                                            errordlg([newLabel{1} ' is an invalid name in MATLAB! Change this and create files again! Name must not include special characters and must not start with a number!']);
                                            invalidFieldNames = 1;
                                            break;
                                        end
                                    end
                                end
                                exportEMG.time = app.EMG.time;
                                if invalidFieldNames == 1
                                    exportEMG = app.EMG;
                                end
                            else
                                errordlg('EMG labels were not renamed! No CSV was selected! Exporting was continued with default names');
                                exportEMG = app.EMG;
                                exportEMG.time = app.EMG.time;
                            end
                        else
                            exportEMG = app.EMG;
                            exportEMG.time = app.EMG.time;
                        end

                        if isfield(exportEMG, 'frequency')
                            exportEMG = rmfield(exportEMG, 'frequency');
                        end

                        fields = fieldnames(exportEMG);
                        fields = fields(~contains(fields, 'time'));
                        for f = 1 : numel(fields)
                            exportEMG.(fields{f}) = exportEMG.(fields{f}) * app.EMGfactorSpinner.Value;
                        end

                        if app.croptovalidfootstrikesbufferCheckBox.Value == 1
                            fields = fieldnames(exportEMG);
                            keepIdx = and(exportEMG.time >= exportStartTime, exportEMG.time <= exportEndTime);
                            for i = 1 : numel(fields)
                                exportEMG.(fields{i}) = exportEMG.(fields{i})(keepIdx);
                            end
                        end

                        write_sto_file(exportEMG, fullfile(output_folder, 'EMG_raw.sto'))
                    end

                    if app.ClimbingCheckBox.Value % EMG-channel Fz_left Fy_left Fz_right Fy_right are stored in seperate Force file and special XML will be written to apply forces to body
                        forces = struct;

                        marker = load_marker_trc(fullfile(output_folder, 'marker_experimental.trc'));
                        markerLength = length(marker.RFIN_X);
                        forceLength = length(app.EMG.time);

                        tmp_fields = fieldnames(app.EMG);
                        forces.time = app.EMG.time;

                        forces.ground_force_1_vx = zeros(length(forces.time), 1);
                        if xor(app.invertforcesensorsCheckBox.Value, app.mirrorCheckBox.Value)
                            forces.ground_force_1_vy = app.EMG.(tmp_fields{contains(tmp_fields, 'Fz_right')});
                            if app.mirrorCheckBox.Value
                                forces.ground_force_1_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_right')});
                            else
                                forces.ground_force_1_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_right')}) * -1;
                            end
                        else
                            forces.ground_force_1_vy = app.EMG.(tmp_fields{contains(tmp_fields, 'Fz_left')});
                            if app.mirrorCheckBox.Value
                                forces.ground_force_1_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_left')});
                            else
                                forces.ground_force_1_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_left')}) * -1;
                            end
                        end
                        forces.ground_force_1_px = interp1(1:markerLength, cell2mat(marker.LFIN_X), linspace(1,markerLength,forceLength))';
                        forces.ground_force_1_py = interp1(1:markerLength, cell2mat(marker.LFIN_Y), linspace(1,markerLength,forceLength))';
                        forces.ground_force_1_pz = interp1(1:markerLength, cell2mat(marker.LFIN_Z), linspace(1,markerLength,forceLength))';
                        forces.ground_moment_1_mx = zeros(length(forces.time), 1);
                        forces.ground_moment_1_my = zeros(length(forces.time), 1);
                        forces.ground_moment_1_mz = zeros(length(forces.time), 1);
                        forces.ground_force_2_vx = zeros(length(forces.time), 1);
                        if xor(app.invertforcesensorsCheckBox.Value, app.mirrorCheckBox.Value)
                            forces.ground_force_2_vy = app.EMG.(tmp_fields{contains(tmp_fields, 'Fz_left')});
                            if app.mirrorCheckBox.Value
                                forces.ground_force_2_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_left')});
                            else
                                forces.ground_force_2_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_left')}) * -1;
                            end
                        else
                            forces.ground_force_2_vy = app.EMG.(tmp_fields{contains(tmp_fields, 'Fz_right')});
                            if app.mirrorCheckBox.Value
                                forces.ground_force_2_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_right')});
                            else
                                forces.ground_force_2_vz = app.EMG.(tmp_fields{contains(tmp_fields, 'Fy_right')}) * -1;
                            end
                        end
                        forces.ground_force_2_px = interp1(1:markerLength, cell2mat(marker.RFIN_X), linspace(1,markerLength,forceLength))';
                        forces.ground_force_2_py = interp1(1:markerLength, cell2mat(marker.RFIN_Y), linspace(1,markerLength,forceLength))';
                        forces.ground_force_2_pz = interp1(1:markerLength, cell2mat(marker.RFIN_Z), linspace(1,markerLength,forceLength))';
                        forces.ground_moment_2_mx = zeros(length(forces.time), 1);
                        forces.ground_moment_2_my = zeros(length(forces.time), 1);
                        forces.ground_moment_2_mz = zeros(length(forces.time), 1);


                        if app.croptovalidfootstrikesbufferCheckBox.Value == 1
                            fields = fieldnames(forces);
                            keepIdx = and(forces.time >= exportStartTime, forces.time <= exportEndTime);
                            for i = 1 : numel(fields)
                                forces.(fields{i}) = forces.(fields{i})(keepIdx);
                            end
                        end

                        write_sto_file(forces, fullfile(output_folder, 'forces.mot'));

                        copyfile('climbing_forces.xml', fullfile(output_folder, 'GRF.xml'));
                    end
                end

                app.croptovalidfootstrikesbufferCheckBox.Value = cropCheckBoxOriginal;

                d.Value = 1;
                d.Message = 'Finished';
                writeLog(app, ['Files for ' app.c3dFileName ' sucessfully written']);
                close(d);
            end
        end

        % Button pushed function: openselectedfilewithdefaultprogramButton
        function openselectedfilewithdefaultprogramButtonPushed(app, event)
            if app.fileLoaded
                winopen(fullfile(app.path, app.c3dFileName));
            end
        end

        % Value changed function: leftListBox, rightListBox
        function updateCycles(app, event)
            if app.Ignorec3deventsCheckBox.Value == 0
                [~,I] = sort(str2double(app.leftListBox.Value));
                app.leftListBox.Value = app.leftListBox.Value(I);
                [~,I] = sort(str2double(app.rightListBox.Value));
                app.rightListBox.Value = app.rightListBox.Value(I);
                if app.fileLoaded && app.c3d.getNumForces() > 0
                    selectedFrames = cellfun(@str2double, app.rightListBox.Value);
                    app.rightStepForcePlateAssignment = app.rightStepForcePlateAssignmentAll;
                    if ~isempty(app.rightStepForcePlateAssignment)
                        isStillSelected = ismember(app.rightStepForcePlateAssignment(:, 1), selectedFrames);
                    else
                        isStillSelected = [];
                    end
                    app.rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(isStillSelected, :);

                    selectedFrames = cellfun(@str2double, app.leftListBox.Value);
                    app.leftStepForcePlateAssignment = app.leftStepForcePlateAssignmentAll;

                    if ~isempty(app.leftStepForcePlateAssignment)
                        isStillSelected = ismember(app.leftStepForcePlateAssignment(:, 1), selectedFrames);
                    else
                        isStillSelected = [];
                    end
                    app.leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(isStillSelected, :);

                    updateCyclePlot(app);
                    updateGRF_Plot(app, 0);
                end
            end
        end

        % Value changed function: cutofffrequencySpinner, 
        % ...and 2 other components
        function updateGRF_Plot(app, event)
            cla(app.UIAxesGRFs);

            try
                if app.fileLoaded && app.c3d.getNumForces() > 0
                    for i = 1 : (numel(fieldnames(app.grf_forces)) - 1) / 9
                        if app.AMTITandemTreadmillCheckBox.Value == 0 || (app.AMTITandemTreadmillCheckBox.Value == 1 && (i == 7 || i ==8))
                            plot(app.UIAxesGRFs, app.grf_forces.(['ground_force_' num2str(i) '_vy']));
                        end
                    end

                    filterGRFs(app)

                    verticalForceFields = fieldnames(app.grf_forces_adjusted);
                    verticalForceFields = verticalForceFields(contains(verticalForceFields, '_vy'));
                    if app.AMTITandemTreadmillCheckBox.Value == 1
                        verticalForceFields = verticalForceFields(or(contains(verticalForceFields, '_7'), contains(verticalForceFields, '_8')));
                    end

                    for i = 1 : numel(verticalForceFields)
                        plot(app.UIAxesGRFs, app.grf_forces_adjusted.(verticalForceFields{i}));
                    end
                    ylim(app.UIAxesGRFs, 'auto');
                    xlimitSteps = xlim(app.UIAxesSteps);
                    %                 maxLim = (xlimitSteps(2) - xlimitSteps(1)) / app.frequency_ * app.forceFrequency;
                    maxLim = (xlimitSteps(2) - app.firstFrame_) / app.frequency_ * app.forceFrequency;
                    minLim = (xlimitSteps(1) - app.firstFrame_) / app.frequency_ * app.forceFrequency;

                    xlim(app.UIAxesGRFs, [minLim maxLim])
                end
            catch e
                disp(e.message);
            end
        end

        % Value changed function: detectFPbyMarker, 
        % ...and 1 other component
        function detectforceplatesautomaticallyCheckBoxValueChanged(app, event)
            app.leftStepForcePlateAssignment = [];
            app.rightStepForcePlateAssignment = [];
            if app.AMTITandemTreadmillCheckBox.Value
                % if size(app.leftListBox.Value, 1) >= 1
                %     app.leftStepForcePlateAssignment(end+1, 1) = str2double(app.leftListBox.Value{1});
                %     app.leftStepForcePlateAssignment(end, 2) = 7; % FP7 is always left leg
                % end
                % if size(app.rightListBox.Value, 1) >= 1
                %     app.rightStepForcePlateAssignment(end+1, 1) = str2double(app.rightListBox.Value{1});
                %     app.rightStepForcePlateAssignment(end, 2) = 8; % FP8 is always right leg
                % end

                for i = 1 : size(app.leftListBox.Value, 2)
                    app.leftStepForcePlateAssignment(end+1, 1) = str2double(app.leftListBox.Value{i});
                    app.leftStepForcePlateAssignment(end, 2) = 7; % FP7 is always left leg
                end
                for i = 1 : size(app.rightListBox.Value, 2)
                    app.rightStepForcePlateAssignment(end+1, 1) = str2double(app.rightListBox.Value{i});
                    app.rightStepForcePlateAssignment(end, 2) = 8; % FP8 is always right leg
                end

                app.rightStepForcePlateAssignmentAll = app.rightStepForcePlateAssignment;
                app.leftStepForcePlateAssignmentAll = app.leftStepForcePlateAssignment;
                updateCyclePlot(app);
            else
                if app.detectforceplatesautomaticallyCheckBox.Value
                    if app.fileLoaded
                        app.c3d.rotateData('x', 90);
                        temp_markers = osimTableToStruct(app.c3d.getTable_markers);
                        if isfield(temp_markers, app.detectFPbyMarker.Value)
                            app.detectFPbyMarker.BackgroundColor = 'white';
                            app.detectFPbyMarker.FontColor = 'black';
                            app.leftStepForcePlateAssignment = [];
                            app.rightStepForcePlateAssignment = [];
                            [forceplates, ~] = btkGetForcePlatforms(app.acq);

                            for i = 1 : size(app.leftListBox.Value, 2)
                                if strcmp(app.MovementTypeDropDown.Value, 'SideStep')
                                    app.leftStepForcePlateAssignment(end+1, 1) = str2double(app.leftListBox.Value{i});
                                    app.leftStepForcePlateAssignment(end, 2) = 3;
                                else
                                    % find closest forceplate to heel marker
                                    heelLocation = temp_markers.(app.detectFPbyMarker.Value)(str2double(app.leftListBox.Value{i}) - app.firstFrame_ + 20, :);
                                    if max(heelLocation) < 10
                                        heelLocation = heelLocation * 1000;
                                    end
                                    isset = 0;
                                    for j = 1 : size(forceplates, 1)
                                        fpXlimits = sort(forceplates(j).corners(1, :));
                                        fpYlimits = sort(forceplates(j).corners(2, :));
                                        if heelLocation(1) > fpXlimits(1) && heelLocation(1) < fpXlimits(3) ...
                                                && heelLocation(2) > fpYlimits(1) && heelLocation(2) < fpYlimits(3)
                                            isset = 1;
                                            break;
                                        end
                                    end
                                    if isset
                                        app.leftStepForcePlateAssignment(end+1, 1) = str2double(app.leftListBox.Value{i});
                                        app.leftStepForcePlateAssignment(end, 2) = j;
                                    end
                                end
                            end
                            for i = 1 : size(app.rightListBox.Value, 2)
                                if strcmp(app.MovementTypeDropDown.Value, 'SideStep')
                                    app.rightStepForcePlateAssignment(end+1, 1) = str2double(app.rightListBox.Value{i});
                                    app.rightStepForcePlateAssignment(end, 2) = 4;
                                else

                                    if strcmp(app.detectFPbyMarker.Value(1), 'L')
                                        rightMarkerName = ['R' app.detectFPbyMarker.Value(2:end)];
                                    else
                                        rightMarkerName = app.detectFPbyMarker.Value;
                                    end

                                    % find closest forceplate to heel marker
                                    heelLocation = temp_markers.(rightMarkerName)(str2double(app.rightListBox.Value{i}) - app.firstFrame_ + 20, :);
                                    if max(heelLocation) < 10
                                        heelLocation = heelLocation * 1000;
                                    end
                                    isset = 0;
                                    for j = 1 : size(forceplates, 1)
                                        fpXlimits = sort(forceplates(j).corners(1, :));
                                        fpYlimits = sort(forceplates(j).corners(2, :));
                                        if heelLocation(1) > fpXlimits(1) && heelLocation(1) < fpXlimits(3) ...
                                                && heelLocation(2) > fpYlimits(1) && heelLocation(2) < fpYlimits(3)
                                            isset = 1;
                                            break;
                                        end
                                    end
                                    if isset
                                        app.rightStepForcePlateAssignment(end+1, 1) = str2double(app.rightListBox.Value{i});
                                        app.rightStepForcePlateAssignment(end, 2) = j;
                                    end
                                end
                            end
                        else
                            app.detectFPbyMarker.Value = 'SET THIS!';
                            app.detectFPbyMarker.BackgroundColor = 'red';
                            app.detectFPbyMarker.FontColor = 'white';
                            app.detectforceplatesautomaticallyCheckBox.Value = 0;
                        end

                        app.c3d.rotateData('x', -90);

                        enfFiles = dir(fullfile(app.path, strrep(app.c3dFileName, '.c3d', '.*enf')));
                        if numel(enfFiles) == 1
                            filename = fullfile(enfFiles(1).folder, enfFiles(1).name);
                            enfText = fileread(filename);
                            keepIdx = [];
                            for i = 1 : size(app.rightStepForcePlateAssignment, 1)
                                if contains(lower(enfText), ['fp' num2str(app.rightStepForcePlateAssignment(i, 2)) '=right'])
                                    keepIdx(end+1) = i;
                                end
                            end
                            if ~isempty(keepIdx)
                                keepValues = app.rightStepForcePlateAssignment(keepIdx, 1);
                                keepIdxListBox = ismember(str2double(app.rightListBox.Value), keepValues);
                                app.rightListBox.Value = app.rightListBox.Value(keepIdxListBox);
                                app.rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(keepIdx, :);
                            else
                                app.rightListBox.Value = {};
                                app.rightStepForcePlateAssignment = [];
                            end
                            keepIdx = [];
                            for i = 1 : size(app.leftStepForcePlateAssignment, 1)
                                if contains(lower(enfText), ['fp' num2str(app.leftStepForcePlateAssignment(i, 2)) '=left'])
                                    keepIdx(end+1) = i;
                                end
                            end
                            if ~isempty(keepIdx)
                                keepValues = app.leftStepForcePlateAssignment(keepIdx, 1);
                                keepIdxListBox = ismember(str2double(app.leftListBox.Value), keepValues);
                                app.leftListBox.Value = app.leftListBox.Value(keepIdxListBox);
                                app.leftStepForcePlateAssignment = app.leftStepForcePlateAssignment(keepIdx, :);
                            else
                                app.leftListBox.Value = {};
                                app.leftStepForcePlateAssignment = [];
                            end
                        end

                        app.rightStepForcePlateAssignmentAll = app.rightStepForcePlateAssignment;
                        app.leftStepForcePlateAssignmentAll = app.leftStepForcePlateAssignment;
                        updateCyclePlot(app);
                    end
                else
                    app.leftStepForcePlateAssignment = [];
                    app.rightStepForcePlateAssignment = [];
                    app.leftStepForcePlateAssignmentAll = [];
                    app.rightStepForcePlateAssignmentAll = [];
                    updateCyclePlot(app);
                end
            end
        end

        % Value changed function: filterEMGCheckBox
        function filterEMGCheckBoxValueChanged(app, event)
            value = app.filterEMGCheckBox.Value;
            if value
                readEMG(app);
                filterEMG(app);
            else
                app.EMG = [];
                app.EMG_lowPass = [];
            end
        end

        % Button pushed function: SelectexistingCSVButton
        function SelectexistingCSVButtonPushed(app, event)
            if ~isempty(app.path)
                [filename, tmp_path] = uigetfile(fullfile(app.path, '*.csv'));
            else
                [filename, tmp_path] = uigetfile('*.csv');
            end
            app.emgLabelCSVPath = fullfile(tmp_path, filename);
        end

        % Button pushed function: createCSVandopenButton
        function createCSVandopenButtonPushed(app, event)
            if ~isempty(app.EMG_lowPass)
                [filename, tmp_path] = uiputfile(fullfile(app.path, 'emg_labels.csv'));
                channels = fieldnames(app.EMG_lowPass);
                channels = sort(channels);
                for i = 1 : size(channels, 1)
                    channels{i, 2} = channels{i, 1};
                end
                writetable(cell2table(channels, 'VariableNames', {'original label', 'new label'}), fullfile(tmp_path, filename), 'Delimiter', ';');
                app.emgLabelCSVPath = fullfile(tmp_path, filename);
                winopen(app.emgLabelCSVPath);
            end
        end

        % Value changed function: Ignorec3deventsCheckBox
        function Ignorec3deventsCheckBoxValueChanged(app, event)
            value = app.Ignorec3deventsCheckBox.Value;
            if value
                app.forceplatecontactsLabel.Visible = "on";
                app.validfootstrikeeventsLabel.Visible = "off";
                app.detectforceplatesautomaticallyCheckBox.Value = 0;
            else
                app.forceplatecontactsLabel.Visible = "off";
                app.validfootstrikeeventsLabel.Visible = "on";
                app.detectforceplatesautomaticallyCheckBox.Value = 1;
                answer = MFquestdlg([0.4, 0.4], 'Reprocess the c3d file to identify events and foot contacts?', 'Reprocess?', 'Yes', 'No', 'Yes');
                if strcmp(answer, 'Yes')
                    app.processC3Dfile();
                end
            end
        end

        % Value changed function: AMTITandemTreadmillCheckBox
        function AMTITandemTreadmillCheckBoxValueChanged(app, event)
            % value = app.AMTITandemTreadmillCheckBox.Value;
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
            app.rotatearoundyaxisDropDown.Value = '270°';
            app.rotatearoundyaxisDropDown.Enable = "on";
            app.detectforceplatesautomaticallyCheckBox.Value = 0;
        end

        % Value changed function: removeblanksetcfromfilenameCheckBox
        function removeblanksetcfromfilenameCheckBoxValueChanged(app, event)
            value = app.removeblanksetcfromfilenameCheckBox.Value;
            app.processC3Dfile()
        end
    end

    % Component initialization
    methods (Access = public)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create PrepareC3DfileUIFigure and hide until all components are created
            app.PrepareC3DfileUIFigure = uifigure('Visible', 'off');
            app.PrepareC3DfileUIFigure.Position = [10 10 1517 917];
            app.PrepareC3DfileUIFigure.Name = 'PrepareC3Dfile';
            app.PrepareC3DfileUIFigure.Interruptible = 'off';

            % Create UIAxesSteps
            app.UIAxesSteps = uiaxes(app.PrepareC3DfileUIFigure);
            title(app.UIAxesSteps, 'Steps')
            ylabel(app.UIAxesSteps, '\color{red}left \color{black} and \color{green} right \color{black} cycles')
            app.UIAxesSteps.YTick = [-1 -0.5 0 0.5 1];
            app.UIAxesSteps.YTickLabel = {'-1'; '-0.5'; '0'; '0.5'; '1'};
            app.UIAxesSteps.Position = [16 327 1490 158];

            % Create UIAxesGRFs
            app.UIAxesGRFs = uiaxes(app.PrepareC3DfileUIFigure);
            title(app.UIAxesGRFs, 'GRFs')
            xlabel(app.UIAxesGRFs, 'frame')
            ylabel(app.UIAxesGRFs, '[N]')
            app.UIAxesGRFs.YGrid = 'on';
            app.UIAxesGRFs.Position = [16 24 1490 304];

            % Create Selectc3dfileButton
            app.Selectc3dfileButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.Selectc3dfileButton.ButtonPushedFcn = createCallbackFcn(app, @Selectc3dfileButtonPushed, true);
            app.Selectc3dfileButton.Position = [39 733 134 22];
            app.Selectc3dfileButton.Text = 'Select c3d file';

            % Create PrepareC3DfileforOpenSimsimulationsLabel
            app.PrepareC3DfileforOpenSimsimulationsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.PrepareC3DfileforOpenSimsimulationsLabel.HorizontalAlignment = 'center';
            app.PrepareC3DfileforOpenSimsimulationsLabel.FontSize = 20;
            app.PrepareC3DfileforOpenSimsimulationsLabel.FontWeight = 'bold';
            app.PrepareC3DfileforOpenSimsimulationsLabel.Position = [1 848 1283 70];
            app.PrepareC3DfileforOpenSimsimulationsLabel.Text = {'Prepare C3D file for OpenSim simulations'; 'creating marker_experimental.trc / grf.mot / grf.xml / settings.mat and a copy of original c3d-file'};

            % Create SelectaC3DfileLabel
            app.SelectaC3DfileLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.SelectaC3DfileLabel.Position = [197 733 682 22];
            app.SelectaC3DfileLabel.Text = 'Select a C3D file';

            % Create detectforceplatesautomaticallyCheckBox
            app.detectforceplatesautomaticallyCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.detectforceplatesautomaticallyCheckBox.ValueChangedFcn = createCallbackFcn(app, @detectforceplatesautomaticallyCheckBoxValueChanged, true);
            app.detectforceplatesautomaticallyCheckBox.Text = 'detect forceplates automatically';
            app.detectforceplatesautomaticallyCheckBox.Position = [53 655 191 22];
            app.detectforceplatesautomaticallyCheckBox.Value = true;

            % Create filterGRFsCheckBox
            app.filterGRFsCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.filterGRFsCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.filterGRFsCheckBox.Text = 'filter GRFs';
            app.filterGRFsCheckBox.Position = [52 594 80 22];
            app.filterGRFsCheckBox.Value = true;

            % Create detectwalkingdirectionautomaticallybymarkerCheckBox
            app.detectwalkingdirectionautomaticallybymarkerCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.ValueChangedFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Text = 'detect walking direction automatically by marker';
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Position = [53 676 280 22];
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = true;

            % Create rotatearoundyaxisDropDownLabel
            app.rotatearoundyaxisDropDownLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.rotatearoundyaxisDropDownLabel.HorizontalAlignment = 'right';
            app.rotatearoundyaxisDropDownLabel.Enable = 'off';
            app.rotatearoundyaxisDropDownLabel.Position = [404 676 109 22];
            app.rotatearoundyaxisDropDownLabel.Text = 'rotate around y axis';

            % Create rotatearoundyaxisDropDown
            app.rotatearoundyaxisDropDown = uidropdown(app.PrepareC3DfileUIFigure);
            app.rotatearoundyaxisDropDown.Items = {'0°', '90°', '180°', '270°'};
            app.rotatearoundyaxisDropDown.Enable = 'off';
            app.rotatearoundyaxisDropDown.Position = [528 676 63 22];
            app.rotatearoundyaxisDropDown.Value = '0°';

            % Create CreatefilesButton
            app.CreatefilesButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.CreatefilesButton.ButtonPushedFcn = createCallbackFcn(app, @CreatefilesButtonPushed, true);
            app.CreatefilesButton.Enable = 'off';
            app.CreatefilesButton.Position = [983 532 139 34];
            app.CreatefilesButton.Text = 'Create files';

            % Create leftListBoxLabel
            app.leftListBoxLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.leftListBoxLabel.HorizontalAlignment = 'right';
            app.leftListBoxLabel.Position = [609 542 25 22];
            app.leftListBoxLabel.Text = 'left ';

            % Create leftListBox
            app.leftListBox = uilistbox(app.PrepareC3DfileUIFigure);
            app.leftListBox.Items = {};
            app.leftListBox.Multiselect = 'on';
            app.leftListBox.ValueChangedFcn = createCallbackFcn(app, @updateCycles, true);
            app.leftListBox.Position = [649 492 100 74];
            app.leftListBox.Value = {};

            % Create validfootstrikeeventsLabel
            app.validfootstrikeeventsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.validfootstrikeeventsLabel.FontWeight = 'bold';
            app.validfootstrikeeventsLabel.Position = [694 581 136 22];
            app.validfootstrikeeventsLabel.Text = 'valid foot strike events';

            % Create rightListBoxLabel
            app.rightListBoxLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.rightListBoxLabel.HorizontalAlignment = 'right';
            app.rightListBoxLabel.Position = [779 542 29 22];
            app.rightListBoxLabel.Text = 'right';

            % Create rightListBox
            app.rightListBox = uilistbox(app.PrepareC3DfileUIFigure);
            app.rightListBox.Items = {};
            app.rightListBox.Multiselect = 'on';
            app.rightListBox.ValueChangedFcn = createCallbackFcn(app, @updateCycles, true);
            app.rightListBox.Position = [823 492 100 74];
            app.rightListBox.Value = {};

            % Create openselectedfilewithdefaultprogramButton
            app.openselectedfilewithdefaultprogramButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.openselectedfilewithdefaultprogramButton.ButtonPushedFcn = createCallbackFcn(app, @openselectedfilewithdefaultprogramButtonPushed, true);
            app.openselectedfilewithdefaultprogramButton.Visible = 'off';
            app.openselectedfilewithdefaultprogramButton.Position = [39 707 222 22];
            app.openselectedfilewithdefaultprogramButton.Text = 'open selected file with default program';

            % Create LogTextAreaLabel
            app.LogTextAreaLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.LogTextAreaLabel.HorizontalAlignment = 'right';
            app.LogTextAreaLabel.Position = [1126 825 26 22];
            app.LogTextAreaLabel.Text = 'Log';

            % Create LogTextArea
            app.LogTextArea = uitextarea(app.PrepareC3DfileUIFigure);
            app.LogTextArea.Editable = 'off';
            app.LogTextArea.Position = [1167 484 338 365];

            % Create filterorderSpinnerLabel
            app.filterorderSpinnerLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.filterorderSpinnerLabel.HorizontalAlignment = 'right';
            app.filterorderSpinnerLabel.Position = [163 594 59 22];
            app.filterorderSpinnerLabel.Text = 'filter order';

            % Create filterorderSpinner
            app.filterorderSpinner = uispinner(app.PrepareC3DfileUIFigure);
            app.filterorderSpinner.Step = 2;
            app.filterorderSpinner.Limits = [2 10];
            app.filterorderSpinner.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.filterorderSpinner.Position = [233 594 53 22];
            app.filterorderSpinner.Value = 4;

            % Create cutofffrequencySpinnerLabel
            app.cutofffrequencySpinnerLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinnerLabel.HorizontalAlignment = 'right';
            app.cutofffrequencySpinnerLabel.Position = [300 594 91 22];
            app.cutofffrequencySpinnerLabel.Text = 'cutoff frequency';

            % Create cutofffrequencySpinner
            app.cutofffrequencySpinner = uispinner(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner.Limits = [4 100];
            app.cutofffrequencySpinner.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.cutofffrequencySpinner.Position = [403 594 53 22];
            app.cutofffrequencySpinner.Value = 20;

            % Create detectWalkingByMarker
            app.detectWalkingByMarker = uieditfield(app.PrepareC3DfileUIFigure, 'text');
            app.detectWalkingByMarker.ValueChangedFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectWalkingByMarker.ValueChangingFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectWalkingByMarker.Position = [332 676 73 22];
            app.detectWalkingByMarker.Value = 'LTOE';

            % Create detectFPbyMarker
            app.detectFPbyMarker = uieditfield(app.PrepareC3DfileUIFigure, 'text');
            app.detectFPbyMarker.ValueChangedFcn = createCallbackFcn(app, @detectforceplatesautomaticallyCheckBoxValueChanged, true);
            app.detectFPbyMarker.Position = [253 655 73 22];
            app.detectFPbyMarker.Value = 'LTOE';

            % Create v20Label
            app.v20Label = uilabel(app.PrepareC3DfileUIFigure);
            app.v20Label.Position = [1467 872 28 22];
            app.v20Label.Text = 'v2.0';

            % Create forceplatecontactsLabel
            app.forceplatecontactsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.forceplatecontactsLabel.FontWeight = 'bold';
            app.forceplatecontactsLabel.Visible = 'off';
            app.forceplatecontactsLabel.Position = [698 581 136 22];
            app.forceplatecontactsLabel.Text = 'force plate contacts';

            % Create filtermarkersCheckBox
            app.filtermarkersCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.filtermarkersCheckBox.Text = 'filter markers';
            app.filtermarkersCheckBox.Position = [52 629 92 22];

            % Create filterorderSpinnerLabel_Marker
            app.filterorderSpinnerLabel_Marker = uilabel(app.PrepareC3DfileUIFigure);
            app.filterorderSpinnerLabel_Marker.HorizontalAlignment = 'right';
            app.filterorderSpinnerLabel_Marker.Position = [163 629 59 22];
            app.filterorderSpinnerLabel_Marker.Text = 'filter order';

            % Create filterorderSpinner_Marker
            app.filterorderSpinner_Marker = uispinner(app.PrepareC3DfileUIFigure);
            app.filterorderSpinner_Marker.Step = 2;
            app.filterorderSpinner_Marker.Limits = [2 10];
            app.filterorderSpinner_Marker.Position = [233 629 53 22];
            app.filterorderSpinner_Marker.Value = 2;

            % Create cutofffrequencySpinner_2Label
            app.cutofffrequencySpinner_2Label = uilabel(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner_2Label.HorizontalAlignment = 'right';
            app.cutofffrequencySpinner_2Label.Position = [300 629 91 22];
            app.cutofffrequencySpinner_2Label.Text = 'cutoff frequency';

            % Create cutofffrequencySpinner_Marker
            app.cutofffrequencySpinner_Marker = uispinner(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner_Marker.Limits = [4 30];
            app.cutofffrequencySpinner_Marker.Position = [403 629 53 22];
            app.cutofffrequencySpinner_Marker.Value = 5;

            % Create ClimbingPanel
            app.ClimbingPanel = uipanel(app.PrepareC3DfileUIFigure);
            app.ClimbingPanel.Title = 'Climbing';
            app.ClimbingPanel.FontWeight = 'bold';
            app.ClimbingPanel.Position = [804 787 288 60];

            % Create ClimbingCheckBox
            app.ClimbingCheckBox = uicheckbox(app.ClimbingPanel);
            app.ClimbingCheckBox.Text = 'Climbing';
            app.ClimbingCheckBox.Position = [8 8 69 22];

            % Create invertforcesensorsCheckBox
            app.invertforcesensorsCheckBox = uicheckbox(app.ClimbingPanel);
            app.invertforcesensorsCheckBox.Text = 'invert force sensors';
            app.invertforcesensorsCheckBox.Position = [91 8 126 22];

            % Create mirrorCheckBox
            app.mirrorCheckBox = uicheckbox(app.ClimbingPanel);
            app.mirrorCheckBox.Tooltip = {'EMG stays the same!'};
            app.mirrorCheckBox.Text = 'mirror';
            app.mirrorCheckBox.Position = [222 8 53 22];

            % Create SetbeforeopeningC3DPanel
            app.SetbeforeopeningC3DPanel = uipanel(app.PrepareC3DfileUIFigure);
            app.SetbeforeopeningC3DPanel.Title = 'Set before opening C3D';
            app.SetbeforeopeningC3DPanel.FontWeight = 'bold';
            app.SetbeforeopeningC3DPanel.Position = [41 766 749 81];

            % Create removeblanksetcfromfilenameCheckBox
            app.removeblanksetcfromfilenameCheckBox = uicheckbox(app.SetbeforeopeningC3DPanel);
            app.removeblanksetcfromfilenameCheckBox.ValueChangedFcn = createCallbackFcn(app, @removeblanksetcfromfilenameCheckBoxValueChanged, true);
            app.removeblanksetcfromfilenameCheckBox.Text = 'remove blanks etc from filename';
            app.removeblanksetcfromfilenameCheckBox.Position = [8 29 196 22];
            app.removeblanksetcfromfilenameCheckBox.Value = true;

            % Create deleteoriginalfileafterrenamingCheckBox
            app.deleteoriginalfileafterrenamingCheckBox = uicheckbox(app.SetbeforeopeningC3DPanel);
            app.deleteoriginalfileafterrenamingCheckBox.Text = 'delete original file after renaming';
            app.deleteoriginalfileafterrenamingCheckBox.Position = [213 29 196 22];
            app.deleteoriginalfileafterrenamingCheckBox.Value = true;

            % Create AMTITandemTreadmillCheckBox
            app.AMTITandemTreadmillCheckBox = uicheckbox(app.SetbeforeopeningC3DPanel);
            app.AMTITandemTreadmillCheckBox.ValueChangedFcn = createCallbackFcn(app, @AMTITandemTreadmillCheckBoxValueChanged, true);
            app.AMTITandemTreadmillCheckBox.Tooltip = {'EMG stays the same!'};
            app.AMTITandemTreadmillCheckBox.Text = 'AMTI Tandem Treadmill';
            app.AMTITandemTreadmillCheckBox.Position = [418 30 148 22];

            % Create Ignorec3deventsCheckBox
            app.Ignorec3deventsCheckBox = uicheckbox(app.SetbeforeopeningC3DPanel);
            app.Ignorec3deventsCheckBox.ValueChangedFcn = createCallbackFcn(app, @Ignorec3deventsCheckBoxValueChanged, true);
            app.Ignorec3deventsCheckBox.Text = 'Ignore c3d events and set forceplate contacts manually';
            app.Ignorec3deventsCheckBox.Position = [8 5 320 22];

            % Create MovementTypeDropDownLabel
            app.MovementTypeDropDownLabel = uilabel(app.SetbeforeopeningC3DPanel);
            app.MovementTypeDropDownLabel.HorizontalAlignment = 'right';
            app.MovementTypeDropDownLabel.Position = [397 6 90 22];
            app.MovementTypeDropDownLabel.Text = 'Movement Type';

            % Create MovementTypeDropDown
            app.MovementTypeDropDown = uidropdown(app.SetbeforeopeningC3DPanel);
            app.MovementTypeDropDown.Items = {'Gait', 'Run', 'Counter Movement Jump', 'Stand-Sit-Stand', 'SideStep'};
            app.MovementTypeDropDown.Position = [502 6 238 22];
            app.MovementTypeDropDown.Value = 'Gait';

            % Create automaticMovementtypeCheckBox
            app.automaticMovementtypeCheckBox = uicheckbox(app.SetbeforeopeningC3DPanel);
            app.automaticMovementtypeCheckBox.Text = 'automatic Movement type';
            app.automaticMovementtypeCheckBox.Position = [582 31 159 22];
            app.automaticMovementtypeCheckBox.Value = true;

            % Create TimingSettingsPanel
            app.TimingSettingsPanel = uipanel(app.PrepareC3DfileUIFigure);
            app.TimingSettingsPanel.Title = 'Timing Settings';
            app.TimingSettingsPanel.FontWeight = 'bold';
            app.TimingSettingsPanel.Position = [939 587 214 133];

            % Create exporteachstepseparatelyCheckBox
            app.exporteachstepseparatelyCheckBox = uicheckbox(app.TimingSettingsPanel);
            app.exporteachstepseparatelyCheckBox.Text = 'export each step separately';
            app.exporteachstepseparatelyCheckBox.Position = [10 81 169 22];
            app.exporteachstepseparatelyCheckBox.Value = true;

            % Create croptovalidfootstrikesbufferCheckBox
            app.croptovalidfootstrikesbufferCheckBox = uicheckbox(app.TimingSettingsPanel);
            app.croptovalidfootstrikesbufferCheckBox.Text = 'crop to valid foot strikes + buffer';
            app.croptovalidfootstrikesbufferCheckBox.Position = [11 57 193 22];
            app.croptovalidfootstrikesbufferCheckBox.Value = true;

            % Create BufferStartsLabel
            app.BufferStartsLabel = uilabel(app.TimingSettingsPanel);
            app.BufferStartsLabel.HorizontalAlignment = 'right';
            app.BufferStartsLabel.Position = [8 30 81 22];
            app.BufferStartsLabel.Text = 'Buffer Start [s]';

            % Create BufferStartsSpinner
            app.BufferStartsSpinner = uispinner(app.TimingSettingsPanel);
            app.BufferStartsSpinner.Step = 0.1;
            app.BufferStartsSpinner.Position = [154 30 50 22];
            app.BufferStartsSpinner.Value = 0.2;

            % Create BufferEndsSpinnerLabel
            app.BufferEndsSpinnerLabel = uilabel(app.TimingSettingsPanel);
            app.BufferEndsSpinnerLabel.HorizontalAlignment = 'right';
            app.BufferEndsSpinnerLabel.Position = [8 9 77 22];
            app.BufferEndsSpinnerLabel.Text = 'Buffer End [s]';

            % Create BufferEndsSpinner
            app.BufferEndsSpinner = uispinner(app.TimingSettingsPanel);
            app.BufferEndsSpinner.Step = 0.1;
            app.BufferEndsSpinner.Position = [154 9 50 22];
            app.BufferEndsSpinner.Value = 0.2;

            % Create EMGSettingsPanel
            app.EMGSettingsPanel = uipanel(app.PrepareC3DfileUIFigure);
            app.EMGSettingsPanel.Title = 'EMG Settings';
            app.EMGSettingsPanel.FontWeight = 'bold';
            app.EMGSettingsPanel.Position = [42 482 530 106];

            % Create filterEMGCheckBox
            app.filterEMGCheckBox = uicheckbox(app.EMGSettingsPanel);
            app.filterEMGCheckBox.ValueChangedFcn = createCallbackFcn(app, @filterEMGCheckBoxValueChanged, true);
            app.filterEMGCheckBox.Text = 'filter EMG';
            app.filterEMGCheckBox.Position = [8 54 176 22];
            app.filterEMGCheckBox.Value = true;

            % Create SelectexistingCSVButton
            app.SelectexistingCSVButton = uibutton(app.EMGSettingsPanel, 'push');
            app.SelectexistingCSVButton.ButtonPushedFcn = createCallbackFcn(app, @SelectexistingCSVButtonPushed, true);
            app.SelectexistingCSVButton.Position = [308 5 121 22];
            app.SelectexistingCSVButton.Text = 'Select existing CSV';

            % Create createCSVandopenButton
            app.createCSVandopenButton = uibutton(app.EMGSettingsPanel, 'push');
            app.createCSVandopenButton.ButtonPushedFcn = createCallbackFcn(app, @createCSVandopenButtonPushed, true);
            app.createCSVandopenButton.Position = [154 5 130 22];
            app.createCSVandopenButton.Text = 'create CSV and open';

            % Create renameEMGlabelsCheckBox
            app.renameEMGlabelsCheckBox = uicheckbox(app.EMGSettingsPanel);
            app.renameEMGlabelsCheckBox.Text = 'rename EMG labels';
            app.renameEMGlabelsCheckBox.Position = [8 5 128 22];

            % Create normalizeallmagnitudesto1CheckBox
            app.normalizeallmagnitudesto1CheckBox = uicheckbox(app.EMGSettingsPanel);
            app.normalizeallmagnitudesto1CheckBox.Text = 'normalize all magnitudes to 1';
            app.normalizeallmagnitudesto1CheckBox.Position = [114 55 177 22];

            % Create exportrawEMGCheckBox
            app.exportrawEMGCheckBox = uicheckbox(app.EMGSettingsPanel);
            app.exportrawEMGCheckBox.Text = 'export raw EMG';
            app.exportrawEMGCheckBox.Position = [7 30 108 22];

            % Create EMGfactorSpinnerLabel
            app.EMGfactorSpinnerLabel = uilabel(app.EMGSettingsPanel);
            app.EMGfactorSpinnerLabel.HorizontalAlignment = 'right';
            app.EMGfactorSpinnerLabel.Position = [319 56 66 22];
            app.EMGfactorSpinnerLabel.Text = 'EMG factor';

            % Create EMGfactorSpinner
            app.EMGfactorSpinner = uispinner(app.EMGSettingsPanel);
            app.EMGfactorSpinner.Limits = [1 10000];
            app.EMGfactorSpinner.Tooltip = {'Use this factor if output signal has steps.'; 'Of course, you need to have the same factor for all trials of the same participant!'};
            app.EMGfactorSpinner.Position = [396 56 54 22];
            app.EMGfactorSpinner.Value = 1;

            % Create Image
            app.Image = uiimage(app.PrepareC3DfileUIFigure);
            app.Image.Position = [622 614 289 105];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'v1_transparent_upscaled.png');

            % Show the figure after all components are created
            app.PrepareC3DfileUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = processC3D_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PrepareC3DfileUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PrepareC3DfileUIFigure)
        end
    end
end