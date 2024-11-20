classdef processC3D_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PrepareC3DfileUIFigure          matlab.ui.Figure
        mirrorCheckBox                  matlab.ui.control.CheckBox
        invertforcesensorsCheckBox      matlab.ui.control.CheckBox
        normalizeEMGCheckBox            matlab.ui.control.CheckBox
        cutofffrequencySpinner_Marker   matlab.ui.control.Spinner
        cutofffrequencySpinner_2Label   matlab.ui.control.Label
        filterorderSpinner_Marker       matlab.ui.control.Spinner
        filterorderSpinnerLabel_Marker  matlab.ui.control.Label
        filtermarkersCheckBox           matlab.ui.control.CheckBox
        ClimbingCheckBox                matlab.ui.control.CheckBox
        Ignorec3deventsCheckBox         matlab.ui.control.CheckBox
        forceplatecontactsLabel         matlab.ui.control.Label
        EMGfactorSpinner                matlab.ui.control.Spinner
        EMGfactorSpinnerLabel           matlab.ui.control.Label
        renameEMGlabelsCheckBox         matlab.ui.control.CheckBox
        createCSVandopenButton          matlab.ui.control.Button
        SelectexistingCSVButton         matlab.ui.control.Button
        filterEMGandexporttostoCheckBox  matlab.ui.control.CheckBox
        v16Label                        matlab.ui.control.Label
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
        setGRFstozerooutsideregionoffootcontactCheckBox  matlab.ui.control.CheckBox
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
        grf_forces_zero_leveled % Description
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
            drawnow; % Force text area to update
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
                disp('clicked')
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

        %         function fillCallback(app, evnt)
        %             disp('clicked')
        %         end
        %                 disp('clicked')
        %                 src.Vertices(1, 1);
        %                 src.Vertices(4, 1);
        %                 if src.Vertices(2, 2) > 0 % left
        %
        %                 else % right
        %
        %                 end
        %                 leftStepForcePlateAssignment
        %                 if evnt.Button == 1
        %                     app.states{s} = 5;
        %                 end
        %             end
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
            if app.filterGRFsCheckBox.Value && app.c3d.getNumForces() > 0
                app.grf_forces_adjusted = [];
                Fs = app.forceFrequency;  % Sampling Frequency
                N  = app.filterorderSpinner.Value;  % Filter Order
                Fc = app.cutofffrequencySpinner.Value;  % Cutoff Frequency
                [filter_b, filter_a] = butter(N/2,Fc/(Fs/2), 'low'); % divide order by 2 because filtfilt doubles it again
                forcesFields = fieldnames(app.grf_forces_zero_leveled);

                forcesFields = sort(forcesFields);
                for i = 1 : numel(forcesFields)
                    if char(forcesFields(i)) ~= "time" && ~contains(forcesFields(i), '_p')
                        fieldForThreshold = strrep(forcesFields(i), '_moment', '_force');
                        fieldForThreshold{1}(end) = 'y';
                        fieldForThreshold{1}(end-1) = 'v';
                        app.grf_forces_adjusted.(char(forcesFields(i))) = filtfilt(filter_b, filter_a, app.grf_forces_zero_leveled.(char(forcesFields(i))));
                        if app.setGRFstozerooutsideregionoffootcontactCheckBox.Value == 1
                            firstIndexOverThreshold = find(abs(app.grf_forces_zero_leveled.(fieldForThreshold{1})) > threshold, 1);
                            lastIndexOverThreshold = find(abs(app.grf_forces_zero_leveled.(fieldForThreshold{1})(firstIndexOverThreshold + 100 : end)) < threshold, 1) + firstIndexOverThreshold + 99;
                            app.grf_forces_adjusted.(char(forcesFields(i)))(1 : firstIndexOverThreshold-1) = 0;
                            app.grf_forces_adjusted.(char(forcesFields(i)))(lastIndexOverThreshold : end) = 0;
                        end
                    else % other filtering for point of application 
                        if contains(forcesFields(i), '_p') && app.setGRFstozerooutsideregionoffootcontactCheckBox.Value == 1
                            fieldForThreshold = strrep(forcesFields(i), '_moment', '_force');
                            fieldForThreshold{1}(end) = 'y';
                            fieldForThreshold{1}(end-1) = 'v';
                            firstIndexOverThreshold = find(abs(app.grf_forces_zero_leveled.(fieldForThreshold{1})) > threshold, 1);

                            lastIndexOverThreshold = find(abs(app.grf_forces_zero_leveled.(fieldForThreshold{1})(firstIndexOverThreshold + 100 : end)) < threshold, 1) + firstIndexOverThreshold + 100;
                            tmp_maxTimeThresholdStart = min(timeThresholdFPs, firstIndexOverThreshold);
                            tmp_step = app.grf_forces_zero_leveled.(char(forcesFields(i)))(firstIndexOverThreshold + tmp_maxTimeThresholdStart : lastIndexOverThreshold - timeThresholdFPs);
                            tmp_step_filtered = filtfilt(filter_b, filter_a, tmp_step);
                            if size(tmp_step_filtered, 1) > 1
                                tmp_interpolated = interp1(1 : size(tmp_step_filtered, 1), tmp_step_filtered, linspace(1, size(tmp_step_filtered, 1), lastIndexOverThreshold - firstIndexOverThreshold + tmp_maxTimeThresholdStart + timeThresholdFPs));
                                if lastIndexOverThreshold >= size(app.grf_forces_zero_leveled.(fieldForThreshold{1}), 1)
                                    tmp_interpolated = tmp_interpolated(1 : end - 30);
                                end
                            app.grf_forces_adjusted.(char(forcesFields(i))) = [zeros(1, firstIndexOverThreshold - tmp_maxTimeThresholdStart), tmp_interpolated, zeros(1, length(app.grf_forces_zero_leveled.(char(forcesFields(i)))) - lastIndexOverThreshold - timeThresholdFPs)]';
                            else
                                app.grf_forces_adjusted.(char(forcesFields(i))) = app.grf_forces_zero_leveled.(char(forcesFields(i)));
                            end
                        else
                            % do not filter time
                            app.grf_forces_adjusted.(char(forcesFields(i))) = app.grf_forces_zero_leveled.(char(forcesFields(i)));
                        end
                    end
                end
            else
                app.grf_forces_adjusted = app.grf_forces_zero_leveled;
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
                    app.filterEMGandexporttostoCheckBox.Value = 0;
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
                    if app.normalizeEMGCheckBox.Value == 1
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

                if app.Ignorec3deventsCheckBox.Value == 0
                    app.validfootstrikeeventsLabel.Visible = "on";
                    app.forceplatecontactsLabel.Visible = "off";
                    if isfield(c3devents, 'Left_Foot_Strike')
                        app.leftFootStrikeList = sprintfc('%.0f', c3devents.Left_Foot_Strike * app.frequency_);
                        app.leftListBox.Items = app.leftFootStrikeList(1 : end-1);
                        app.leftListBox.Value = app.leftFootStrikeList(1 : end-1);
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
                        app.rightListBox.Value = app.rightFootStrikeList(1 : end-1);
                    else
                        app.rightListBox.Items = {};
                        app.rightFootStrikeList = {};
                    end
                    if isfield(c3devents, 'Right_Foot_Off')
                        app.rightFootOffList = floor(c3devents.Right_Foot_Off * app.frequency_);
                    else
                        app.rightFootOffList = [];
                    end                
                    if isempty(app.rightListBox.Items) && isempty(app.leftListBox.Items)
                        answer = 'Yes'; % MFquestdlg([0.4, 0.4], 'No events in the c3d file! Therefore, start and end of the trial as well as the foot contact with the force plates cannot be calculated automatically! Do you want to set the force plates manually?', 'Define manually?', 'Yes', 'No, I''ll experimentally check what the result of this is or select another file', 'Yes');
                        if ~strcmp(answer, 'Yes')
                            return;
                        end
                        if strcmp(answer, 'Yes')
                            app.Ignorec3deventsCheckBox.Value = 1;
                            app.detectforceplatesautomaticallyCheckBox.Value = 0;
                            app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;
                            app.validfootstrikeeventsLabel.Visible = "off";
                            app.forceplatecontactsLabel.Visible = "on";
                        end
                    end
                else
%                     % not good if you open another trial with the same
%                     % settings...
%                     app.rightListBox.Items = {};
%                     app.rightFootStrikeList = {};
%                     app.leftListBox.Items = {};
%                     app.leftFootStrikeList = {};
                end

                drawnow;
                figure(app.PrepareC3DfileUIFigure);

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

                app.c3d = osimC3D(c3dpath, 1);
                try
                    app.forceFrequency = app.c3d.getRate_force();
                catch ME
                    app.forceFrequency = 1000;
                end
                app.c3d.rotateData('x', -90);
                app.c3d.convertMillimeters2Meters();
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
                d.Message = 'process and filter forces...';
                if app.c3d.getNumForces() ~= 0
                    setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged(app, 0);
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

                if app.filterEMGandexporttostoCheckBox.Value
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


            drawnow;
            figure(app.PrepareC3DfileUIFigure);

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

                        distanceTravelled = app.markers.(markerToTrack)(end, :) - app.markers.(markerToTrack)(1, :);
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

                % create folder with the same name as c3d file
                if ~exist(output_folder, 'dir')
                    mkdir(output_folder)
                end

                d.Value = d.Value + 0.1;
                d.Message = 'rotating data ...';
                if ~app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value
                    app.rotationAngle = str2double(strrep(app.rotatearoundyaxisDropDown.Value, '°', ''));
                end
                c3d_tmp = app.c3d;
                c3d_tmp.rotateData('y', app.rotationAngle);
                d.Value = d.Value + 0.1;
                d.Message = 'writing marker_experimental.trc file with marker locations...';
                c3d_tmp.writeTRC(fullfile(output_folder, 'marker_experimental.trc'));

                if app.filtermarkersCheckBox.Value
                    disp('filtering marker data...')
                    d.Message = 'filtering marker data...';
                    markerData = load_marker_trc(fullfile(output_folder, 'marker_experimental.trc'));
                    markerNames = fieldnames(markerData);

                    Fs = app.frequency_;  % Sampling Frequency
                    N  = app.filterorderSpinner_Marker.Value;  % Filter Order
                    Fc = app.cutofffrequencySpinner_Marker.Value;  % Cutoff Frequency
                    [filter_b, filter_a] = butter(N/2,Fc/(Fs/2), 'low'); % divide order by 2 because filtfilt doubles it again

                    for i = 1 : numel(markerNames)
                        if strcmp(markerNames{i}(end-1:end), '_X') || strcmp(markerNames{i}(end-1:end), '_Y') || strcmp(markerNames{i}(end-1:end), '_Z')
                            try
                                markerData.(char(markerNames(i))) = filtfilt(filter_b, filter_a, cell2mat(markerData.(char(markerNames(i)))));
                            catch
                                markerData = rmfield(markerData, (char(markerNames(i))));
                            end
                        else
                            try
                                markerData.(char(markerNames(i))) = cell2mat(markerData.(char(markerNames(i))));
                            catch
                                markerData = rmfield(markerData, (char(markerNames(i))));
                            end
                        end
                    end

                    text = fileread(fullfile(output_folder, 'marker_experimental.trc'));
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
                            filteredTraj = filtfilt(filter_b, filter_a, traj);
                            for row = 7 : numel(data)-1 % skip 6 header rows
                                data{row}{col} = filteredTraj(row-6);
                                % traj(row-6) = cell2mat();
                            end
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


                    status = movefile(fullfile(output_folder, 'marker_experimental.trc'), fullfile(output_folder, 'marker_experimental_unfiltered.trc'));
                    if status
                        fid = fopen(fullfile(output_folder, 'marker_experimental.trc'),'wt');
                    else
                        fid = fopen(fullfile(output_folder, 'marker_experimental_filtered.trc'),'wt');
                    end
                    fprintf(fid, finalStr);
                    fclose(fid);
                end


                if app.mirrorCheckBox.Value
                    mirrorMarkerData(app, fullfile(output_folder, 'marker_experimental.trc'));
                end

                nForces = app.c3d.getNumForces();

                d.Value = d.Value + 0.1;
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

                    setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged(app, 0);
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

                    d.Value = d.Value + 0.1;
                    d.Message = 'writing forces to grf.mot file...';

                    write_sto_file(app.grf_forces_adjusted, grfFileName);
                end

                grforces = xml_read('GRF_file_all.xml');
                grforces_generated = xml_read('GRF_file_empty.xml');
                counter = 1;
                if app.Ignorec3deventsCheckBox.Value == 0
                    if ~app.mirrorCheckBox.Value
                        for i = 1 : size(app.rightStepForcePlateAssignment, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(app.rightStepForcePlateAssignment(i, 2));
                            counter = counter + 1;
                        end
                        for i = 1 : size(app.leftStepForcePlateAssignment, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(9 + app.leftStepForcePlateAssignment(i, 2));
                            counter = counter + 1;
                        end
                    else
                        for i = 1 : size(app.rightStepForcePlateAssignment, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(9 + app.rightStepForcePlateAssignment(i, 2));
                            counter = counter + 1;
                        end
                        for i = 1 : size(app.leftStepForcePlateAssignment, 1)
                            grforces_generated.ExternalLoads.objects.ExternalForce(counter) = grforces.ExternalLoads.objects.ExternalForce(app.leftStepForcePlateAssignment(i, 2));
                            counter = counter + 1;
                        end
                    end
                    Pref = struct;
                    Pref.StructItem = false;
    
                    d.Value = d.Value + 0.1;
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
    
                    d.Value = d.Value + 0.1;
                    d.Message = 'writing GRF.XML ...';
                    xml_write(fullfile(output_folder, 'GRF.xml'), grforces_generated, 'OpenSimDocument', Pref);
                end

                d.Value = d.Value + 0.1;
                cycle = struct;
                if app.Ignorec3deventsCheckBox.Value == 0
                    d.Message = 'preparing information about cycles ...';
                    cycle = struct;
                    if ~isempty(app.leftListBox.Value)
                        cycle.left = struct;
                        cycle.left.start = str2double(app.leftListBox.Value) - app.firstFrame_;
                        footStrikes = str2double(app.leftFootStrikeList) - app.firstFrame_;
                        for i = 1 : size(cycle.left.start, 2)
                            cycle.left.end(i) = footStrikes( find( footStrikes > cycle.left.start(i), 1 ) );
                            if ~isempty(app.leftFootOffList)
                                cycle.left.footOff(i) = app.leftFootOffList( find( app.leftFootOffList  > cycle.left.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                            end
                        end
                    end
                    if ~isempty(app.rightListBox.Value)
                        cycle.right = struct;
                        cycle.right.start = str2double(app.rightListBox.Value) - app.firstFrame_;
                        footStrikes = str2double(app.rightFootStrikeList) - app.firstFrame_;
                        for i = 1 : size(cycle.right.start, 2)
                            cycle.right.end(i) = footStrikes( find( footStrikes > cycle.right.start(i), 1 ) );
                            if ~isempty(app.rightFootOffList)
                                cycle.right.footOff(i) = app.rightFootOffList( find( app.rightFootOffList > cycle.right.start(i) + app.firstFrame_, 1 ) ) - app.firstFrame_;
                            end
                        end
                    end

                    cycle_tmp = cycle;
                    cycle = struct;
                    if app.mirrorCheckBox.Value
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
                duration = (btkGetPointFrameNumber(app.acq) - 1) / frequency;

                d.Value = d.Value + 0.1;
                d.Message = 'write settings.mat and copy c3d file ...';
                save(fullfile(output_folder, 'settings.mat'), 'cycle', 'firstFrame', 'forceFrequency', 'frequency', 'duration', '-mat');
                copyfile(fullfile(app.path, app.c3dFileName), fullfile(output_folder, 'c3dfile.c3d'));

                d.Value = d.Value + 0.1;
                invalidFieldNames = 0;
                if app.filterEMGandexporttostoCheckBox.Value
                    if app.renameEMGlabelsCheckBox.Value
                        if ~isempty(app.emgLabelCSVPath)
                            exportEMG = struct;
                            channels = fieldnames(app.EMG_lowPass);
                            channels = sort(channels);
                            channels = channels(~strcmp(channels, 'time'));
                            emgLabels = readtable(app.emgLabelCSVPath);
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
                    write_sto_file(exportEMG, fullfile(output_folder, 'EMG_filtered.sto'))
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

                    write_sto_file(forces, fullfile(output_folder, 'forces.mot'));

                    copyfile('climbing_forces.xml', fullfile(output_folder, 'GRF.xml'));
                end

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
                    %                 detectforceplatesautomaticallyCheckBoxValueChanged(app, 0);
                    %                 setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged(app, 0);


                    selectedFrames = cellfun(@str2double, app.rightListBox.Value);
                    app.rightStepForcePlateAssignment = app.rightStepForcePlateAssignmentAll;
                    isStillSelected = ismember(app.rightStepForcePlateAssignment(:, 1), selectedFrames);
                    app.rightStepForcePlateAssignment = app.rightStepForcePlateAssignment(isStillSelected, :);

                    selectedFrames = cellfun(@str2double, app.leftListBox.Value);
                    app.leftStepForcePlateAssignment = app.leftStepForcePlateAssignmentAll;
                    isStillSelected = ismember(app.leftStepForcePlateAssignment(:, 1), selectedFrames);
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
                        plot(app.UIAxesGRFs, app.grf_forces.(['ground_force_' num2str(i) '_vy']));
                    end

                    filterGRFs(app)

                    verticalForceFields = fieldnames(app.grf_forces_adjusted);
                    verticalForceFields = verticalForceFields(contains(verticalForceFields, '_vy'));
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
            end
        end

        % Value changed function: detectFPbyMarker, 
        % ...and 1 other component
        function detectforceplatesautomaticallyCheckBoxValueChanged(app, event)

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
                            % find closest forceplate to heel marker
                            heelLocation = temp_markers.(app.detectFPbyMarker.Value)(str2double(app.leftListBox.Value{i}) - app.firstFrame_ + 20, :) * 1000;
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
                        for i = 1 : size(app.rightListBox.Value, 2)
                            if isfield(temp_markers, (strrep(app.detectFPbyMarker.Value, 'L', 'R')))
                                rightMarkerName = (strrep(app.detectFPbyMarker.Value, 'L', 'R'));
                            else
                                rightMarkerName = app.detectFPbyMarker.Value;
                            end
                            % find closest forceplate to heel marker
                            heelLocation = temp_markers.(rightMarkerName)(str2double(app.rightListBox.Value{i}) - app.firstFrame_ + 20, :) * 1000;
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
                    else
                        app.detectFPbyMarker.Value = 'SET THIS!';
                        app.detectFPbyMarker.BackgroundColor = 'red';
                        app.detectFPbyMarker.FontColor = 'white';
                        app.detectforceplatesautomaticallyCheckBox.Value = 0;
                    end

                    app.c3d.rotateData('x', -90);

                    app.rightStepForcePlateAssignmentAll = app.rightStepForcePlateAssignment;
                    app.leftStepForcePlateAssignmentAll = app.leftStepForcePlateAssignment;
                    updateCyclePlot(app);
                end
                app.setGRFstozerooutsideregionoffootcontactCheckBox.Enable = "on";
            else
                app.leftStepForcePlateAssignment = [];
                app.rightStepForcePlateAssignment = [];
                app.leftStepForcePlateAssignmentAll = [];
                app.rightStepForcePlateAssignmentAll = [];
                updateCyclePlot(app);
                app.setGRFstozerooutsideregionoffootcontactCheckBox.Enable = "off";
                app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;

                if ~isempty(app.c3d) && app.c3d.getNumForces() ~= 0
                    setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged(app, 0);
                end
            end
        end

        % Value changed function: 
        % setGRFstozerooutsideregionoffootcontactCheckBox
        function setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged(app, event)
            if app.fileLoaded
                app.grf_forces_zero_leveled = [];
                if app.setGRFstozerooutsideregionoffootcontactCheckBox.Value

                    forcesFields = fieldnames(app.grf_forces);
                    footStrikes = str2double(app.leftListBox.Value);
                    for i = 1 : size(app.leftStepForcePlateAssignment, 1)
                        forcesFieldsOfThisFP = forcesFields(contains(forcesFields, num2str(app.leftStepForcePlateAssignment(i, 2))));
                        forcesFieldsOfThisFP = forcesFieldsOfThisFP(~contains(forcesFieldsOfThisFP, '_p'));
                        for j = 1 : numel(forcesFieldsOfThisFP)
                            nextStrike = footStrikes( find( footStrikes > app.leftStepForcePlateAssignment(i, 1), 1 ) );
                            nextStrike = nextStrike / app.frequency_ * app.forceFrequency - app.forceOffset;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) = app.grf_forces.(char(forcesFieldsOfThisFP(j)));
                            footStrikeIndex = app.leftStepForcePlateAssignment(i, 1) / app.frequency_ * app.forceFrequency - app.forceOffset;
                            if footStrikeIndex > 11
                                meanOffset = 10;
                            else
                                meanOffset = footStrikeIndex - 1;
                            end
                            meanVal = mean(forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(1 : footStrikeIndex - meanOffset));
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) = forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) - meanVal;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(1 : footStrikeIndex) = 0;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(nextStrike : end) = 0;
                            %                             if contains(forcesFieldsOfThisFP(j), '_vy')
                            %                                 forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) < 0) = 0;
                            %                             end
                        end
                        pointFieldsOfThisFP = forcesFields(contains(forcesFields, num2str(app.leftStepForcePlateAssignment(i, 2))));
                        pointFieldsOfThisFP = pointFieldsOfThisFP(contains(pointFieldsOfThisFP, '_p'));
                        for j = 1 : numel(pointFieldsOfThisFP)
                            forces_zero_leveled.(char(pointFieldsOfThisFP(j))) = app.grf_forces.(char(pointFieldsOfThisFP(j)));
                        end
                    end
                    footStrikes = str2double(app.rightListBox.Value);
                    for i = 1 : size(app.rightStepForcePlateAssignment, 1)
                        forcesFieldsOfThisFP = forcesFields(contains(forcesFields, num2str(app.rightStepForcePlateAssignment(i, 2))));
                        forcesFieldsOfThisFP = forcesFieldsOfThisFP(~contains(forcesFieldsOfThisFP, '_p'));
                        for j = 1 : numel(forcesFieldsOfThisFP)
                            nextStrike = footStrikes( find( footStrikes > app.rightStepForcePlateAssignment(i, 1), 1 ) );
                            nextStrike = nextStrike / app.frequency_ * app.forceFrequency - app.forceOffset;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) = app.grf_forces.(char(forcesFieldsOfThisFP(j)));
                            footStrikeIndex = app.rightStepForcePlateAssignment(i, 1) / app.frequency_ * app.forceFrequency - app.forceOffset;
                            if footStrikeIndex > 11
                                meanOffset = 10;
                            else
                                meanOffset = footStrikeIndex - 1;
                            end
                            meanVal = mean(forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(1 : footStrikeIndex - meanOffset));
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) = forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) - meanVal;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(1 : footStrikeIndex) = 0;
                            forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(nextStrike : end) = 0;
                            %                             if contains(forcesFieldsOfThisFP(j), '_vy')
                            %                                 forces_zero_leveled.(char(forcesFieldsOfThisFP(j)))(forces_zero_leveled.(char(forcesFieldsOfThisFP(j))) < 0) = 0;
                            %                             end
                        end
                        pointFieldsOfThisFP = forcesFields(contains(forcesFields, num2str(app.rightStepForcePlateAssignment(i, 2))));
                        pointFieldsOfThisFP = pointFieldsOfThisFP(contains(pointFieldsOfThisFP, '_p'));
                        for j = 1 : numel(pointFieldsOfThisFP)
                            forces_zero_leveled.(char(pointFieldsOfThisFP(j))) = app.grf_forces.(char(pointFieldsOfThisFP(j)));
                        end
                    end
                    forces_zero_leveled.time = app.grf_forces.time;
                    app.grf_forces_zero_leveled = forces_zero_leveled;
                else
                    app.grf_forces_zero_leveled = app.grf_forces;
                end
                updateGRF_Plot(app, 0);
            end
        end

        % Value changed function: filterEMGandexporttostoCheckBox
        function filterEMGandexporttostoCheckBoxValueChanged(app, event)
            value = app.filterEMGandexporttostoCheckBox.Value;
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
                app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;
            else
                app.forceplatecontactsLabel.Visible = "off";
                app.validfootstrikeeventsLabel.Visible = "on";
                app.detectforceplatesautomaticallyCheckBox.Value = 1;
                app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 1;
                answer = MFquestdlg([0.4, 0.4], 'Reprocess the c3d file to identify events and foot contacts?', 'Reprocess?', 'Yes', 'No', 'Yes');
                if strcmp(answer, 'Yes')
                    app.processC3Dfile();
                end
            end
        end
    end

    % Component initialization
    methods (Access = public)

        % Create UIFigure and components
        function createComponents(app)

            % Create PrepareC3DfileUIFigure and hide until all components are created
            app.PrepareC3DfileUIFigure = uifigure('Visible', 'off');
            app.PrepareC3DfileUIFigure.Position = [10 10 1294 829];
            app.PrepareC3DfileUIFigure.Name = 'PrepareC3Dfile';
            app.PrepareC3DfileUIFigure.Interruptible = 'off';

            % Create UIAxesSteps
            app.UIAxesSteps = uiaxes(app.PrepareC3DfileUIFigure);
            title(app.UIAxesSteps, 'Steps')
            ylabel(app.UIAxesSteps, '\color{red}left \color{black} and \color{green} right \color{black} cycles')
            app.UIAxesSteps.YTick = [-1 -0.5 0 0.5 1];
            app.UIAxesSteps.YTickLabel = {'-1'; '-0.5'; '0'; '0.5'; '1'};
            app.UIAxesSteps.Position = [15 329 1269 158];

            % Create UIAxesGRFs
            app.UIAxesGRFs = uiaxes(app.PrepareC3DfileUIFigure);
            title(app.UIAxesGRFs, 'GRFs')
            xlabel(app.UIAxesGRFs, 'frame')
            ylabel(app.UIAxesGRFs, '[N]')
            app.UIAxesGRFs.YGrid = 'on';
            app.UIAxesGRFs.Position = [15 26 1269 304];

            % Create Selectc3dfileButton
            app.Selectc3dfileButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.Selectc3dfileButton.ButtonPushedFcn = createCallbackFcn(app, @Selectc3dfileButtonPushed, true);
            app.Selectc3dfileButton.Position = [45 716 134 22];
            app.Selectc3dfileButton.Text = 'Select c3d file';

            % Create PrepareC3DfileforOpenSimsimulationsLabel
            app.PrepareC3DfileforOpenSimsimulationsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.PrepareC3DfileforOpenSimsimulationsLabel.HorizontalAlignment = 'center';
            app.PrepareC3DfileforOpenSimsimulationsLabel.FontSize = 20;
            app.PrepareC3DfileforOpenSimsimulationsLabel.FontWeight = 'bold';
            app.PrepareC3DfileforOpenSimsimulationsLabel.Position = [1 760 1283 70];
            app.PrepareC3DfileforOpenSimsimulationsLabel.Text = {'Prepare C3D file for OpenSim simulations'; 'creating marker_experimental.trc / grf.mot / grf.xml / settings.mat and a copy of original c3d-file'};

            % Create SelectaC3DfileLabel
            app.SelectaC3DfileLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.SelectaC3DfileLabel.Position = [203 716 682 22];
            app.SelectaC3DfileLabel.Text = 'Select a C3D file';

            % Create detectforceplatesautomaticallyCheckBox
            app.detectforceplatesautomaticallyCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.detectforceplatesautomaticallyCheckBox.ValueChangedFcn = createCallbackFcn(app, @detectforceplatesautomaticallyCheckBoxValueChanged, true);
            app.detectforceplatesautomaticallyCheckBox.Text = 'detect forceplates automatically';
            app.detectforceplatesautomaticallyCheckBox.Position = [45 619 191 22];
            app.detectforceplatesautomaticallyCheckBox.Value = true;

            % Create filterGRFsCheckBox
            app.filterGRFsCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.filterGRFsCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.filterGRFsCheckBox.Text = 'filter GRFs';
            app.filterGRFsCheckBox.Position = [45 577 80 22];
            app.filterGRFsCheckBox.Value = true;

            % Create setGRFstozerooutsideregionoffootcontactCheckBox
            app.setGRFstozerooutsideregionoffootcontactCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.setGRFstozerooutsideregionoffootcontactCheckBox.ValueChangedFcn = createCallbackFcn(app, @setGRFstozerooutsideregionoffootcontactCheckBoxValueChanged, true);
            app.setGRFstozerooutsideregionoffootcontactCheckBox.Text = 'set GRFs to zero outside region of foot contact and do post zero levelling';
            app.setGRFstozerooutsideregionoffootcontactCheckBox.Position = [45 598 417 22];
            app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = true;

            % Create detectwalkingdirectionautomaticallybymarkerCheckBox
            app.detectwalkingdirectionautomaticallybymarkerCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.ValueChangedFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Text = 'detect walking direction automatically by marker';
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Position = [45 640 280 22];
            app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = true;

            % Create rotatearoundyaxisDropDownLabel
            app.rotatearoundyaxisDropDownLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.rotatearoundyaxisDropDownLabel.HorizontalAlignment = 'right';
            app.rotatearoundyaxisDropDownLabel.Enable = 'off';
            app.rotatearoundyaxisDropDownLabel.Position = [396 640 109 22];
            app.rotatearoundyaxisDropDownLabel.Text = 'rotate around y axis';

            % Create rotatearoundyaxisDropDown
            app.rotatearoundyaxisDropDown = uidropdown(app.PrepareC3DfileUIFigure);
            app.rotatearoundyaxisDropDown.Items = {'0°', '90°', '180°', '270°'};
            app.rotatearoundyaxisDropDown.Enable = 'off';
            app.rotatearoundyaxisDropDown.Position = [520 640 63 22];
            app.rotatearoundyaxisDropDown.Value = '0°';

            % Create CreatefilesButton
            app.CreatefilesButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.CreatefilesButton.ButtonPushedFcn = createCallbackFcn(app, @CreatefilesButtonPushed, true);
            app.CreatefilesButton.Enable = 'off';
            app.CreatefilesButton.Position = [634 496 139 34];
            app.CreatefilesButton.Text = 'Create files';

            % Create leftListBoxLabel
            app.leftListBoxLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.leftListBoxLabel.HorizontalAlignment = 'right';
            app.leftListBoxLabel.Position = [558 601 25 22];
            app.leftListBoxLabel.Text = 'left ';

            % Create leftListBox
            app.leftListBox = uilistbox(app.PrepareC3DfileUIFigure);
            app.leftListBox.Items = {};
            app.leftListBox.Multiselect = 'on';
            app.leftListBox.ValueChangedFcn = createCallbackFcn(app, @updateCycles, true);
            app.leftListBox.Position = [598 551 100 74];
            app.leftListBox.Value = {};

            % Create validfootstrikeeventsLabel
            app.validfootstrikeeventsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.validfootstrikeeventsLabel.FontWeight = 'bold';
            app.validfootstrikeeventsLabel.Position = [647 640 136 22];
            app.validfootstrikeeventsLabel.Text = 'valid foot strike events';

            % Create rightListBoxLabel
            app.rightListBoxLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.rightListBoxLabel.HorizontalAlignment = 'right';
            app.rightListBoxLabel.Position = [728 601 29 22];
            app.rightListBoxLabel.Text = 'right';

            % Create rightListBox
            app.rightListBox = uilistbox(app.PrepareC3DfileUIFigure);
            app.rightListBox.Items = {};
            app.rightListBox.Multiselect = 'on';
            app.rightListBox.ValueChangedFcn = createCallbackFcn(app, @updateCycles, true);
            app.rightListBox.Position = [772 551 100 74];
            app.rightListBox.Value = {};

            % Create openselectedfilewithdefaultprogramButton
            app.openselectedfilewithdefaultprogramButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.openselectedfilewithdefaultprogramButton.ButtonPushedFcn = createCallbackFcn(app, @openselectedfilewithdefaultprogramButtonPushed, true);
            app.openselectedfilewithdefaultprogramButton.Visible = 'off';
            app.openselectedfilewithdefaultprogramButton.Position = [45 678 222 22];
            app.openselectedfilewithdefaultprogramButton.Text = 'open selected file with default program';

            % Create LogTextAreaLabel
            app.LogTextAreaLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.LogTextAreaLabel.HorizontalAlignment = 'right';
            app.LogTextAreaLabel.Position = [895 714 26 22];
            app.LogTextAreaLabel.Text = 'Log';

            % Create LogTextArea
            app.LogTextArea = uitextarea(app.PrepareC3DfileUIFigure);
            app.LogTextArea.Editable = 'off';
            app.LogTextArea.Position = [936 486 338 252];

            % Create filterorderSpinnerLabel
            app.filterorderSpinnerLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.filterorderSpinnerLabel.HorizontalAlignment = 'right';
            app.filterorderSpinnerLabel.Position = [156 577 59 22];
            app.filterorderSpinnerLabel.Text = 'filter order';

            % Create filterorderSpinner
            app.filterorderSpinner = uispinner(app.PrepareC3DfileUIFigure);
            app.filterorderSpinner.Step = 2;
            app.filterorderSpinner.Limits = [2 10];
            app.filterorderSpinner.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.filterorderSpinner.Position = [226 577 53 22];
            app.filterorderSpinner.Value = 4;

            % Create cutofffrequencySpinnerLabel
            app.cutofffrequencySpinnerLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinnerLabel.HorizontalAlignment = 'right';
            app.cutofffrequencySpinnerLabel.Position = [293 577 91 22];
            app.cutofffrequencySpinnerLabel.Text = 'cutoff frequency';

            % Create cutofffrequencySpinner
            app.cutofffrequencySpinner = uispinner(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner.Limits = [4 100];
            app.cutofffrequencySpinner.ValueChangedFcn = createCallbackFcn(app, @updateGRF_Plot, true);
            app.cutofffrequencySpinner.Position = [396 577 53 22];
            app.cutofffrequencySpinner.Value = 20;

            % Create detectWalkingByMarker
            app.detectWalkingByMarker = uieditfield(app.PrepareC3DfileUIFigure, 'text');
            app.detectWalkingByMarker.ValueChangedFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectWalkingByMarker.ValueChangingFcn = createCallbackFcn(app, @detectwalkingdirectionautomaticallybymarkerCheckBoxValueChanged, true);
            app.detectWalkingByMarker.Position = [324 640 73 22];
            app.detectWalkingByMarker.Value = 'LTOE';

            % Create detectFPbyMarker
            app.detectFPbyMarker = uieditfield(app.PrepareC3DfileUIFigure, 'text');
            app.detectFPbyMarker.ValueChangedFcn = createCallbackFcn(app, @detectforceplatesautomaticallyCheckBoxValueChanged, true);
            app.detectFPbyMarker.Position = [245 619 73 22];
            app.detectFPbyMarker.Value = 'LTOE';

            % Create v16Label
            app.v16Label = uilabel(app.PrepareC3DfileUIFigure);
            app.v16Label.Position = [1246 784 28 22];
            app.v16Label.Text = 'v1.6';

            % Create filterEMGandexporttostoCheckBox
            app.filterEMGandexporttostoCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.filterEMGandexporttostoCheckBox.ValueChangedFcn = createCallbackFcn(app, @filterEMGandexporttostoCheckBoxValueChanged, true);
            app.filterEMGandexporttostoCheckBox.Text = 'filter EMG and export to *.sto';
            app.filterEMGandexporttostoCheckBox.Position = [45 521 176 22];
            app.filterEMGandexporttostoCheckBox.Value = true;

            % Create SelectexistingCSVButton
            app.SelectexistingCSVButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.SelectexistingCSVButton.ButtonPushedFcn = createCallbackFcn(app, @SelectexistingCSVButtonPushed, true);
            app.SelectexistingCSVButton.Position = [371 496 121 22];
            app.SelectexistingCSVButton.Text = 'Select existing CSV';

            % Create createCSVandopenButton
            app.createCSVandopenButton = uibutton(app.PrepareC3DfileUIFigure, 'push');
            app.createCSVandopenButton.ButtonPushedFcn = createCallbackFcn(app, @createCSVandopenButtonPushed, true);
            app.createCSVandopenButton.Position = [217 496 130 22];
            app.createCSVandopenButton.Text = 'create CSV and open';

            % Create renameEMGlabelsCheckBox
            app.renameEMGlabelsCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.renameEMGlabelsCheckBox.Text = 'rename EMG labels';
            app.renameEMGlabelsCheckBox.Position = [71 496 128 22];
            app.renameEMGlabelsCheckBox.Value = true;

            % Create EMGfactorSpinnerLabel
            app.EMGfactorSpinnerLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.EMGfactorSpinnerLabel.HorizontalAlignment = 'right';
            app.EMGfactorSpinnerLabel.Position = [253 521 66 22];
            app.EMGfactorSpinnerLabel.Text = 'EMG factor';

            % Create EMGfactorSpinner
            app.EMGfactorSpinner = uispinner(app.PrepareC3DfileUIFigure);
            app.EMGfactorSpinner.Limits = [1 10000];
            app.EMGfactorSpinner.Tooltip = {'Use this factor if output signal has steps.'; 'Of course, you need to have the same factor for all trials of the same participant!'};
            app.EMGfactorSpinner.Position = [330 521 119 22];
            app.EMGfactorSpinner.Value = 1;

            % Create forceplatecontactsLabel
            app.forceplatecontactsLabel = uilabel(app.PrepareC3DfileUIFigure);
            app.forceplatecontactsLabel.FontWeight = 'bold';
            app.forceplatecontactsLabel.Visible = 'off';
            app.forceplatecontactsLabel.Position = [647 640 136 22];
            app.forceplatecontactsLabel.Text = 'force plate contacts';

            % Create Ignorec3deventsCheckBox
            app.Ignorec3deventsCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.Ignorec3deventsCheckBox.ValueChangedFcn = createCallbackFcn(app, @Ignorec3deventsCheckBoxValueChanged, true);
            app.Ignorec3deventsCheckBox.Text = 'Ignore c3d events and set forceplate contacts manually';
            app.Ignorec3deventsCheckBox.Position = [582 678 320 22];

            % Create ClimbingCheckBox
            app.ClimbingCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.ClimbingCheckBox.Text = 'Climbing';
            app.ClimbingCheckBox.Position = [277 678 69 22];

            % Create filtermarkersCheckBox
            app.filtermarkersCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.filtermarkersCheckBox.Text = 'filter markers';
            app.filtermarkersCheckBox.Position = [45 551 92 22];

            % Create filterorderSpinnerLabel_Marker
            app.filterorderSpinnerLabel_Marker = uilabel(app.PrepareC3DfileUIFigure);
            app.filterorderSpinnerLabel_Marker.HorizontalAlignment = 'right';
            app.filterorderSpinnerLabel_Marker.Position = [156 551 59 22];
            app.filterorderSpinnerLabel_Marker.Text = 'filter order';

            % Create filterorderSpinner_Marker
            app.filterorderSpinner_Marker = uispinner(app.PrepareC3DfileUIFigure);
            app.filterorderSpinner_Marker.Step = 2;
            app.filterorderSpinner_Marker.Limits = [2 10];
            app.filterorderSpinner_Marker.Position = [226 551 53 22];
            app.filterorderSpinner_Marker.Value = 2;

            % Create cutofffrequencySpinner_2Label
            app.cutofffrequencySpinner_2Label = uilabel(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner_2Label.HorizontalAlignment = 'right';
            app.cutofffrequencySpinner_2Label.Position = [293 551 91 22];
            app.cutofffrequencySpinner_2Label.Text = 'cutoff frequency';

            % Create cutofffrequencySpinner_Marker
            app.cutofffrequencySpinner_Marker = uispinner(app.PrepareC3DfileUIFigure);
            app.cutofffrequencySpinner_Marker.Limits = [4 30];
            app.cutofffrequencySpinner_Marker.Position = [396 551 53 22];
            app.cutofffrequencySpinner_Marker.Value = 5;

            % Create normalizeEMGCheckBox
            app.normalizeEMGCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.normalizeEMGCheckBox.Text = 'normalize EMG';
            app.normalizeEMGCheckBox.Position = [461 521 105 22];
            app.normalizeEMGCheckBox.Value = true;

            % Create invertforcesensorsCheckBox
            app.invertforcesensorsCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.invertforcesensorsCheckBox.Text = 'invert force sensors';
            app.invertforcesensorsCheckBox.Position = [360 678 126 22];

            % Create mirrorCheckBox
            app.mirrorCheckBox = uicheckbox(app.PrepareC3DfileUIFigure);
            app.mirrorCheckBox.Tooltip = {'EMG stays the same!'};
            app.mirrorCheckBox.Text = 'mirror';
            app.mirrorCheckBox.Position = [491 678 53 22];

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