folder = 'C:\Users\Willi\ucloud\PhD\Toewalkers_Andreas_GaitPhases\Projekt Andreas\ID_2\Static FB Anterior - CAST 1';

markerData = load_marker_trc(fullfile(folder, 'marker_experimental.trc'));
markerNames = fieldnames(markerData);

% markersOfInterest = {'R_IAS', 'L_IAS', 'SACR', 'L_FLE', 'R_FLE', 'LKNM', 'RKNM', 'L_FAL', 'R_FAL', 'LANM', 'RANM'};
markersOfInterest = {'R_IAS', 'L_IAS', 'L_IPS', 'R_IPS', 'L_FLE', 'R_FLE', 'L_FME', 'R_FME', 'L_FAL', 'R_FAL', 'L_TAM', 'R_TAM'};
markerDataOfInterest = [];
for i = 1 :numel(markersOfInterest)
    markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
    markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
    markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
end

LKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_FLE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_FME'), :, :)))/2;
RKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_FLE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_FME'), :, :)))/2;
LAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_FAL'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_TAM'), :, :)))/2;
RAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_FAL'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_TAM'), :, :)))/2;
JointCenters.LKJC=LKJC;
JointCenters.RKJC=RKJC;

% add AJC_offset marker based on Bruening et al., 2008
% right AJC_offset
SL=RKJC-RAJC;
SL_distance=sqrt( SL(1,1)*SL(1,1)+SL(1,2)*SL(1,2)+SL(1,3)*SL(1,3) );
offset = 0.027*SL_distance;
for t=1:size(RKJC, 1)
    originRTibia=RAJC;
    [e1TibiaR,e2TibiaR,e3TibiaR]=segmentorientation_r(RKJC(t,:)-RAJC(t,:),squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_FAL'), t, :))'-RAJC(t,:));
    rotationmatrix_1=[e1TibiaR;e2TibiaR;e3TibiaR];
    RAJC_offset_rTibia(t,:)=[-offset 0 0];
    RAJC_offset(t,:)=(mldivide(rotationmatrix_1,(RAJC_offset_rTibia(t,:)'))+originRTibia(t,:)')';
    JointCenters.RAJC(t,:)=RAJC_offset(t,:);
end

% left AJC_offset
SL_left=LKJC-LAJC;
SL_left_distance=sqrt( SL_left(1,1)*SL_left(1,1)+SL_left(1,2)*SL_left(1,2)+SL_left(1,3)*SL_left(1,3) );
offset_left = 0.027*SL_left_distance;
for t=1:size(LKJC, 1)
    originLTibia=LAJC;
    [e1TibiaL,e2TibiaL,e3TibiaL]=segmentorientation_l(LKJC(t,:)-LAJC(t,:),squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_FAL'), t, :))'-LAJC(t,:));
    rotationmatrix_2=[e1TibiaL;e2TibiaL;e3TibiaL];
    LAJC_offset_lTibia(t,:)=[-offset_left 0 0];
    LAJC_offset(t,:)=(mldivide(rotationmatrix_2,(LAJC_offset_lTibia(t,:)'))+originLTibia(t,:)')';
    JointCenters.LAJC(t,:)=LAJC_offset(t,:);
end

% add Harrington HJC
ASISvector=squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_IAS'), :, :))-squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_IAS'), :, :));
ASISdistance=sqrt((ASISvector(1,1))^2+(ASISvector(1,2))^2+(ASISvector(1,3))^2);
MidASIS = (squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_IAS'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_IAS'), :, :)))/2;
if any(strcmp(markersOfInterest, 'SACR'))
    MidPSIS = squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'SACR'), :, :));
else
    MidPSIS = (squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_IPS'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_IPS'), :, :)))/2;
end

PW=ASISdistance;
%   Define the pelvis origin
originPelvis=MidASIS;
for t=1:size(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_IAS'), :, :)), 1)
    % HJC definition based on modified Harrington equations (global
    % coordinates) Sangeux, 2015
    H_ap=-0.138*PW-10.4;
    H_v=-0.305*PW-10.9;
    H_ml=0.33*PW+7.3;
    H_ml_left=-H_ml;
    RHJC_Harringtion(t,:) = horzcat(H_ap,H_v,H_ml);
    LHJC_Harringtion(t,:) = horzcat(H_ap,H_v,H_ml_left);

    %   Calculate unit vectors of Pelvis segment relative to global
    [e1Pelvis,e2Pelvis,e3Pelvis]=segmentorientation_1Frame(originPelvis(t,:)-MidPSIS(t,:),(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'R_IAS'), t, :))-squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'L_IAS'), t, :)))');
    rotationmatrix=[e1Pelvis;e2Pelvis;e3Pelvis];
    RHJC(t,:)=(mldivide(rotationmatrix,(RHJC_Harringtion(t,:)'))+originPelvis(t,:)')';
    LHJC(t,:)=(mldivide(rotationmatrix,(LHJC_Harringtion(t,:)'))+originPelvis(t,:)')';
    JointCenters.RHJC(t,:)=RHJC(t,:);
    JointCenters.LHJC(t,:)=LHJC(t,:);
    rotationmatrix=[];
end

text = fileread(fullfile(folder, 'marker_experimental.trc'));
lines = strsplit(text, '\n');
data = [];
for i = 1 : numel(lines)
    data{i} = strsplit(lines{i}, '\t', 'CollapseDelimiters', false);
end

jointCenterNames = fieldnames(JointCenters);
origNrOfCols = numel(data{7});
for row = 7 : numel(data) - 1
    for marker = 1 : numel(jointCenterNames)
        col = origNrOfCols + (marker-1)*3;
        if marker < numel(jointCenterNames)
            data{4}{col} = [jointCenterNames{marker} 'WK'];
        else
            data{4}{col} = [jointCenterNames{marker} 'WK '];
        end
        for xyz = 1 : 3
            col = origNrOfCols + (marker-1)*3 + xyz - 1;
            data{row}{col} = JointCenters.(jointCenterNames{marker})(row-6, xyz);
            switch xyz
                case 1
                    data{5}{col} = ['X' num2str((origNrOfCols/ 3) + marker - 1) ' '];
                case 2
                    data{5}{col} = ['Y' num2str((origNrOfCols/ 3) + marker - 1) ' '];
                case 3
                    data{5}{col} = ['Z' num2str((origNrOfCols/ 3) + marker - 1) ' '];
            end
        end
    end
end

data{1}{4} = strrep(data{1}{4}, '\', '\\');
data{3}{4} = str2double(data{3}{4}) + 6;

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

fid = fopen(fullfile(folder, 'marker_experimental_with_JointCenters.trc'),'wt');
fprintf(fid, finalStr);
fclose(fid);