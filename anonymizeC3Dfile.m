clear; 

inputFilePath = 'dynamic03.c3d'; 
outputFilePath = 'output_no_subject_info.c3d';

h = btkReadAcquisition(inputFilePath);
metaData = btkGetMetaData(h);

if isfield(metaData.children, 'SUBJECTS')
    if isfield(metaData.children.SUBJECTS.children, 'NAMES')
        btkRemoveMetaData(h, 'SUBJECTS', 'NAMES');
    end
end

numEvents = btkGetEventNumber(h);
for i = 1 : numEvents
    btkSetEventSubject(h, i, 'ANONYM');
end

btkWriteAcquisition(h, outputFilePath);

btkCloseAcquisition(h);

disp('C3D file saved');
