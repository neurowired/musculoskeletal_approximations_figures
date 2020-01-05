function [ metaMuscleConfig ] = extractMetaMuscle()
%EXTRACTMETAMUSCLE Imports metaMuscle in a comfortable way
    opts = detectImportOptions('arm_metaMuscle.csv');
    opts.LineEnding = {'\r\n'};
    mmt = readtable('arm_metaMuscle.csv', opts);
    
    metaMuscleConfig.nMuscles = height(mmt);
    metaMuscleConfig.U = metaMuscleConfig.nMuscles;
    metaMuscleConfig.sMuscle = mmt.sMuscle;
    metaMuscleConfig.idDOFList = {};
    for iMuscle=1:metaMuscleConfig.U
        metaMuscleConfig.idDOFList{iMuscle} = str2num(mmt.idDOFList{iMuscle});
    end
    metaMuscleConfig.muscle2id = containers.Map(mmt.sMuscle, mmt.idMuscle);
end

