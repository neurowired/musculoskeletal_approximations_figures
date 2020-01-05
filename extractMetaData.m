function [ metaDataConfig ] = extractMetaData(  )
%EXTRACTMETADATA Imports metaDof in a comfortable way
    % Data location
    metaDataConfig.dataDirectory = 'DATA 9 Point/';
    metaDataConfig.testdataDirectory = 'TESTDATA 8 Point/';
    metaDataConfig.totalPlotDirectory = 'totalPlot/';

    % Functions
    metaDataConfig.functionsLocation = 'out/';

    % Presets for plots and figures
    metaDataConfig.xtick        = [-pi;     -pi/2;      0;      pi/2;       pi];
    metaDataConfig.xticklabel   = {'-\pi';  '-\pi/2';   '0';    '\pi/2';    '\pi'};
    metaDataConfig.xtick_wp4        = [-pi;     -pi/2;      -pi/4;      0;      pi/4;       pi/2;       pi];
    metaDataConfig.xticklabel_wp4   = {'-\pi';  '-\pi/2';   '-\pi/4';   '0';    '\pi/4';    '\pi/2';    '\pi'};

    % metaDOF
    md = readtable('arm_metaDOF.csv');

    metaDataConfig.M = height(md);
    metaDataConfig.modelLoaded = metaDataConfig.M;
    metaDataConfig.idDOF = md.idDOF;
    metaDataConfig.sDof = md.sDOF;
    metaDataConfig.rangeMin = md.nRangeMin;
    metaDataConfig.rangeMax = md.nRangeMax;
    metaDataConfig.rangeWidth = md.nRangeMax - md.nRangeMin;
    metaDataConfig.ic = md.nIC;

    metaDataConfig.signal2id = containers.Map(metaDataConfig.sDof, metaDataConfig.idDOF);

end

