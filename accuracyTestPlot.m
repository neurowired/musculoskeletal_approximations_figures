function accuracyTestPlot( filename, bPlotFull, bPlotPerMuscle, bPlotFullNO )
    load(filename);

    % colors
    spliFaceColor = [0 0 0];  % black  % [30 136 229]./255;  % blue
    spliFaceAlpha = 1;
    spliEdgeAlpha = 0;

    funcFaceColor = [117 117 117]./255;  % gray  % [229 57 53]./255;  % red
    funcFaceAlpha = .4;
    funcEdgeAlpha = 0;

    intcFaceColor = [67 160 71]./255;  % green
    intcFaceAlpha = 0;
    intcEdgeAlpha = 1;

    numBins = 35;
    numBinsSpline = 36;

    % data processing
    alph = 0.01;

    muscleIndexMam = [muscleIndexTotal; fitmuscleIndexTotal];
    muscleIndexLen = [muscleIndexLenTotal; fitmuscleIndexLenTotal];
    dofIndex = [dofIndexTotal; fitdofIndexTotal];

    fitlenDiffSNTotal = zeros(size(fitlenDiffFNTotal));
    totLenS = [lenDiffSNTotal; fitlenDiffSNTotal];
    [totLenSNO, numLenSOut, meeLenS, meeLenSNO, totLenSOuMids, ~] = accuracyTestStats(totLenS, alph, muscleIndexLen, []);
    totLenF = [lenDiffFNTotal; fitlenDiffFNTotal];
    [totLenFNO, numLenFOut, meeLenF, meeLenFNO, totLenFOuMids, ~] = accuracyTestStats(totLenF, alph, muscleIndexLen, []);
    totLenC = [lenDiffCNTotal; fitlenDiffCNTotal];
    [totLenCNO, numLenCOut, meeLenC, meeLenCNO, totLenCOuMids, ~] = accuracyTestStats(totLenC, alph, muscleIndexLen, []);

    fitmomarmDiffSNTotal = zeros(size(fitmomarmDiffFNTotal));
    totMamS = [momarmDiffSNTotal; fitmomarmDiffSNTotal];
    [totMamSNO, numMamSOut, meeMamS, meeMamSNO, totMamSOuMids, totMamSOuDids] = accuracyTestStats(totMamS, alph, muscleIndexMam, dofIndex);
    totMamF = [momarmDiffFNTotal; fitmomarmDiffFNTotal];
    [totMamFNO, numMamFOut, meeMamF, meeMamFNO, totMamFOuMids, totMamFOuDids] = accuracyTestStats(totMamF, alph, muscleIndexMam, dofIndex);
    totMamC = [momarmDiffCNTotal; fitmomarmDiffCNTotal];
    [totMamCNO, numMamCOut, meeMamC, meeMamCNO, totMamCOuMids, totMamCOuDids] = accuracyTestStats(totMamC, alph, muscleIndexMam, dofIndex);

    if bPlotFull
        figure(); hold on;
        subplot(2, 1, 1); hold on;
        histogram(100*totLenS, numBinsSpline, ...
                  'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                  'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totLenF, numBins, ...
                  'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                  'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totLenC, numBins, ...
                  'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                  'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
        ylim([0, 1]);
        legend('Splines', 'Functions', 'Consistent Functions');
        title('length')
        xlabel('Error, %')
        ylabel('Portion of data points')

        subplot(2, 1, 2); hold on;
        histogram(100*totMamS, numBinsSpline, ...
                  'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                  'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totMamF, numBins, ...
                  'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                  'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totMamC, numBins, ...
                  'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                  'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
        ylim([0, 1]);
        legend('Splines', 'Functions', 'Consistent Functions');
        title('moment arm')
        xlabel('Error, %')
        ylabel('Portion of data points')
    end

    if bPlotFullNO
        xlims = 0.6;
        xlimsSpline = 0.3;
        bins = linspace(-xlims, xlims, numBins);
        binsSpline = linspace(-xlimsSpline, xlimsSpline, numBinsSpline);
        figure(); hold on
        axs(1) = subplot(2, 1, 1); hold on;
        histogram(100*totLenSNO, binsSpline, ...
                  'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                  'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totLenFNO, bins, ...
                  'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                  'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totLenCNO, bins, ...
                  'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                  'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
        xlim([-xlims; xlims]);
        ylim([0, 1]);
        legend('Splines', 'Functions', 'Consistent Functions');
        title('length')
        xlabel('Error, %')
        ylabel('Portion of data points')

        xlims = 8;
        xlimsSpline = 4;
        bins = linspace(-xlims, xlims, numBins);
        binsSpline = linspace(-xlimsSpline, xlimsSpline, numBinsSpline);
        axs(2) = subplot(2, 1, 2); hold on;
        histogram(100*totMamSNO, binsSpline, ...
                  'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                  'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totMamFNO, bins, ...
                  'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                  'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
        histogram(100*totMamCNO, bins, ...
                  'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                  'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
        xlim([-xlims; xlims]);
        ylim([0, 1]);
        legend('Splines', 'Functions', 'Consistent Functions');
        legend('Splines', 'Functions', 'Consistent Functions');
        title('moment arm')
        xlabel('Error, %')
    end

    if bPlotPerMuscle
        uMuscleIds = unique(muscleIndexLen);
        xnSubplots = ceil(sqrt(length(uMuscleIds)));
        ynSubplots = ceil(length(uMuscleIds) / xnSubplots);

        % len
        figure(); hold on;
        axs = zeros(2*length(uMuscleIds));
        suptitle('Muscle length errors per muscle on test data.');
        for iMuscle=1:length(uMuscleIds)
            axs(iMuscle) = subplot(xnSubplots, ynSubplots, iMuscle); hold on;
            histogram(lenDiffSNTotal(muscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBinsSpline, ...
                      'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                      'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
            histogram(lenDiffFNTotal(muscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                      'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
            histogram(lenDiffCNTotal(muscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                      'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
            title(num2str(uMuscleIds(iMuscle)));
        end

        figure(); hold on;
        suptitle('Muscle length errors per muscle on fit data.');
        for iMuscle=1:length(uMuscleIds)
            axs(length(uMuscleIds)+iMuscle) = subplot(xnSubplots, ynSubplots, iMuscle); hold on;
            histogram(fitlenDiffSNTotal(fitmuscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBinsSpline, ...
                      'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                      'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
            histogram(fitlenDiffFNTotal(fitmuscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                      'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
            histogram(fitlenDiffCNTotal(fitmuscleIndexLenTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                      'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
            title(num2str(uMuscleIds(iMuscle)));
        end
        linkaxes(axs, 'x');

        % mam
        figure(); hold on;
        suptitle('Moment arms errors per muscle on test data.');
        for iMuscle=1:length(uMuscleIds)
            axs(iMuscle) = subplot(xnSubplots, ynSubplots, iMuscle); hold on;
            histogram(momarmDiffSNTotal(muscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBinsSpline, ...
                      'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                      'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
            histogram(momarmDiffFNTotal(muscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                      'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
            histogram(momarmDiffCNTotal(muscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                      'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
            title(num2str(uMuscleIds(iMuscle)));
        end

        figure(); hold on;
        suptitle('Moment arms errors per muscle on fit data.');
        for iMuscle=1:length(uMuscleIds)
            axs(iMuscle+length(uMuscleIds)) = subplot(xnSubplots, ynSubplots, iMuscle); hold on;
            histogram(fitmomarmDiffSNTotal(fitmuscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBinsSpline, ...
                      'FaceColor', spliFaceColor, 'FaceAlpha', spliFaceAlpha, ...
                      'EdgeAlpha', spliEdgeAlpha, 'Normalization', 'probability');
            histogram(fitmomarmDiffFNTotal(fitmuscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', funcFaceColor, 'FaceAlpha', funcFaceAlpha, ...
                      'EdgeAlpha', funcEdgeAlpha, 'Normalization', 'probability');
            histogram(fitmomarmDiffCNTotal(fitmuscleIndexTotal==uMuscleIds(iMuscle)), ...
                      numBins, ...
                      'FaceColor', intcFaceColor, 'FaceAlpha', intcFaceAlpha, ...
                      'EdgeAlpha', intcEdgeAlpha, 'Normalization', 'probability');
            title(num2str(uMuscleIds(iMuscle)));
        end
        linkaxes(axs, 'x');
    end
end
