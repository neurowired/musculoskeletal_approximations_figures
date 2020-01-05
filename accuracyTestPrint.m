function accuracyTestPrint()
    filename = 'datatest33.mat';
    bPrintDetails = 0;
    load(filename);

    fprintf('NTotal: %d; NTotal fit: %d.\n', NTotal, length(fitlenDiffFNTotal));
    fprintf('Times using tic-toc (unreliable):\n');
    fprintf('\tTotal time of spline eval:                %f s\n', tsTotal);
    fprintf('\tTotal time of function eval:              %f s\n', tfTotal);
    fprintf('\tTotal time of consistent function eval:   %f s\n', tcTotal);
    fprintf('\tAverage time of single spline eval:                %f ms\n', tsTotal/NTotal*1000);
    fprintf('\tAverage time of single function eval:              %f ms\n', tfTotal/NTotal*1000);
    fprintf('\tAverage time of single consistent function eval:   %f ms\n', tcTotal/NTotal*1000);
    fprintf('Times using profiler:\n');
    fprintf('\tSplines        mean: %f ms, std: %f ms.\n', mean(tProfsTotal)*1000, std(tProfsTotal)*1000);
    fprintf('\tFunction       mean: %f ms, std: %f ms.\n', mean(tProffTotal)*1000, std(tProffTotal)*1000);
    fprintf('\tConsfunction   mean: %f ms, std: %f ms.\n', mean(tProfcTotal)*1000, std(tProfcTotal)*1000);

    alph = 0.01;

    muscleIndexMam = [muscleIndexTotal; fitmuscleIndexTotal];
    muscleIndexLen = [muscleIndexLenTotal; fitmuscleIndexLenTotal];
    dofIndex = [dofIndexTotal; fitdofIndexTotal];

    fprintf('\nLengths:\n');
    totLenS = [lenDiffSNTotal; zeros(size(fitlenDiffFNTotal))];
    [totLenSNO, numLenSOut, meeLenS, meeLenSNO, totLenSOuMids, ~] = accuracyTestStats(totLenS, alph, muscleIndexLen, []);

    totLenF = [lenDiffFNTotal; fitlenDiffFNTotal];
    [totLenFNO, numLenFOut, meeLenF, meeLenFNO, totLenFOuMids, ~] = accuracyTestStats(totLenF, alph, muscleIndexLen, []);

    totLenC = [lenDiffCNTotal; fitlenDiffCNTotal];
    [totLenCNO, numLenCOut, meeLenC, meeLenCNO, totLenCOuMids, ~] = accuracyTestStats(totLenC, alph, muscleIndexLen, []);

    [~, ttestSC] = ttest2(totLenSNO, totLenCNO, 'Alpha', alph, 'Tail', 'both', 'Vartype', 'unequal');
    [~, ttestFC] = ttest2(totLenFNO, totLenCNO, 'Alpha', alph, 'Tail', 'both', 'Vartype', 'unequal');

    if bPrintDetails
        fprintf('\tSplines DETAILED\n')
        printPerMuscleLenReport(totLenS, muscleIndexLen);
        fprintf('\tFunctions DETAILED\n')
        printPerMuscleLenReport(totLenF, muscleIndexLen);
        fprintf('\tConsistent Functions DETAILED\n')
        printPerMuscleLenReport(totLenC, muscleIndexLen);
    end

    fprintf('\tSplines\n');
    printReport(totLenS, totLenSNO, numLenSOut, meeLenS, meeLenSNO, totLenSOuMids, []);
    fprintf('\tFunctions\n');
    printReport(totLenF, totLenFNO, numLenFOut, meeLenF, meeLenFNO, totLenFOuMids, []);
    fprintf('\tConsistent Functions\n');
    printReport(totLenC, totLenCNO, numLenCOut, meeLenC, meeLenCNO, totLenCOuMids, []);
    fprintf('\tP-value t-test between splines and constrained polynomials: %f.\n', ttestSC);
    fprintf('\tP-value t-test between unconstrained and constrained polynomials: %f.\n', ttestFC);


    fprintf('\nMoment arms:\n');
    totMamS = [momarmDiffSNTotal; zeros(size(fitmomarmDiffFNTotal))];
    [totMamSNO, numMamSOut, meeMamS, meeMamSNO, totMamSOuMids, totMamSOuDids] = accuracyTestStats(totMamS, alph, muscleIndexMam, dofIndex);

    totMamF = [momarmDiffFNTotal; fitmomarmDiffFNTotal];
    [totMamFNO, numMamFOut, meeMamF, meeMamFNO, totMamFOuMids, totMamFOuDids] = accuracyTestStats(totMamF, alph, muscleIndexMam, dofIndex);

    totMamC = [momarmDiffCNTotal; fitmomarmDiffCNTotal];
    [totMamCNO, numMamCOut, meeMamC, meeMamCNO, totMamCOuMids, totMamCOuDids] = accuracyTestStats(totMamC, alph, muscleIndexMam, dofIndex);

    [~, ttestSC] = ttest2(totMamSNO, totMamCNO, 'Alpha', alph, 'Tail', 'both', 'Vartype', 'unequal');
    [~, ttestFC] = ttest2(totMamFNO, totMamCNO, 'Alpha', alph, 'Tail', 'both', 'Vartype', 'unequal');

    fprintf('\tSplines\n');
    printReport(totMamS, totMamSNO, numMamSOut, meeMamS, meeMamSNO, totMamSOuMids, totMamSOuDids);
    fprintf('\tFunctions\n');
    printReport(totMamF, totMamFNO, numMamFOut, meeMamF, meeMamFNO, totMamFOuMids, totMamFOuDids);
    fprintf('\tConsistent Functions\n');
    printReport(totMamC, totMamCNO, numMamCOut, meeMamC, meeMamCNO, totMamCOuMids, totMamCOuDids);
    fprintf('\tP-value distance between splines and constrained polynomials: %f.\n', ttestSC);
    fprintf('\tP-value distance between unconstrained and constrained polynomials: %f.\n\n', ttestFC);

end

function printReport(tot, totNO, numOut, mee, meeNO, totOuMids, totOuDids)
    fprintf('\t\t              | mean(abs), %%  |  std(abs), %%  |  max(abs), %%  |  MEE(abs), %%  |  rms(abs), %%  |  std(rms(abs)), %%  |\n');
    fprintf('\t\t  full data   |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f       |\n', ...
            100*mean(abs(tot)), 100*std(abs(tot)), 100*max(abs(tot)), ...
            mee*100, 100*sqrt(mean(tot.^2)), 100*sqrt(std(tot.^2)));
    fprintf('\t\t w/o outliers |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f  |  % 11.8f       |\n', ...
            100*mean(abs(totNO)), 100*std(abs(totNO)), 100*max(abs(totNO)), ...
            meeNO*100, 100*sqrt(mean(totNO.^2)), 100*sqrt(std(totNO.^2)));
    fprintf('\t\tPortion of outliers: %f%%, cutoff val (MEE(abs)): %f%%.\n', numOut/length(tot)*100, mee*100);
    fprintf('\t\tOutlier muscle ids:');
    for totOuMid=totOuMids
        fprintf(' %2d', totOuMid);
    end
    fprintf('\n');
    if ~isempty(totOuDids)
        fprintf('\t\tOutlier dof ids:');
        for totOuDid=totOuDids
            fprintf(' %2d', totOuDid);
        end
        fprintf('\n');
    end
end


function printPerMuscleLenReport(tot, muscleIndex)
    fprintf('\t\t  idMuscle    | mean(abs), %%  |  std(abs), %%  |  max(abs), %%  |  rms(abs), %%  |  std(rms(abs)), %%  |\n');
    uMuscles = unique(muscleIndex)';
    for iMuscle=uMuscles
        iMuscleSubset = iMuscle == muscleIndex;
        fprintf('\t\t  %9d   |  % 11.6f  |  % 11.6f  |  % 11.6f  |  % 11.6f  |  % 11.6f       |\n', ...
                iMuscle, 100*mean(abs(tot(iMuscleSubset))), 100*std(abs(tot(iMuscleSubset))), 100*max(abs(tot(iMuscleSubset))), ...
                100*sqrt(mean(tot(iMuscleSubset).^2)), 100*sqrt(std(tot(iMuscleSubset).^2)));
    end
end
