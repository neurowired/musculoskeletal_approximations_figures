function [totNO, numOut, mee, meeNO, totOuMids, totOuDids] = accuracyTestStats(tot, alph, muscleIndexTotal, dofIndexTotal)
    totAbs = abs(tot);
    totAbsMea = mean(totAbs);
    totAbsStd = std(totAbs);
    k = 1./sqrt(alph);
    mee = totAbsMea + k * totAbsStd;  % maximum expected error

    totNO = tot;
    totOu = abs(totNO) > mee;
    totNO(totOu) = [];
    numOut = length(tot) - length(totNO);
    meeNO = mean(abs(totNO)) + k * std(abs(totNO));

    totOuMids = unique(muscleIndexTotal(totOu));
    if ~isempty(dofIndexTotal)
        totOuDids = unique(dofIndexTotal(totOu));
    else
        totOuDids = [];
    end
end
