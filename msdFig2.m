function msdFig2(  )
%MSDFIG1 Summary of this function goes here
%   Detailed explanation goes here
    set(0, 'DefaultTextInterpreter', 'none');
    set(0, 'DefaultLegendInterpreter', 'none');

    metaDataConfig = extractMetaData();
    metaMuscleConfig = extractMetaMuscle();

    % load data
    load('msdFig2.mat');

    % plot
    figure(); hold on;
    suptitle(char(metaMuscleConfig.sMuscle{muscleId}))
    ax(1) = subplot(2, 1, 1); hold on;
        scatter3(dat_x, dat_y, dat_l, 'r');
        surf(app_x, app_y, app_l, 'edgecolor', 'r', 'facecolor', 'none');
        set(gca, 'xtick', metaDataConfig.xtick)
        set(gca, 'xticklabel', metaDataConfig.xticklabel)
        xlabel([sDof1 ', rad']);
        set(gca, 'ytick', metaDataConfig.xtick)
        set(gca, 'yticklabel', metaDataConfig.xticklabel)
        ylabel([sDof2 ', rad']);
        zlabel('Length, cm');
        axis tight;
        xlim([-pi/2, pi/2]);
        ylim([-pi/2, pi/2]);
        zlim([-inf, 34]);
        grid on;
        view(64,31);
    ax(2) = subplot(2, 1, 2); hold on;
        scatter3(dat_x, dat_y, dat_z1, 'g');
        scatter3(dat_x, dat_y, dat_z2, 'b');
        surf(app_x, app_y, app_z1, 'edgecolor', 'g', 'facecolor', 'none');
        surf(app_x, app_y, app_z2, 'edgecolor', 'b', 'facecolor', 'none');
        set(gca, 'xtick', metaDataConfig.xtick)
        set(gca, 'xticklabel', metaDataConfig.xticklabel)
        xlabel([sDof1 ', rad']);
        set(gca, 'ytick', metaDataConfig.xtick)
        set(gca, 'yticklabel', metaDataConfig.xticklabel)
        ylabel([sDof2 ', rad']);
        zlabel('Moment arms, mm');
        axis tight;
        xlim([-pi/2, pi/2]);
        ylim([-pi/2, pi/2]);
        legend(sDof1, sDof2);
        grid on;
        view(64,31);
end
