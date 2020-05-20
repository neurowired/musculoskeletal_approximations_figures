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
    ax(1) = subplot(2, 2, 1); hold on;
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
    ax(2) = subplot(2, 2, 3); hold on;
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
        
    min_l = min([min(dat_l), min(min(app_l))]);
    max_l = max([max(dat_l), max(max(app_l))]);
    l_r0 = min_l+0.9*(max_l-min_l);
    fmax = 122.5;
    fmaxpass = 12.25;
    u = 1;
    fl = @(l_norm) 2.5 * l_norm - 1.25 * l_norm .^2;
    fpass = @(l) (exp(2*(l-l_r0)/(max_l-min_l)) - 1) ./ (exp(1)-1) .* l>l_r0;
    norm = @(l) (l-min_l) ./ (max_l - min_l);
    f = (fmax .* fl(norm(app_l))*u +fmaxpass*fpass(app_l));
    tor_z1 = app_z1 .* f / 1000;
    tor_z2 = app_z2 .* f / 1000;
    
    ax(3) = subplot(2, 2, 2); hold on;
        surf(app_x, app_y, f, 'edgecolor', 'r', 'facecolor', 'none');
        set(gca, 'xtick', metaDataConfig.xtick)
        set(gca, 'xticklabel', metaDataConfig.xticklabel)
        xlabel([sDof1 ', rad']);
        set(gca, 'ytick', metaDataConfig.xtick)
        set(gca, 'yticklabel', metaDataConfig.xticklabel)
        ylabel([sDof2 ', rad']);
        zlabel('Force, N');
        axis tight;
        xlim([-pi/2, pi/2]);
        ylim([-pi/2, pi/2]);
        zlim([0, 160]);
        grid on;
        view(64,31);
    
    ax(4) = subplot(2, 2, 4); hold on;
        surf(app_x, app_y, tor_z1, 'edgecolor', 'g', 'facecolor', 'none');
        surf(app_x, app_y, tor_z2, 'edgecolor', 'b', 'facecolor', 'none');
        set(gca, 'xtick', metaDataConfig.xtick)
        set(gca, 'xticklabel', metaDataConfig.xticklabel)
        xlabel([sDof1 ', rad']);
        set(gca, 'ytick', metaDataConfig.xtick)
        set(gca, 'yticklabel', metaDataConfig.xticklabel)
        ylabel([sDof2 ', rad']);
        zlabel('Torque, Nm');
        axis tight;
        xlim([-pi/2, pi/2]);
        ylim([-pi/2, pi/2]);
        legend(sDof1, sDof2);
        grid on;
        view(64,31);
end
