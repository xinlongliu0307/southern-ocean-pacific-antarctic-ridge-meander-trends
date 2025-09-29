% -----------------------------------------------------------------------
%  Helper: produce one figure (raw series, Sen trend, annotation) -------
% -----------------------------------------------------------------------
function plotWithTrend(timeNum, monthIdx, rawList, colours, labels, ...
                       yLabel, slopeUnitTxt)

    % ---- figure set-up -------------------------------------------------
    figure('Color','w','Position',[120 120 900 480]); hold on;

    % 1. plot raw curves, capture handles for legend
    hRaw = gobjects(numel(rawList),1);
    for i = 1:numel(rawList)
        hRaw(i) = plot(timeNum, rawList{i}, '-', ...
                       'Color', colours{i}, 'LineWidth',0.5);
    end
    xlabel('Year');  ylabel(yLabel);
    legend(hRaw, labels,'Location','best');
    grid on;

    ax = gca;
    ax.XLim  = [timeNum(1) timeNum(end)];
    ax.XTick = timeNum(1):365*5:timeNum(end);  % ~5-year spacing
    datetick('x','yyyy','keepticks','keeplimits');
    ax.FontSize = 9;

    % 2. trend + annotation loop
    for i = 1:numel(rawList)

        y = rawList{i};
        C = colours{i};

        % Sen & M-K
        [~,pMK]  = Mann_Kendall(y);
        senSlope = Sen_Slope(y);        % raw units per month
        senDec   = senSlope*120;        % … per decade
        [~,pMKm] = Modified_Mann_Kendall(y);

        % straight trend line
        b0      = nanmedian(y) - senSlope*nanmedian(monthIdx);
        senLine = b0 + senSlope*monthIdx;
        plot(timeNum, senLine, '--', 'Color',C,'LineWidth',1.2);

        % annotation text (6 months from right edge)
        xText = timeNum(end) - 180;
        yText = senLine(end-6);
        txt   = sprintf('%.2f%s  p=%.3f', senDec, slopeUnitTxt, pMK);
        text(xText, yText, txt, 'Color',C,'FontSize',7, ...
             'Interpreter','tex', 'Horiz','left','Vert','middle');

        % console
        fprintf('%-12s  slope = %+8.4f %s,  pMK=%.3f,  pMKmod=%.3f\n', ...
                labels{i}, senDec, slopeUnitTxt, pMK, pMKm);
    end
end