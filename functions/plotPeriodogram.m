function fig = plotPeriodogram(P_candi,powerArray,confidenceLine,confidence,options)
arguments
    P_candi;
    powerArray;
    confidenceLine;
    confidence;
    options.Visible = 'off'
    options.titleName = 'Periodogram';
    options.xlabel = 'Period/hour';
    options.ylabel = 'Statistic'
end
[m, I] = max(powerArray - confidenceLine);
fig = figure; hold on
fig.Visible = options.Visible;
plot(P_candi,powerArray,'b');
plot(P_candi,confidenceLine,'r');
if m > 0
    scatter(P_candi(I),powerArray(I),50,'*','MarkerEdgeColor','#A2142F')
end
xlabel(options.xlabel); 
ylabel(options.ylabel);
title(options.titleName);
xlim([P_candi(1) P_candi(end)]);
yceil = RoundByMagnitude(min(confidenceLine)*5);
if m+confidenceLine(I)+20 < yceil || isnan(m)
    ylim([0 yceil]);
else
    ylim([0 RoundByMagnitude(m+confidenceLine(I)+50)]);
end
textLabel = num2str(P_candi(I));
text(P_candi(I),powerArray(I)+10,textLabel,'FontSize',12)
legend({'Power',['CI: ' num2str(confidence,8)]});

end