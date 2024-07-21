function [varargout] = actograph(circdata,options)
%% actograph | a paltry actogram plotter for circadian rhythm data
arguments
    circdata;
    options.Ymax = max(circdata,[],'all');
    options.titleName = {'Locomotion Heatmap', 'Locomotion Trendline', 'Actogram'};
    options.Heatmap = true;
    options.ColorPatch = {};
    options.ReverseLastRow = false;
end

if ~options.ReverseLastRow
    % if the number of nan in the last row is more than half of the row length
    % then remove the last row
    if sum(isnan(circdata(end,:))) > size(circdata,2)/2
        circdata = circdata(1:end-1,:);
    end
    len_data = size(circdata,1);
    if ~isempty(options.ColorPatch)
        for i = 1:length(options.ColorPatch)
            if size(options.ColorPatch{i}.PatchMap,1) > len_data
                options.ColorPatch{i}.PatchMap = options.ColorPatch{i}.PatchMap(1:len_data,:);
            end
        end
    end

%% --| PLOT COLORMAP AND BAR GRAPH WITH SMOOTHED SPLINE OVERLAY
if options.Heatmap
    fig1 = figure;
    set(fig1,'OuterPosition',[200 350 800 500],'Color',[1 1 1]);
    hax1 = axes('Position',[.05 .07 .4 .85],'Color','none');
    hax2 = axes('Position',[.55 .07 .4 .85],'Color','none');
    colormap(hax1,'hot'); colormap(hax2,'bone');
    axes(hax1); 
    imagesc(circdata); 
    title(options.titleName{1});
    % axis off;
    
    % smeth = {'moving','lowess','loess','sgolay','rlowess','rloess'};
    % degSm = .1; typSm = 5; muCd = mean(circdata);
    % smUnc = smooth(muCd,degSm,smeth{typSm});
    smUnc = mean(circdata, 1, 'omitnan');
    
    axes(hax2); 
    bar(circdata',1.0); hold on;
    ph = plot(smUnc);
    set(ph,'Color',[.8 .1 .1],'LineWidth',4,'Marker','none'); hold off
    title(options.titleName{2});
    % axis off; 
    FigsOut.heatmap = fig1;
end



%% --| PLOT ACTOGRAM ACTOGRAPH
fig2 = figure;
set(fig2,'OuterPosition',[1000 50 450 800],'Color',[1 1 1]);
rowNum = size(circdata,1);
axPos = fliplr(linspace(.02, 0.97-0.98/rowNum, rowNum));
% rowh = axPos(1) - axPos(2) - .005;%default
rowh = axPos(1) - axPos(2) - .01; %LY
Ymax = options.Ymax;
if isnan(Ymax) || Ymax <= 0
    Ymax = 100;
end

for row_i = 1:rowNum
    axes('Position',[.05 axPos(row_i) .90 rowh],'Color','none','XTick',[],'YTick',[]); axis off; hold on
    if ~isempty(options.ColorPatch)
        for i = 1:length(options.ColorPatch)
            c = options.ColorPatch{i}.Color;
            if row_i > size(options.ColorPatch{i}.PatchMap,1)
                continue;
            end
            pMat = options.ColorPatch{i}.PatchMap(row_i,:);
%             y = [0,0,Ymax,Ymax]; % all patch 
            if i==1 
                y = [0,0,Ymax,Ymax];
            else % band patch
                y1 = 1.1 - i/10;
                y2 = 1.2 - i/10;
                y = [y1,y1,y2,y2] * Ymax;
            end
            for p = 1:length(pMat)
                x = [(p-1),p,p,(p-1)] + 0.5;
                if ~isnan(pMat(p))
                    patch(x,y,c,'EdgeColor','none','FaceAlpha',pMat(p));
                end
            end
        end
    end
    bar(circdata(row_i,:),1.0,'EdgeAlpha',0,'FaceColor',[0 0 0]);
    set(gca,'YLim',[0 Ymax]);
    if row_i == 1; title(options.titleName{3}+" (Ylim="+Ymax+")"); end
end
FigsOut.actogram = fig2;

varargout = {FigsOut, circdata, options.ColorPatch};

end