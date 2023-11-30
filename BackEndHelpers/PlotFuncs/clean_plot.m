function clean_plot(f,varargin)
% clean_plot does basic figure clean-up
%
% It can be called directly after plotting
% e.g. 
%     >>  figure; hold on;
%     >>  plot(randn(10,2),randn(10,2),'.'); 
%     >>  clean_plot; 
%
% INPUTS
%     f (optional) : figure handle
%
% clean_plot.m adjusts a few properties
%     1. Makes figure larger and sets background color to 'w'
%     2. Sets axes tick directions to 'out', makes them larger, and gray, and
%        increases font size
%     3. Sets marker size ('.' to 15, 'o' and others to 5); 
%     4. Sets line width to 1
%     5. Sets EdgeColor to 'none' (for bar or patch) and FaceAlpha to 0.5
%        (for e.g. patch)
%     6. Adjusts legend if it exists
%%
if nargin==0 || isempty(f)
    f = gcf; % get figure handle
end
set(gcf,'Renderer','painters');
if ~isempty(varargin)
    if any(cellfun(@(x) contains(x,'time','IgnoreCase',true), varargin))
        timeplot = true; 
    else
        timeplot = false; 
    end
    if any(cellfun(@(x) contains(x,'timebar','IgnoreCase',true), varargin))
        timebar = true; 
        tbar_i = cellfun(@(x) contains(x,'timebar','IgnoreCase',true), varargin); 
        tstr = varargin{tbar_i}; 
        extras = length(tstr)>7; 
        if extras
            extra_i = strfind(tstr,'_')+1; 
            extrastr = tstr(extra_i:end); 
            
            num_text = cellfun(@(x) ismember(str2double(x),0:9),num2cell(extrastr)); 
            bar_val = str2double(extrastr(num_text)); 
            bar_unit = extrastr(~num_text); 
        else
            bar_val = 1; 
            bar_unit = 's'; 
        end
    else
        timebar = false; 
    end
else
    timeplot = false; 
    timebar = false; 
end

% 1: Resize 
if all(f.Position(3:4)==[560 420])
    f.Position(1:2) = f.Position(1:2).*.8; 
    f.Position(3:4) = [560 420]*1.25; 
end
set(f,'Color','w'); 

is_axishandle = contains(char(class(f)),'Axes'); 

if ~is_axishandle
    subplot_num = length(f.Children); % Number of axes/subplots
else
    subplot_num=1; 
end
for q = 1:subplot_num % loop through axes
    
    if is_axishandle
        continue_condition = 1; 
    else
        continue_condition = contains(char(class(f.Children(q))),'Axes')*2; 
    end
    
    if ismember(continue_condition,[1 2])
        
        if continue_condition == 1
            h = f; 
        elseif continue_condition ==2 
            h = f.Children(q); % axis handle
        end

        % 2: Fix Ticks and Axes
        set(h,'TickDir','out'); 
        box(h,'off'); 
        set(h,'TickLength',[.02 .02]); 
        set(h,'FontSize',14); 
        set(h,'LabelFontSizeMultiplier',1.25); 
        set(h,'XColor',[0.25 0.25 0.25]); 
        set(h,'YColor',[0.25 0.25 0.25]); 

        % 3: Set MarkerSize
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'MarkerSize') 
                if strcmp(h.Children(i).Marker,'.') && h.Children(i).MarkerSize==6
                    h.Children(i).MarkerSize = 15; 
                elseif strcmp(h.Children(i).Marker,'o') && h.Children(i).MarkerSize==6 
                    h.Children(i).MarkerSize = 5; 
                end
            end
        end
        % 4: Set LineWidth
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'LineWidth') && h.Children(i).LineWidth==0.5
                h.Children(i).LineWidth = 1; 
            end
        end

        % 5: Set EdgeColor and FaceAlpha
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'EdgeColor') 
                has_markerface = isprop(h.Children(i),'MarkerFaceColor') && ~strcmp(h.Children(i).MarkerFaceColor,'none'); 
                has_face = isprop(h.Children(i),'FaceColor') && ~(strcmp(h.Children(i).FaceColor,'none') || strcmp(char(h.Children(i).FaceColor),char([0 0 0]))); 
                if has_markerface||has_face
                    h.Children(i).EdgeColor = 'none'; 
                end
            end
        end
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'FaceAlpha') 
                if strcmp(get(h.Children(i),'Type'),'bar')
                    h.Children(i).FaceAlpha = 0.9;
                elseif h.Children(i).FaceAlpha==1
                    h.Children(i).FaceAlpha = 0.5;
                end
            end
        end

        if timeplot || timebar
            for i = 1:length(h.Children)
               if isprop(h.Children(i),'XData') 
                    h.Children(i).XData = h.Children(i).XData*.02; 
               end
            end
            [ylims,xlims] = deal(NaN(length(h.Children),2)); 
            for i = 1:length(h.Children)
                if isprop(h.Children(i),'XData') 
                    xlims(i,:) = [min(h.Children(i).XData),max(h.Children(i).XData)];  
                    ylims(i,:) = [min(h.Children(i).YData),max(h.Children(i).YData)];  
                end
            end
            set(h,'XLim',[min(xlims(:,1)),max(xlims(:,2))]) ;
            h.XLabel.String = 'Time (s)'; 
            

            if timebar
                xl = [min(xlims(:,1)),max(xlims(:,2))]; 
                offset = range(xl)*.025; 
                for i = 1:length(h.Children)
                    if isprop(h.Children(i),'XData')
                        h.Children(i).XData = h.Children(i).XData+offset; 
                    end
                end
                if strcmp(bar_unit,'s') || strcmp(bar_unit,'sec')
                    bbound = [0 bar_val]; 
                elseif strcmp(bar_unit,'m') || strcmp(bar_unit,'min')
                    bbound = [0 bar_val*50]; 
                elseif strcmp(bar_unit,'h') || strcmp(bar_unit,'hrs')
                    bbound = [0 bar_val*50*60]; 
                end
                tb = plot(h,bbound+offset,min(ylims(:))*[1 1] + .5*range(ylims(:)),'Color',[.5 .5 .5],'LineWidth',3) ;
                tb.HandleVisibility = 'off'; 
                h.XAxis.Visible = 'off'; 
                text(h,0+offset,min(ylims(:)) + .5*range(ylims(:)),[num2str(bar_val) '\_' bar_unit],'FontSize',12,'Color',[.5 .5 .5],...
                    'VerticalAlignment','bottom','HorizontalAlignment','left'); 
            end
                
            
        end
        % 6: Check if legend exists and clean it up
        if ~isempty(h.Legend)
            leg = h.Legend; 
            set(leg,'Location','best'); 
            set(leg,'Box','off'); 
            set(leg,'FontSize',14); 
        end

    end

end

