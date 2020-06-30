% this code creates the coop plots in fig2_v2

% load COOP results for specified S value
[file,path_p]=uigetfile({'*.mat';'*.*'},['Select which COOP .mat result files to use '],...
        path,'MultiSelect','on');

if iscell(file) %Determine if the user wants to go through one file or multiple
    NumFiles = length(file);
    for i=1:NumFiles
        temp = file{i};
        FileIndex{i} = temp;
    end
else
    NumFiles = 1;
    FileIndex{1} = file;
end

for count=1:NumFiles
    
    filename_temp_p = FileIndex{count};
    path_and_filename_p = [path_p filename_temp_p];
    filename_p = path_and_filename_p;
    
    %load parameters
    load(filename_p,'S','COOP')

    % Produce the input matrix data
    C_full=COOP;

    % Set [min,max] value of C to scale colors
    clrLim = [0,1]; 

    % Set the  [min,max] of diameter where 1 consumes entire grid square
    diamLim = [1,1];

    % Compute center of each circle 
    % This assumes the x and y values were not entered in imagesc()
    x = 1 : 1 : size(C_full,2); % x edges
    y = 1 : 1 : size(C_full,1); % y edges
    [xAll, yAll] = meshgrid(x,y); 

    % Set color of each rectangle
    cmap = parula(256); 
    Cscaled = (C_full - clrLim(1))/range(clrLim); % always [0:1]
    colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
    colIdx(1)=colIdx(6);

    % Set size of each circle
    diamSize = Cscaled * range(diamLim) + diamLim(1); 
%     diamSize = triu(Cscaled * range(diamLim) + diamLim(1)); % uncomment this instead of above line if using C_mm or C_ee

    % Create figure
    fh = figure(); 
    ax = axes(fh); 
    hold(ax,'on')
    colormap(ax,'parula');

    % Create circles
    theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
    h = arrayfun(@(k)fill(diamSize(k)/2 * cos(theta) + xAll(k), ...
        diamSize(k)/2 * sin(theta) + yAll(k), cmap(colIdx(k),:)),1:numel(xAll)); 
    
    axis(ax,'equal')
    axis(ax,'tight')
    set(ax,'YDir','Reverse')
    colorbar()
    caxis(clrLim);
    title({FileIndex{count},[ ' S=',num2str(S)]})

end

