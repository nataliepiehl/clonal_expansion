%simulation of stem cell clonal expansion
% with memory included
%written by AsahiNakamuta out of Kyoto University
%lightly annotated NP 10/11/2023

tic
clear variables
clc

% Find all timeseries matrices
mats = dir(fullfile("/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/matrices/", "timeseries*.mat"));

% For each matrix...
for i=1:length(mats)
    % Define filename
    mat_filename = mats(i).name;

    % Load matrix
    x_matrix = importdata("/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/matrices/" + mat_filename);

    % Load in variables
    load(replace(strcat("/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/matrices/" + mat_filename), 'timeseries', 'variables'))
    
    % Generate clone size over time plot for each clone
    tmax = 250; %changed to 250 so each simulation has data for all time points (400 open niche case)
    fig = figure('visible','off');
    fig.Position = [10 10 500 800]; 
    for k=1:num_of_clones
        subplot(num_of_clones,1,k)
        plot(x_matrix(k,:))
        xlim([0 tmax*(lambda*n_openniche+epsilon*num_of_clones)])
        ylim([0,n_openniche])
        title(['clone ',num2str(k)])
        if k ~= num_of_clones
            set(gca,'xtick',[])
        end
    end

     % Add title
    sgtitle(strcat("\epsilon = ", num2str(epsilon), ", \lambda = ", num2str(lambda), ", \alpha = ", num2str(memory_strength)))

    % Give common xlabel and ylabel
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Clone size');
    xlabel(han,'time');

    % Save image (clone size across time)
    % fontsize(16, "points")
    fontsize(fig, 16, "points")
    plot_filename = "/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/plots/" + erase(mat_filename, ".mat") + "_size.png";
    exportgraphics(fig, plot_filename, 'Resolution', 300)

    % Generate plot of proportion of each clone over time
    fig = figure('visible','off');
    fig.Position = [10 10 500 400]; 
    b = bar(x_matrix',1,'stacked','FaceColor','flat');
    xlim([0 tmax*(lambda*n_openniche+epsilon*num_of_clones)])
    ax = gca;
    ax.XAxis.Exponent = 4;

    % Add title
    sgtitle(strcat("\epsilon = ", num2str(epsilon), ", \lambda = ", num2str(lambda), ", \alpha = ", num2str(memory_strength)))

    % Add color map (Renoir colors from MetBrewer)
    color_map = [
        25 18 82;
        112 90 163;
        158 154 218;
        255 176 174;
        243 123 108;
        207 37 26;
        242 153 0;
        255 186 49;
        174 165 23;
        40 90 32
        ];
    color_map = color_map / 255;
    for clone=1:10
        b(clone).CData = color_map(clone,:);
    end

    % Give common xlabel and ylabel
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Clone size (stacked)');
    xlabel(han,'time');

    % Save image (clone proportion across time)
    % fontsize(16, "points")
    fontsize(fig, 16, "points")
    plot_filename = "/projects/p31666/nat/clonal_expansion/results/hierarchal_neutral_competition/memory/plots/" + erase(mat_filename, ".mat") + "_proportion.png";
    exportgraphics(fig, plot_filename, 'Resolution', 300)
end

toc