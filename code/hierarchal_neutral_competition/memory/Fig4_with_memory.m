%scaling law in clone size
%now with memory
%written by AsahiNakamuta out of Kyoto University
%lightly annotated NP 10/11/2023

tic
clear variables
clc

% Find all timeseries matrices
mats = dir(fullfile("results/hierarchal_neutral_competition/memory/matrices/", "clone1_timeseries*.mat"));

% For each matrix...
for mat=1:length(mats)
    % Define filename
    mat_filename = mats(mat).name;

    % Load matrix
    x_matrix = importdata("results/hierarchal_neutral_competition/memory/matrices/" + mat_filename);

    % Load in variables
    load(replace(strcat("results/hierarchal_neutral_competition/memory/matrices/" + mat_filename), 'timeseries', 'variables'))

    % Initialize the matrices to store the clone sizes per time point
    temp = zeros(3,n_openniche);
    unscale = zeros(3,n_openniche);
    scale = zeros(3,n_openniche);
    x_axis = zeros(3,n_openniche);
    for t = 1:3
        % 5 because we set tmax to 15 in the simulation
        pickup = round(5*t*(epsilon*num_of_clones+lambda*n_openniche));
        % For each simulation
        for i = 1:m
            if x_clone1(i,pickup)~=0
                temp(t,x_clone1(i,pickup)) = temp(t,x_clone1(i,pickup))+1;
            end
        end
        % For each possible clone size...
        for j=1:n_openniche
            % Divide number in each clone size by number of clones for
            % proportion density
            unscale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup));
            % Multiply by the average clone size
            scale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup))*sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup));
            % Divide by the average clone size
            x_axis(t,j) = j/(sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup)));
        end
    end

    % Define color map (Renoir colors from MetBrewer)
    color_map = [
        25 18 82;
        255 176 174;
        174 165 23;
        ];
    color_map = color_map / 255;
    
    % Produce unscaled figure
    fig = figure('visible','off');
    fig.Position = [10 10 500 400]; 
    semilogy(0,0)
    hold on
    for t = 1:3
        p = semilogy(unscale(t,:), "LineWidth", 2);
        p.Color = color_map(t,:);
    end
    hold off

    % Add formatting
    sgtitle(strcat("\epsilon = ", num2str(epsilon), ", \lambda = ", num2str(lambda), ": unscaled"))
    legend("", "time = 5", "time = 10", "time = 15")
    ylabel('Probability density');
    xlabel('Clone size');

    % Save image (clone size across time)
    fontsize(20, "points")
    plot_filename = "results/hierarchal_neutral_competition/memory/plots/" + erase(mat_filename, ".mat") + "_unscaled.png";
    exportgraphics(fig, plot_filename, 'Resolution', 300)

    % Produce scaled figure
    fig = figure('visible','off');
    fig.Position = [10 10 500 400]; 
    semilogy(0,0)
    hold on
    for t = 1:3
        p = semilogy(x_axis(t,:),scale(t,:), "LineWidth", 2);
        p.Color = color_map(t,:);
    end
    hold off

    % Add formatting
    sgtitle(strcat("\epsilon = ", num2str(epsilon), ", \lambda = ", num2str(lambda), ": scaled"))
    legend("", "time = 5", "time = 10", "time = 15")
    ylabel('Probability density x avg. clone size');
    xlabel('Clone size / avg. clone size');

    % Save image (clone size across time)
    fontsize(20, "points")
    plot_filename = "results/hierarchal_neutral_competition/memory/plots/" + erase(mat_filename, ".mat") + "_scaled.png";
    exportgraphics(fig, plot_filename, 'Resolution', 300)
end

toc