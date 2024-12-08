% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/simulation_gui.m
function simulation_gui_dynamic_geometry()
    % Create the GUI figure
    figWidth = 800;
    figHeight = 700;
    fig = uifigure('Position', [100 100 figWidth figHeight], 'Name', 'Simulation Parameters');

    % Layout parameters
    startY = figHeight - 60;   % start near the top
    labelWidth = 500;
    fieldWidth = 500;
    labelX = 20;
    fieldX = 220;
    height = 22;
    vertGap = 50; % vertical gap between main fields
    explanationOffset = 25; % vertical distance below field for explanation text

    % Function to add a field with label and explanation and return handle
    function [currentY_out, fieldHandle] = addField(labelStr, defaultValue, explanationStr, fieldType, currentY)
        % label
        uilabel(fig, 'Position', [labelX currentY labelWidth height], 'Text', labelStr);

        % field
        switch fieldType
            case 'text'
                fieldHandle = uieditfield(fig, 'text', ...
                    'Position', [fieldX currentY fieldWidth height], 'Value', defaultValue);
            case 'numeric'
                fieldHandle = uieditfield(fig, 'numeric', ...
                    'Position', [fieldX currentY fieldWidth height], 'Value', str2double(defaultValue));
            case 'textarea'
                fieldHandle = uitextarea(fig, ...
                    'Position', [fieldX (currentY - 60) fieldWidth 80], 'Value', defaultValue);
                currentY = currentY - 60; 
            otherwise
                error('Unsupported field type');
        end

        % explanation
        uilabel(fig, 'Position', [fieldX (currentY - explanationOffset) fieldWidth 20], ...
            'Text', explanationStr, 'FontSize',10);

        % Adjust currentY for next field
        if strcmp(fieldType, 'textarea')
            currentY_out = currentY - vertGap - 80; % extra space for textarea
        else
            currentY_out = currentY - vertGap;
        end
    end

    currentY = startY;

    % Nucleus Obstruction Vector
    [currentY, nucRelVecEdit] = addField('Nucleus Obstruction Vector:', '[1]', ...
        'Specify the nucleus obstruction vector. (e.g., [1] as obstruction, [0] as no obstruction)', ...
        'text', currentY);

    % Number of Trials
    [currentY, numTrialsEdit] = addField('Number of Trials:', '1', ...
        'Specify how many times to repeat the simulation. (e.g., 3)', ...
        'numeric', currentY);

    % Stretch Schedule (textarea)
    [currentY, stretchScheduleEdit] = addField('Stretch Schedule:', ...
        {'0.0, 48*3600, 1.0, 1.0;'; '48*3600, 72*3601, 1.0, 2.0;'}, ...
        'Specify the start time, end time, initial width, and final width for stretching. (e.g., 0.0,48*3600,1.0,1.0)', ...
        'textarea', currentY);

    % Filename Save
    [currentY, filenameSaveEdit] = addField('Filename Save:', 'default_filename.mat', ...
        'Specify the .mat filename to save simulation results. (e.g., default_filename.mat)', ...
        'text', currentY);

    % Videoname Save
    [currentY, videonameSaveEdit] = addField('Videoname Save:', 'default_video.mp4', ...
        'Specify the .mp4 filename to save simulation video. (e.g., default_video.mp4)', ...
        'text', currentY);

    % Run Simulation button
    runButtonY = currentY - vertGap;
    uibutton(fig, 'Position', [figWidth/2-50 runButtonY 100 30], 'Text', 'Run Simulation', ...
        'ButtonPushedFcn', @(runButton,event) runSimulation());

    % Callback function to run the simulation
    function runSimulation()
        % Get parameters directly from the handles
        nuc_rel_vec = str2num(nucRelVecEdit.Value); %#ok<ST2NM>
        num_trials = numTrialsEdit.Value;
        stretch_schedule_str = stretchScheduleEdit.Value;
        stretch_schedule = eval(['[', strjoin(stretch_schedule_str', ';'), ']']);

        params.nuc_rel_vec = nuc_rel_vec;
        params.stretch_schedule = stretch_schedule;

        filename_base = filenameSaveEdit.Value;
        videoname_base = videonameSaveEdit.Value;

        % Run simulation
        for k = 1:num_trials
            params.filename_save = replace(filename_base, '.mat', ['_trial=', num2str(k), '.mat']);
            params.videoname_save = replace(videoname_base, '.mp4', ['_trial=', num2str(k), '.mp4']);
            single_cell_units_linked_v3_DynamicGeometry(params);
        end
        uialert(fig, 'Simulation completed.', 'Success');
    end
end
