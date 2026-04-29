classdef FASC_GUI < matlab.apps.AppBase

    % =====================================================================
    % Properties: App Components
    % =====================================================================
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        MainGrid             matlab.ui.container.GridLayout

        % --- Left Sidebar (Controls) ---
        SidePanel            matlab.ui.container.Panel
        SidebarGrid          matlab.ui.container.GridLayout

        % Header
        HeaderPanel          matlab.ui.container.Panel
        AppTitle             matlab.ui.control.Label
        AppSubtitle          matlab.ui.control.Label

        % Import/Save Area
        ImportGrid           matlab.ui.container.GridLayout
        BtnLoad              matlab.ui.control.Button
        BtnSave              matlab.ui.control.Button

        % Parameters Area
        ScrollContainer      matlab.ui.container.Panel
        ScrollLayout         matlab.ui.container.GridLayout

        % Dynamic Components
        LblIdxPos            matlab.ui.control.Label
        EditIdxPos           matlab.ui.control.EditField
        LblIdxNeg            matlab.ui.control.Label
        EditIdxNeg           matlab.ui.control.EditField

        % Inputs
        EditSimInner         matlab.ui.control.NumericEditField
        EditSimInter         matlab.ui.control.NumericEditField
        EditInitLim          matlab.ui.control.NumericEditField
        EditMaxClust         matlab.ui.control.NumericEditField
        EditMaxIter          matlab.ui.control.NumericEditField
        EditMinVol           matlab.ui.control.NumericEditField
        DropStrategy         matlab.ui.control.DropDown
        DropAlgorithm        matlab.ui.control.DropDown

        BtnRun               matlab.ui.control.Button

        % Status
        StatusPanel          matlab.ui.container.Panel
        StatusLamp           matlab.ui.control.Lamp
        StatusLabel          matlab.ui.control.Label

        % --- Right Content (Tabs) ---
        ContentPanel         matlab.ui.container.Panel
        TabGroup             matlab.ui.container.TabGroup
        TabConvergence       matlab.ui.container.Tab
        TabHeatmap           matlab.ui.container.Tab
        TabDistribution      matlab.ui.container.Tab

        % Containers for embedded plots
        PanelConv            matlab.ui.container.Panel
        PanelOutliers        matlab.ui.container.Panel
        PanelHeatmap         matlab.ui.container.Panel
        PanelDist            matlab.ui.container.Panel

        % --- Bottom (Log) ---
        LogPanel             matlab.ui.container.Panel
        LogText              matlab.ui.control.TextArea
    end

    % =====================================================================
    % Properties: Visual Theme & Data
    % =====================================================================
    properties (Access = private)
        % --- Modern Palette ---
        ColorBgFigure  = [0.94 0.94 0.96]; 
        ColorBgPanel   = [0.94 0.94 0.96]; 
        ColorInputBg   = [0.22 0.25 0.30]; 

        ColorLabel     = [0.20 0.20 0.20]; 
        ColorTextMain  = [0.95 0.95 0.95]; 
        ColorTextDim   = [0.50 0.50 0.50]; 

        ColorAccent    = [0.00 0.48 1.00]; 
        ColorSuccess   = [0.16 0.75 0.45]; 
        ColorWarning   = [1.00 0.75 0.00]; 
        ColorError     = [0.90 0.25 0.25]; 

        ColorContentBg = [0.98 0.98 1.00]; 

        DataMatrix
        Results
        FutureObj
        LogTimer
        LogFileName

        % Rows for dynamic hiding
        RowIdxPos
        RowIdxNeg
    end

    % =====================================================================
    % Methods: Constructor (Entry Point)
    % =====================================================================
    methods (Access = public)
        function app = FASC_GUI
            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app, app.UIFigure);

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.UIFigure);
        end
    end

    % =====================================================================
    % Methods: Application Logic
    % =====================================================================
    methods (Access = private)

        function runAnalysis(app)
            % 1. Validation
            if isempty(app.DataMatrix)
                app.log("Error: No data loaded.", 'err');
                uialert(app.UIFigure, "Please load data (MAT or CSV) first.", "No Data");
                return;
            end

            % 2. Get Parameters
            try
                % Only evaluate Idx if visible/relevant, otherwise default (safe fallback)
                if strcmp(app.DropAlgorithm.Value, 'dual-cosine')
                    idx_pos = eval(app.EditIdxPos.Value);
                    idx_neg = eval(app.EditIdxNeg.Value);
                else
                    idx_pos = [];
                    idx_neg = [];
                end

                p.simInner = app.EditSimInner.Value;
                p.simInter = app.EditSimInter.Value;
                p.initLim = app.EditInitLim.Value;
                p.maxClust = app.EditMaxClust.Value;
                p.maxIter = app.EditMaxIter.Value;
                p.minVol = app.EditMinVol.Value;
                p.strat = app.DropStrategy.Value;
                p.algo = app.DropAlgorithm.Value;
            catch
                app.log("Error: Invalid parameter input syntax.", 'err');
                uialert(app.UIFigure, "Check syntax for Index Positive/Negative.", "Input Error");
                return;
            end

            app.Results.Params = p;
            app.Results.idx_pos = idx_pos;
            app.Results.idx_neg = idx_neg;

            % 3. UI State -> Running
            app.toggleUI('off');
            app.setStatus('running');
            app.log("------------------------------------------");
            app.log("System: Initializing Background Worker...");

            app.LogFileName = [tempname '.txt'];
            fid = fopen(app.LogFileName, 'w'); fclose(fid);

            % 4. Submit to Parallel Pool
            try
                pool = gcp('nocreate');
                if isempty(pool)
                    app.log("System: Starting Parallel Pool (One-time setup)...");
                end

                args = {app.DataMatrix, idx_pos, idx_neg, ...
                    p.simInter, p.simInner, p.initLim, p.maxClust, ...
                    p.maxIter, p.strat, p.minVol, p.algo};

                app.FutureObj = parfeval(@FASC_GUI.runFASC_Wrapper, 4, ...
                    app.LogFileName, args);

                % Setup Timer
                if ~isempty(app.LogTimer) && isvalid(app.LogTimer)
                    stop(app.LogTimer); delete(app.LogTimer);
                end
                app.LogTimer = timer('ExecutionMode', 'fixedRate', ...
                    'Period', 0.5, ...
                    'TimerFcn', @(~,~) app.updateLogFromTimer());
                start(app.LogTimer);

                afterEach(app.FutureObj, @(~) app.checkStatus(), 0);

            catch ME
                app.log("CRITICAL ERROR: " + ME.message, 'err');
                app.toggleUI('on');
                app.setStatus('error');
            end
        end

        function updateLogFromTimer(app)
            if exist(app.LogFileName, 'file')
                try
                    txt = fileread(app.LogFileName);
                    if strlength(txt) > 0
                        lines = split(string(txt), newline);
                        lines = lines(strlength(lines)>0);
                        app.LogText.Value = lines;
                        scroll(app.LogText, 'bottom');
                    end
                catch
                end
            end

            if ~isempty(app.FutureObj) && strcmp(app.FutureObj.State, 'finished')
                stop(app.LogTimer);
                app.finalizeResults();
            end
        end

        function checkStatus(app)
            % Placeholder
        end

        function finalizeResults(app)
            try
                [cCenters, cCounts, cIdx, iInfo] = fetchOutputs(app.FutureObj);

                app.Results.cCenters = cCenters;
                app.Results.cCounts = cCounts;
                app.Results.cIdx = cIdx;
                app.Results.iInfo = iInfo;

                app.log("System: Processing visualizations...");

                if exist(app.LogFileName, 'file')
                    delete(app.LogFileName);
                end

                app.embedPlots();
                app.log("System: Analysis Completed Successfully.");
                app.setStatus('ready');

                % Enable Save Button
                app.BtnSave.Enable = 'on';

            catch ME
                app.log("Error in Execution: " + ME.message, 'err');
                if ~isempty(app.FutureObj.Error)
                    app.log("Worker Trace: " + app.FutureObj.Error.message, 'err');
                end
                app.setStatus('error');
            end
            app.toggleUI('on');
        end

        function saveData(app)
            % Check if results exist
            if isempty(app.Results) || ~isfield(app.Results, 'cCenters')
                uialert(app.UIFigure, "No clustering results to save.", "Error");
                return;
            end

            % Open File Dialog
            [file, path] = uiputfile('*.mat', 'Save Clustering Results');
            if file == 0
                return;
            end

            app.log("System: Saving data and exporting plots...");
            app.toggleUI('off');
            drawnow;

            % 临时关闭图形可见性
            prevVis = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');

            try
                % 1. Save .mat file
                saveStruct.cC = app.Results.cCenters;
                saveStruct.cN = app.Results.cCounts;
                saveStruct.ptcIdx = app.Results.cIdx;
                saveStruct.iInfo = app.Results.iInfo;

                saveFull = fullfile(path, file);
                save(saveFull, '-struct', 'saveStruct');

                % 2. Export Plots as SVG
                [~, baseName, ~] = fileparts(file);

                % Helper 1: Convergence
                Clustering_iterInfoAnalyzer(app.Results.iInfo, 0, 1, baseName, path);

                % Helper 2: Heatmap
                Clustering_clusterListCountsHister(...
                    app.Results.cCounts, app.Results.cCenters, ...
                    app.Results.Params.simInter, app.Results.Params.algo, ...
                    app.Results.idx_pos, app.Results.idx_neg, ...
                    0, 1, baseName, path);

                set(0, 'DefaultFigureVisible', prevVis);

                app.log("System: Saved successfully to " + path);
                uialert(app.UIFigure, "Data and plots saved successfully.", "Success");

            catch ME
                set(0, 'DefaultFigureVisible', prevVis);
                app.log("Error saving data: " + ME.message, 'err');
                uialert(app.UIFigure, "Error saving data. See log.", "Save Error");
            end

            app.toggleUI('on');
        end

        function embedPlots(app)
            % 1. Convergence
            figs = app.captureFigures(@() Clustering_iterInfoAnalyzer(...
                app.Results.iInfo, 1, 0, "", ""));

            if length(figs) >= 1, app.reparentFigure(figs(1), app.PanelConv); end
            if length(figs) >= 2, app.reparentFigure(figs(2), app.PanelOutliers); end

            % 2. Heatmap
            figs = app.captureFigures(@() Clustering_clusterListCountsHister(...
                app.Results.cCounts, app.Results.cCenters, ...
                app.Results.Params.simInter, app.Results.Params.algo, ...
                app.Results.idx_pos, app.Results.idx_neg, ...
                1, 0, "", ""));

            if length(figs) >= 1, app.reparentFigure(figs(1), app.PanelHeatmap); end
            if length(figs) >= 2, app.reparentFigure(figs(2), app.PanelDist); end
        end

        function updateLayout(app)
            % Dynamic layout adjustment based on Algorithm selection
            algo = app.DropAlgorithm.Value;
            isDual = strcmp(algo, 'dual-cosine');

            rPos = app.RowIdxPos;
            rNeg = app.RowIdxNeg;

            % Safety check in case layout hasn't initialized fully
            if isempty(rPos) || isempty(rNeg), return; end

            if isDual
                % Show Index inputs
                app.LblIdxPos.Visible = 'on'; app.EditIdxPos.Visible = 'on';
                app.LblIdxNeg.Visible = 'on'; app.EditIdxNeg.Visible = 'on';
                % Restore height (UPDATED to 24)
                app.ScrollLayout.RowHeight{rPos} = 24;
                app.ScrollLayout.RowHeight{rNeg} = 24;
            else
                % Hide Index inputs
                app.LblIdxPos.Visible = 'off'; app.EditIdxPos.Visible = 'off';
                app.LblIdxNeg.Visible = 'off'; app.EditIdxNeg.Visible = 'off';
                % Collapse height to 0 so other elements slide up
                app.ScrollLayout.RowHeight{rPos} = 0;
                app.ScrollLayout.RowHeight{rNeg} = 0;
            end
        end

        function reparentFigure(~, figHandle, targetPanel)
            delete(targetPanel.Children);
            axesHandles = findall(figHandle, 'Type', 'axes');
            colorbars = findall(figHandle, 'Type', 'colorbar');

            if isempty(axesHandles)
                close(figHandle); return;
            end

            set(axesHandles, 'Parent', targetPanel);
            set(axesHandles, 'Units', 'normalized');
            set(axesHandles, 'OuterPosition', [0.05 0.05 0.9 0.9]);

            if ~isempty(colorbars)
                set(colorbars, 'Parent', targetPanel);
                set(colorbars, 'Units', 'normalized');
            end

            close(figHandle);
        end

        function newFigs = captureFigures(~, funcHandle)
            existingFigs = findall(0, 'Type', 'figure');
            funcHandle();
            allFigs = findall(0, 'Type', 'figure');
            newFigs = setdiff(allFigs, existingFigs);
            [~, idx] = sort([newFigs.Number]);
            newFigs = newFigs(idx);
        end

        function log(app, msg, type)
            if nargin < 3, type = 'info'; end
            ts = string(datestr(now, 'HH:MM:SS'));
            prefix = " ";
            if strcmp(type, 'err'), prefix = " [ERR] "; end

            app.LogText.Value = [app.LogText.Value; " [" + ts + "]" + prefix + msg];
            scroll(app.LogText, 'bottom');
        end

        function toggleUI(app, state)
            isOn = strcmp(state, 'on');
            app.BtnRun.Enable = state;
            app.BtnLoad.Enable = state;

            % Handle Save Button
            if isOn
                if ~isempty(app.Results) && isfield(app.Results, 'cCenters')
                    app.BtnSave.Enable = 'on';
                else
                    app.BtnSave.Enable = 'off';
                end
            else
                app.BtnSave.Enable = 'off';
            end

            if isOn
                app.ScrollContainer.Enable = 'on';
                app.BtnRun.BackgroundColor = app.ColorAccent;
                app.BtnRun.Text = 'RUN ANALYSIS';
            else
                app.ScrollContainer.Enable = 'off';
                app.BtnRun.BackgroundColor = [0.4 0.4 0.4];
                app.BtnRun.Text = 'PROCESSING...';
            end
        end

        function setStatus(app, state)
            switch state
                case 'ready'
                    app.StatusLamp.Color = app.ColorSuccess;
                    app.StatusLabel.Text = "READY";
                    app.StatusLabel.FontColor = app.ColorSuccess;
                case 'running'
                    app.StatusLamp.Color = app.ColorWarning;
                    app.StatusLabel.Text = "PROCESSING";
                    app.StatusLabel.FontColor = app.ColorWarning;
                case 'error'
                    app.StatusLamp.Color = app.ColorError;
                    app.StatusLabel.Text = "ERROR";
                    app.StatusLabel.FontColor = app.ColorError;
            end
        end

        function loadData(app)
            filter = {'*.mat;*.csv;*.txt', 'Data Files (*.mat, *.csv, *.txt)'};
            [f, p] = uigetfile(filter);

            if f == 0, return; end

            fullPath = fullfile(p,f);
            [~, ~, ext] = fileparts(f);

            app.log("IO: Reading " + f + "...");
            app.toggleUI('off');
            drawnow;

            try
                if strcmpi(ext, '.mat')
                    d = load(fullPath);
                    fns = fieldnames(d);
                    maxSize = 0; selectedVar = '';
                    for i=1:length(fns)
                        if isnumeric(d.(fns{i})) && numel(d.(fns{i})) > maxSize
                            maxSize = numel(d.(fns{i}));
                            selectedVar = fns{i};
                        end
                    end
                    if isempty(selectedVar)
                        error("No numeric variables found in .mat file.");
                    end
                    app.DataMatrix = single(d.(selectedVar));
                    msg = sprintf("Loaded MAT: %s [%dx%d]", f, size(app.DataMatrix,1), size(app.DataMatrix,2));

                else
                    rawData = readmatrix(fullPath);
                    if isempty(rawData)
                        error("File empty or unreadable.");
                    end
                    app.DataMatrix = single(rawData);
                    msg = sprintf("Loaded CSV: %s [%dx%d]", f, size(app.DataMatrix,1), size(app.DataMatrix,2));
                end

                app.log(msg);
                app.StatusLabel.Text = "DATA LOADED";

            catch ME
                app.log("Error loading data: " + ME.message, 'err');
                uialert(app.UIFigure, "Failed to load file. See log.", "Load Error");
            end
            app.toggleUI('on');
        end

        % Helper to style inputs uniformly
        function styleInput(app, component)
            component.BackgroundColor = app.ColorInputBg;
            component.FontColor = app.ColorTextMain;
            if isprop(component, 'PlaceholderColor')
                component.PlaceholderColor = app.ColorTextDim;
            end
        end
    end

    % =====================================================================
    % Static Methods
    % =====================================================================
    methods (Static)
        function [cC, cN, ptcIdx, iInfo] = runFASC_Wrapper(logFile, args)
            if exist(logFile, 'file'), delete(logFile); end
            diary(logFile);
            try
                [cC, cN, ptcIdx, iInfo] = FASC(args{:});
            catch ME
                fprintf('\nCRITICAL ERROR: %s\n', ME.message);
                diary off;
                rethrow(ME);
            end
            diary off;
        end
    end

    % =====================================================================
    % Methods: UI Construction
    % =====================================================================
    methods (Access = private)

        function createComponents(app)
            % --- Main Window ---
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1280 800];
            app.UIFigure.Name = 'FASC | v1.0b';
            app.UIFigure.Color = app.ColorBgFigure;
            app.UIFigure.CloseRequestFcn = @(s,e) app.delete();

            % --- Main Grid Layout ---
            app.MainGrid = uigridlayout(app.UIFigure, [2, 2]);
            app.MainGrid.ColumnWidth = {320, '1x'};
            app.MainGrid.RowHeight = {'1x', 160};
            app.MainGrid.Padding = [0 0 0 0];
            app.MainGrid.ColumnSpacing = 0;
            app.MainGrid.RowSpacing = 0;

            % =============================================================
            % 1. SIDEBAR
            % =============================================================
            app.SidePanel = uipanel(app.MainGrid);
            app.SidePanel.Layout.Row = [1 2];
            app.SidePanel.Layout.Column = 1;
            app.SidePanel.BackgroundColor = app.ColorBgPanel;
            app.SidePanel.BorderType = 'none';

            app.SidebarGrid = uigridlayout(app.SidePanel, [5, 1]);
            app.SidebarGrid.RowHeight = {130, 60, '1x', 50, 60};
            app.SidebarGrid.Padding = [20 20 20 20];
            app.SidebarGrid.RowSpacing = 15;

            % --- A. Branding ---
            app.HeaderPanel = uipanel(app.SidebarGrid);
            app.HeaderPanel.BackgroundColor = app.ColorBgPanel;
            app.HeaderPanel.BorderType = 'none';

            tGrid = uigridlayout(app.HeaderPanel, [2,1]);
            tGrid.Padding = [0 0 0 0];
            tGrid.RowSpacing = 0;
            tGrid.RowHeight = {'1x', 25};

            app.AppTitle = uilabel(tGrid);
            app.AppTitle.Layout.Row = 1;
            app.AppTitle.Interpreter = 'html';
            app.AppTitle.FontSize = 36;
            app.AppTitle.FontColor = [0.24 0.24 0.24];
            app.AppTitle.VerticalAlignment = 'bottom';

            app.AppTitle.Text = [...
                '<b style="color:rgb(0,122,255)">F</b>lexible ' ...
                '<b style="color:rgb(0,122,255)">A</b>daptive<br>' ...
                '<b style="color:rgb(0,122,255)">S</b>table ' ...
                '<b style="color:rgb(0,122,255)">C</b>lustering'];

            app.AppSubtitle = uilabel(tGrid);
            app.AppSubtitle.Layout.Row = 2;
            app.AppSubtitle.Text = "Shao SHI @SUSTech";
            app.AppSubtitle.FontColor = app.ColorTextDim;
            app.AppSubtitle.FontSize = 11;
            app.AppSubtitle.FontWeight = 'bold';
            app.AppSubtitle.VerticalAlignment = 'top';

            % --- B. Import & Save ---
            app.ImportGrid = uigridlayout(app.SidebarGrid, [1, 2]);
            app.ImportGrid.Padding = [0 0 0 0];
            app.ImportGrid.ColumnSpacing = 10;

            app.BtnLoad = uibutton(app.ImportGrid, 'Text', ' LOAD DATA', ...
                'Icon', '', 'FontWeight', 'bold');
            app.BtnLoad.BackgroundColor = [0.25 0.27 0.32];
            app.BtnLoad.FontColor = [1 1 1];
            app.BtnLoad.FontSize = 12;
            app.BtnLoad.ButtonPushedFcn = @(s,e) app.loadData();

            app.BtnSave = uibutton(app.ImportGrid, 'Text', ' SAVE RESULTS', ...
                'Icon', '', 'FontWeight', 'bold');
            app.BtnSave.BackgroundColor = [0.25 0.27 0.32]; 
            app.BtnSave.FontColor = [1 1 1];
            app.BtnSave.FontSize = 12;
            app.BtnSave.Enable = 'off';
            app.BtnSave.ButtonPushedFcn = @(s,e) app.saveData();

            % --- C. Parameters ---
            app.ScrollContainer = uipanel(app.SidebarGrid);
            app.ScrollContainer.BackgroundColor = app.ColorBgPanel;
            app.ScrollContainer.BorderType = 'none';

            % 16 rows to accommodate spacing
            app.ScrollLayout = uigridlayout(app.ScrollContainer, [16, 2]);
            
            % [UPDATED] Reduced row height to 24px for compact look
            app.ScrollLayout.RowHeight = repmat({24}, 1, 16);
            app.ScrollLayout.ColumnWidth = {'1x', '1.1x'};
            app.ScrollLayout.Scrollable = 'on';
            app.ScrollLayout.Padding = [0 0 10 0];
            
            % [UPDATED] Reduced row spacing to 4px
            app.ScrollLayout.RowSpacing = 4;

            r=1;

            function addSeparator(titleText)
                lbl = uilabel(app.ScrollLayout, 'Text', titleText);
                lbl.Layout.Row = r; lbl.Layout.Column = [1 2];
                lbl.FontColor = app.ColorAccent;
                lbl.FontWeight = 'bold';
                lbl.FontSize = 10;
                lbl.VerticalAlignment = 'bottom';
                r=r+1;
            end

            function [lbl, comp] = addParam(text, widget, tooltip)
                lbl = uilabel(app.ScrollLayout, 'Text', text);
                lbl.Layout.Row = r; lbl.Layout.Column = 1;
                lbl.FontColor = app.ColorLabel;
                lbl.FontSize = 11;
                lbl.FontWeight = 'bold';

                widget.Layout.Row = r; widget.Layout.Column = 2;
                if nargin > 2, widget.Tooltip = tooltip; end
                app.styleInput(widget); 
                comp = widget;
                r=r+1;
            end

            % --- PARAMETER LAYOUT START ---
            
            addSeparator('CORE ALGORITHM');

            % 1. Metric Selection
            app.DropAlgorithm = uidropdown(app.ScrollLayout, ...
                'Items', {'dual-cosine', 'cosine', 'euclidean', 'l1-norm', 'minimum', ...
                'maximum', 'algebraic', 'logarithmic', 'geometric', ...
                'harmonic', 'enhanced harmonic', 'entropy', 'weighted entropy'}, ...
                'Value', 'dual-cosine');
            addParam('Similarity Algorithm', app.DropAlgorithm, 'Distance Metric');
            % Trigger layout update when changed
            app.DropAlgorithm.ValueChangedFcn = @(s,e) app.updateLayout();

            % 2. Idx Inputs (Moved Here)
            app.RowIdxPos = r;
            app.EditIdxPos = uieditfield(app.ScrollLayout, 'text', 'Value', '1:300');
            [app.LblIdxPos, ~] = addParam('Idx Pos', app.EditIdxPos);

            app.RowIdxNeg = r;
            app.EditIdxNeg = uieditfield(app.ScrollLayout, 'text', 'Value', '301:600');
            [app.LblIdxNeg, ~] = addParam('Idx Neg', app.EditIdxNeg);

            % 3. Optimizer
            app.DropStrategy = uidropdown(app.ScrollLayout, 'Items', {'DASS','SF'}, 'Value', 'DASS');
            addParam('Optimizer', app.DropStrategy, 'Optimization Strategy');

            addSeparator('THRESHOLDS');

            app.EditSimInner = uieditfield(app.ScrollLayout, 'numeric', 'Value', 0.7);
            addParam('Sim (intra)', app.EditSimInner);

            app.EditSimInter = uieditfield(app.ScrollLayout, 'numeric', 'Value', 0.7);
            addParam('Sim (inter)', app.EditSimInter);

            app.EditMinVol = uieditfield(app.ScrollLayout, 'numeric', 'Value', 2);
            addParam('Min Support', app.EditMinVol);

            addSeparator('CONSTRAINTS');

            app.EditInitLim = uieditfield(app.ScrollLayout, 'numeric', 'Value', 8);
            addParam('Seed Numbers', app.EditInitLim);

            app.EditMaxClust = uieditfield(app.ScrollLayout, 'numeric', 'Value', 50);
            addParam('Max Clusters', app.EditMaxClust);

            app.EditMaxIter = uieditfield(app.ScrollLayout, 'numeric', 'Value', 200);
            addParam('Max Iterations', app.EditMaxIter);
            
            % --- PARAMETER LAYOUT END ---

            % --- D. Status Area ---
            app.StatusPanel = uipanel(app.SidebarGrid);
            app.StatusPanel.BackgroundColor = [0.85 0.85 0.88]; 
            app.StatusPanel.BorderType = 'none';
            sGrid = uigridlayout(app.StatusPanel, [1,2]);
            sGrid.Padding=[10 0 10 0]; sGrid.ColumnWidth={'fit','1x'};

            app.StatusLamp = uilamp(sGrid);
            app.StatusLamp.Color = app.ColorSuccess;

            app.StatusLabel = uilabel(sGrid);
            app.StatusLabel.Text = "SYSTEM READY";
            app.StatusLabel.FontColor = [0.3 0.3 0.3]; 
            app.StatusLabel.FontWeight = 'bold';
            app.StatusLabel.FontSize = 10;

            % --- E. Run Button ---
            app.BtnRun = uibutton(app.SidebarGrid, 'Text', 'RUN ANALYSIS', ...
                'FontSize', 14, 'FontWeight', 'bold');
            app.BtnRun.BackgroundColor = app.ColorAccent;
            app.BtnRun.FontColor = [1 1 1];
            app.BtnRun.ButtonPushedFcn = @(s,e) app.runAnalysis();

            % =============================================================
            % 2. RIGHT CONTENT AREA
            % =============================================================

            app.ContentPanel = uipanel(app.MainGrid);
            app.ContentPanel.Layout.Row = 1; app.ContentPanel.Layout.Column = 2;
            app.ContentPanel.BackgroundColor = app.ColorContentBg;
            app.ContentPanel.BorderType = 'none';

            contentWrapper = uigridlayout(app.ContentPanel, [1,1]);
            contentWrapper.Padding = [20 20 20 10];

            app.TabGroup = uitabgroup(contentWrapper);

            function t = createTab(title)
                t = uitab(app.TabGroup, 'Title', title);
                t.BackgroundColor = app.ColorContentBg;
            end

            % -- Tab 1: Convergence --
            app.TabConvergence = createTab('Convergence');
            g1 = uigridlayout(app.TabConvergence, [2,1]);
            g1.RowSpacing = 20; g1.Padding = [10 10 10 10];

            app.PanelConv = uipanel(g1, 'Title', 'ITERATION METRICS', ...
                'BackgroundColor','w', 'FontWeight', 'bold', 'FontSize', 11, ...
                'ForegroundColor', [0.4 0.4 0.4]);
            app.PanelConv.BorderType = 'line';
            app.PanelConv.HighlightColor = [0.9 0.9 0.9];

            app.PanelOutliers = uipanel(g1, 'Title', 'OUTLIER STATISTICS', ...
                'BackgroundColor','w', 'FontWeight', 'bold', 'FontSize', 11, ...
                'ForegroundColor', [0.4 0.4 0.4]);
            app.PanelOutliers.BorderType = 'line';
            app.PanelOutliers.HighlightColor = [0.9 0.9 0.9];

            % -- Tab 2: Heatmap --
            app.TabHeatmap = createTab('Similarity Matrix');
            g2 = uigridlayout(app.TabHeatmap, [1,1]); g2.Padding = [5 5 5 5];
            app.PanelHeatmap = uipanel(g2, 'BackgroundColor','w', 'BorderType','none');

            % -- Tab 3: Distribution --
            app.TabDistribution = createTab('Cluster Distribution');
            g3 = uigridlayout(app.TabDistribution, [1,1]); g3.Padding = [5 5 5 5];
            app.PanelDist = uipanel(g3, 'BackgroundColor','w', 'BorderType','none');

            % =============================================================
            % 3. LOG AREA
            % =============================================================
            app.LogPanel = uipanel(app.MainGrid);
            app.LogPanel.Layout.Row = 2; app.LogPanel.Layout.Column = 2;
            app.LogPanel.BackgroundColor = [0.15 0.17 0.21]; 
            app.LogPanel.BorderType = 'none';

            logWrapper = uigridlayout(app.LogPanel, [2,1]);
            logWrapper.RowHeight = {25, '1x'};
            logWrapper.Padding = [10 5 10 5];
            logWrapper.RowSpacing = 0;

            lblLog = uilabel(logWrapper, 'Text', ' EXECUTION LOG');
            lblLog.FontName = 'Consolas';
            lblLog.FontSize = 10;
            lblLog.FontColor = [0.6 0.6 0.6];

            app.LogText = uitextarea(logWrapper);
            app.LogText.Editable = 'off';
            app.LogText.FontName = 'Consolas';
            app.LogText.FontSize = 11;
            app.LogText.BackgroundColor = [0.15 0.17 0.21];
            app.LogText.FontColor = [0.7 0.9 0.7];
            app.LogText.Value = "System Ready. Please load a dataset to begin.";

            % Apply initial visibility state
            app.updateLayout();

            app.UIFigure.Visible = 'on';
        end
    end
end