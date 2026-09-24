function center = drawVesselCenter(tsIm, opts)
% DRAWVESSELCENTER  Manually pick the center pixel of one or more vessels by clicking on the
% time-averaged image of a full-image ts-track (loadNiftiTs.m's output). Returns the clicked pixel
% positions, one row per vessel, in the SAME pixel convention every crop in this pipeline uses:
%
%   center : [nVessel x 2] double, 1-based (x, y) = (column, row) of tsIm.im{r}.
%
% CONTROLS (figure window):
%   left click            add a vessel at the clicked pixel (rounded to the nearest pixel center)
%   backspace / delete    remove the last vessel added
%   Enter                 accept -- returns the list (an EMPTY list is a valid answer, see CACHE)
%   Escape, or closing    cancel -- errors (drawVesselCenter:cancelled); nothing is cached
% The figure toolbar's zoom/pan tools work as usual: while a zoom/pan tool is active, clicks go to that
% tool, so toggle it off again to place a vessel. opts.zoomXlim/.zoomYlim set the initial view.
%
% CACHE -- this is manual work, so the accepted result is written to disk and reused (the convention
% .bass/pkm/patterns/manual-work-caching.md describes): ONE file per opts.id under opts.cacheDir,
% <cacheDir>/<id>.mat holding `center`. It records the human's RULING, so an accepted EMPTY list is
% cached too (and is then a cache HIT, not "never asked"). opts.force is the standard three-level ladder:
%   0  (default)  cache hit -> return it WITHOUT opening a figure; miss -> prompt.
%   1             prompt ALWAYS, seeded from the cached list (edit an earlier answer).
%   2             prompt ALWAYS, starting blank (the cache file is left untouched until you accept).
% Cancel never writes. force>2 errors. The cache directory should be git-tracked by the calling
% project (see that same pattern page): a click is not regenerable by re-running an algorithm.
%
% INPUT
%   tsIm : loadNiftiTs.m output (.im {1 x nRun} of [ny x nx x 1 x T]); the display is the mean over
%          every frame of every run.
%   opts : struct --
%       .id        char, the cache key (REQUIRED, no default -- name the dataset/subject/session so two
%                  images can never share an entry)                                          ('')
%       .cacheDir  char                                    (<pwd>/drawVesselCenterCache)
%       .force     0 | 1 | 2 -- see CACHE                                                     (0)
%       .clim      [] = [1 99.9] percentiles of the displayed image | [lo hi]                 ([])
%       .zoomXlim/.zoomYlim  [] = full image | [lo hi] pixel range for the initial view       ([])
%       .markerSize  positive scalar                                                          (12)
%
% NO-ARG CALL -- drawVesselCenter() prints this opts listing and returns opts populated with defaults.
%
%   opts = drawVesselCenter; opts.id = 'sub01_ses1';
%   center = drawVesselCenter(tsIm, opts);         % click, Enter
%   vessel = makeVessel(tsIm, center);
%
% See also loadNiftiTs, makeVessel.

    if nargin==0
        printOptsHelp();
        center = defaultOpts();
        return
    end
    if nargin<2 || isempty(opts); opts = struct(); end
    opts = fillMissing(opts, defaultOpts());
    assert(~isempty(opts.id) && (ischar(opts.id) || isstring(opts.id)), 'drawVesselCenter:noId', ...
        'opts.id is required (the cache key) -- name the dataset this image belongs to.');
    opts.id = char(opts.id);
    force = double(opts.force);
    assert(isscalar(force) && ismember(force, [0 1 2]), 'drawVesselCenter:badForce', ...
        'opts.force must be 0, 1 or 2 (got %s).', mat2str(opts.force));
    assert(isstruct(tsIm) && isfield(tsIm,'im') && ~isempty(tsIm.im), 'drawVesselCenter:noIm', ...
        'tsIm must be a loadNiftiTs.m output with a non-empty .im.');

    cacheFile = fullfile(opts.cacheDir, [opts.id '.mat']);
    cached = [];
    if force < 2 && exist(cacheFile, 'file')==2
        S = load(cacheFile, 'center');
        cached = S.center;
        if force == 0
            fprintf('drawVesselCenter: %s -- %d vessel(s) from cache <- %s\n', opts.id, size(cached,1), cacheFile);
            center = cached;
            return
        end
        fprintf('drawVesselCenter: %s -- prompting, seeded from cache (%d vessel(s)) <- %s\n', opts.id, size(cached,1), cacheFile);
    elseif force == 2
        fprintf('drawVesselCenter: %s -- prompting fresh (opts.force=2, cache not loaded)\n', opts.id);
    else
        fprintf('drawVesselCenter: %s -- no cache entry, prompting\n', opts.id);
    end
    seed = zeros(0, 2);
    if ~isempty(cached); seed = cached; end

    % Mean over every frame of every run -- accumulated per run so no run is ever duplicated in memory.
    im = tsImTimeMean(tsIm);

    % ---- interactive part: nested callbacks share `pts`, `hMark`, `accepted` -----------------------
    pts = seed;
    accepted = false;
    fig = figure('Name', ['drawVesselCenter -- ' opts.id], 'NumberTitle', 'off', 'Color', 'k');
    ax  = axes('Parent', fig);
    hIm = imagesc(ax, im);
    axis(ax, 'image'); colormap(ax, 'gray'); hold(ax, 'on');
    cl = opts.clim;
    if isempty(cl); cl = prctile(im(:), [1 99.9]); end
    clim(ax, cl);
    if ~isempty(opts.zoomXlim); xlim(ax, opts.zoomXlim); end
    if ~isempty(opts.zoomYlim); ylim(ax, opts.zoomYlim); end
    ax.XColor = 'w'; ax.YColor = 'w';
    title(ax, {opts.id, 'click = add vessel | backspace = remove last | Enter = accept | Esc = cancel'}, 'Color', 'w');
    hMark = gobjects(0);
    redraw();

    hIm.ButtonDownFcn  = @onClick;
    fig.KeyPressFcn    = @onKey;
    fig.CloseRequestFcn = @onClose;
    uiwait(fig);

    if ~accepted
        error('drawVesselCenter:cancelled', 'drawVesselCenter: %s -- cancelled, nothing cached.', opts.id);
    end
    center = pts;
    if ~exist(opts.cacheDir, 'dir'); mkdir(opts.cacheDir); end
    save(cacheFile, 'center');
    fprintf('drawVesselCenter: %s -- %d vessel(s) accepted, cached -> %s\n', opts.id, size(center,1), cacheFile);

    % ---- nested callbacks ---------------------------------------------------------------------
    function onClick(~, ~)
        if ~strcmp(fig.SelectionType, 'normal'); return; end   % left click only
        p = ax.CurrentPoint;
        x = round(p(1,1)); y = round(p(1,2));
        if x < 1 || y < 1 || x > size(im,2) || y > size(im,1); return; end
        pts(end+1, :) = [x y];
        redraw();
    end
    function onKey(~, evt)
        switch evt.Key
            case {'backspace', 'delete'}
                if ~isempty(pts); pts(end, :) = []; redraw(); end
            case 'return'
                accepted = true;
                delete(fig);
            case 'escape'
                delete(fig);
        end
    end
    function onClose(~, ~)
        delete(fig);          % accepted stays false -> cancel
    end
    function redraw()
        delete(hMark(isgraphics(hMark)));
        hMark = gobjects(0);
        for k = 1:size(pts, 1)
            hMark(end+1) = plot(ax, pts(k,1), pts(k,2), 'o', 'Color', [1 0.8 0], ...
                'MarkerSize', opts.markerSize, 'LineWidth', 1.5); %#ok<AGROW>
            hMark(end+1) = text(ax, pts(k,1) + 2, pts(k,2) - 2, sprintf('%d', k), ...
                'Color', [1 0.8 0], 'FontWeight', 'bold'); %#ok<AGROW>
        end
        drawnow;
    end
end

% ---------------------------------------------------------------------------
function im = tsImTimeMean(tsIm)
    acc = 0; n = 0;
    for r = 1:numel(tsIm.im)
        acc = acc + sum(double(tsIm.im{r}(:,:,1,:)), 4);
        n   = n + size(tsIm.im{r}, 4);
    end
    im = acc / n;
end

function opts = defaultOpts()
    opts.id         = '';
    opts.cacheDir   = fullfile(pwd, 'drawVesselCenterCache');
    opts.force      = 0;
    opts.clim       = [];
    opts.zoomXlim   = [];
    opts.zoomYlim   = [];
    opts.markerSize = 12;
end

function opts = fillMissing(opts, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}); opts.(fn{i}) = d.(fn{i}); end
    end
end

function printOptsHelp()
    fprintf('drawVesselCenter opts -- allowed value(s) per field (default in parentheses):\n');
    fprintf('  opts.id         : char, cache key -- REQUIRED, no default                         ('''')\n');
    fprintf('  opts.cacheDir   : char                                          (<pwd>/drawVesselCenterCache)\n');
    fprintf('  opts.force      : 0 = reuse cache | 1 = re-prompt seeded from cache | 2 = re-prompt blank   (0)\n');
    fprintf('  opts.clim       : [] = [1 99.9] percentiles | [lo hi]                            ([])\n');
    fprintf('  opts.zoomXlim   : [] = full image | [lo hi] pixel range of the initial view       ([])\n');
    fprintf('  opts.zoomYlim   : [] = full image | [lo hi] pixel range of the initial view       ([])\n');
    fprintf('  opts.markerSize : positive scalar                                                (12)\n\n');
end
