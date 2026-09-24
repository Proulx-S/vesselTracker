function [cmap, alphaProfile] = resolveColormapFcn(method, name, N)
% RESOLVECOLORMAPFCN  Resolves a 'method'/'name' colormap spec into an actual [N x 3] colormap matrix.
% Promoted (2026-08-24) out of resolveStatMapColor.m's own resolveColormapFcnLocal once a SECOND caller
% (plotGaussianFitPanels.m's own opts.residualColormap, xc's residual colormap) needed the identical
% 'fname'/'fname(args...)' spec grammar, this codebase's established "promote once 2+ callers need it"
% convention (see resolveStatMapColor.m's own file header for the FIRST instance of that convention).
%
% INPUT
%   method : '' (plain MATLAB/on-path builtin -- resolved via name, see below) | 'fname' | 'fname(args...)'
%            (2026-08-22, Seb's own explicit ask: "a way to pass all arguments to colormap_divergingHue
%            from opts.colormap.UL/OL") -- e.g. 'colormap_divergingHue', or
%            'colormap_divergingHue({[200 250],[10 40]}, 0.3, 0.1, [], [])' to pass that function's own
%            hues/wNeutral/wTransition positional args through. SAME 'fname(args...)' spec convention
%            resolveStatMapSignificance.m's own opts.fdr already uses.
%
%            ONLY 'colormap_divergingHue' is special-cased below (that function's own positional
%            signature needs hues/wNeutral/wTransition/lNeutral/cOuter BEFORE N, see its own header --
%            calling it with a bare N would misread N as its own first positional arg, hues). When args
%            ARE given, they are passed through VERBATIM as colormap_divergingHue's own FIRST positional
%            args (hues, wNeutral, wTransition, lNeutral, cOuter, in that order), with N ALWAYS appended
%            by this function as a PLAIN trailing token -- NOT position-aware, so there is no
%            padding/skipping on this end. To reach a LATER arg (e.g. cOuter, position 5), ALL earlier
%            ones must be supplied too (use [] for any of those left at that function's own default,
%            matching its own "pass [] to keep a default" convention exactly) -- e.g.
%            'colormap_divergingHue([], [], [], [], 0.4)' to set ONLY cOuter, NOT the wrong (but easy to
%            mistakenly write) 'colormap_divergingHue(0.4)' (sets hues=0.4, not cOuter). N is NOT exposed
%            for override here, since it also drives this function's own callers' color-index
%            quantization (a caller-chosen N would need threading back into that step too, not just the
%            colormap call) -- nor is plotProfiles (7th positional arg, a side-effecting PNG save), which
%            has no place in a quiet color-resolution call. Any OTHER method name (not
%            colormap_divergingHue) falls back to a plain fcn(N) call when no args are given, or
%            fcn(args..., N) when they are -- a reasonable default for a future custom method that DOES
%            accept a bare N as its own last positional arg, extended with its own special case here
%            (like colormap_divergingHue's own) if that turns out not to hold.
%   name   : builtin/on-path colormap function name (e.g. 'gray', 'parula') -- used when method=='',
%            ignored otherwise.
%   N      : colormap length.
%
% OUTPUT
%   cmap        : [N x 3] resolved colormap.
%   alphaProfile: [N x 1] per-color alpha/opacity carried through from colormap_divergingHue's own
%                 info.alpha (2026-08-21, e.g. lNeutral='trans') when method resolves to that function --
%                 all-ones otherwise (every builtin, and colormap_divergingHue itself unless
%                 lNeutral='trans'). Only colormap_divergingHue is called with 2 outputs (cmap, info) to
%                 recover this -- any other method's alphaProfile is always ones(N,1), since only that
%                 function's own info struct is known to carry an alpha field.
%
%   cmap = resolveColormapFcn('', 'gray', 256);
%   cmap = resolveColormapFcn('colormap_divergingHue', '', 256);
%   cmap = resolveColormapFcn('colormap_divergingHue({[200 250],[10 40]}, 0.3, 0.1, [], [])', '', 256);
%
% See also resolveStatMapColor, resolveStatMapSignificance.

    if isempty(method)
        assert(ischar(name) && ~isempty(name), 'resolveColormapFcn:nameRequired', ...
            'name is required (e.g. ''gray'', ''parula'') when method is empty.');
        cmapFcn = str2func(name);
        cmap = cmapFcn(N);
        alphaProfile = ones(N,1);   % builtins carry no alpha concept
        return
    end
    % 2026-08-22 BUG FIX: MATLAB's regexp('tokens') silently DROPS a capturing group nested inside
    % ANOTHER capturing group -- caught HERE first (a direct test comparing against a hand-built
    % reference colormap, not just a boolean outcome that could coincidentally match either way), then
    % confirmed to ALSO affect resolveStatMapSignificance.m's own identical former pattern -- see that
    % file's own applyFdrLocal for the fuller incident writeup. The former pattern here,
    % '(\((.*)\))?' (an outer optional group wrapping an inner one), silently returned NO 3rd token at
    % all, so argsStr was ALWAYS '' and 'fname(args...)' methods (e.g. a custom hue/wNeutral spec) fell
    % back to a bare fname(N) call, args completely ignored. Fixed by capturing the WHOLE '(...)' suffix
    % (parens included) as ONE group and stripping the parens manually via indexing.
    tok = regexp(strtrim(method), '^([A-Za-z]\w*)\s*(\(.*\))?$', 'tokens', 'once');
    assert(~isempty(tok), 'resolveColormapFcn:badMethodSpec', ...
        'method ''%s'' is not a valid ''fname'' or ''fname(args...)'' spec.', method);
    fname = tok{1};
    argsStr = ''; if numel(tok)>=2 && ~isempty(tok{2}); argsStr = tok{2}(2:end-1); end
    if strcmp(fname, 'colormap_divergingHue')
        % colormap_divergingHue's own info.alpha (2026-08-21, all-ones unless lNeutral='trans') is
        % threaded straight through as alphaProfile -- this is the ONLY method special-cased here, so
        % it's also the only one that can actually produce a non-trivial alpha.
        if isempty(argsStr)
            [cmap, info] = colormap_divergingHue([], [], [], [], [], N);
        else
            [cmap, info] = eval(sprintf('colormap_divergingHue(%s, N)', argsStr));
        end
        if isfield(info, 'alpha'); alphaProfile = info.alpha; else; alphaProfile = ones(N,1); end
    elseif isempty(argsStr)
        cmapFcn = str2func(fname);
        cmap = cmapFcn(N);
        alphaProfile = ones(N,1);
    else
        cmap = eval(sprintf('%s(%s, N)', fname, argsStr));
        alphaProfile = ones(N,1);
    end
end
