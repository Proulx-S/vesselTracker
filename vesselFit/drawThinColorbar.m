function cb = drawThinColorbar(ax, cmap, cLim, width, side, gap)
% DRAWTHINCOLORBAR  A thin, vertical colorbar flush against one side of ax, spanning its own full
% height -- for an axes whose own colormap is ALREADY spoken for by something else (e.g. a grayscale
% underlay, or a truecolor RGB overlay), so a bare `colorbar(ax)` would show the WRONG scale. Promoted
% (2026-08-22) out of drawStatMapOverlay.m once a SECOND caller (plotVessels.m's own vesselPatchFrame
% route, giving the BACKGROUND role its own independently-configurable colorbar -- opts.colormap.UL --
% not just the overlay's opts.colormap.OL) needed the identical mechanism, this codebase's established
% "promote once 2+ callers need it" convention.
%
% MECHANISM: a small, dedicated, invisible-except-for-its-colorbar axes stacked exactly on ax's own
% position (Units forced to 'normalized' on both so the two positions are comparable), given cmap/cLim
% itself. Does NOT need linkaxes/addlistener (unlike a second axes carrying actual image DATA) -- a
% colorbar draws no spatial data of its own, so there is nothing to keep in pixel-for-pixel sync across
% pan/zoom.
%
% INPUT
%   ax    : the axes this colorbar describes the color scale FOR (not the axes it's drawn IN).
%   cmap  : [N x 3] colormap.
%   cLim  : [lo hi] the color range cmap spans.
%   width : ACROSS-the-bar thickness, in the SAME normalized figure units as ax.Position -- i.e. the
%           bar's width for 'right'/'left' and its height for 'bottom'/'top'.
%   side  : 'right' | 'left' | 'bottom' | 'top' -- which side of ax to sit flush against. The two
%           HORIZONTAL sides set the colorbar's own Location as well as its Position (a colorbar
%           defaults to a vertical Location, and a wide-and-short Position on a vertical bar yields a
%           squashed vertical bar, not a horizontal one).
%   gap   : OPTIONAL clearance between the axes and the bar, normalized figure units. Default 0.01,
%           which is right for 'right'/'left'; a 'bottom'/'top' bar on an axes that keeps its own tick
%           labels needs more (they are drawn in exactly that space).
%
% OUTPUT
%   cb : the colorbar object.
%
%   cb = drawThinColorbar(ax, cmap, [lo hi], 0.02, 'right');
%
% See also resolveStatMapColor, drawStatMapOverlay.

    fig = ancestor(ax,'figure');
    origUnits = ax.Units; ax.Units = 'normalized'; axPos = ax.Position; ax.Units = origUnits;

    axCB = axes('Parent',fig, 'Units','normalized', 'Position',axPos, ...
        'Visible','off', 'HitTest','off', 'Colormap',cmap, 'CLim',cLim);
    if nargin<6 || isempty(gap); gap = 0.01; end   % normalized-units clearance from the axes
    % gap is an ARGUMENT (2026-09-03) rather than the fixed 0.01 it was, because a 'bottom' bar has to
    % clear whatever the axes draws below itself -- x tick labels and an xlabel -- whereas a 'right'
    % bar only ever clears the axes box. 0.01 put the bar straight through the tick labels of a panel
    % that keeps its x ruler (verified visually).
    % 'bottom'/'top' added 2026-09-03 (fitVesselTimeSeriesDiag.m's own correlation panel wanted its
    % scale underneath rather than beside it). Purely additive -- 'right'/'left' behave exactly as
    % before -- but note the HORIZONTAL cases need Location set on the colorbar itself, not just a
    % Position: a colorbar defaults to a vertical Location, and assigning a wide-and-short Position to
    % a vertical bar gives a squashed vertical bar rather than a horizontal one.
    switch side
        case 'right'
            cb = colorbar(axCB);
            cb.Units = 'normalized';
            cb.Position = [axPos(1)+axPos(3)+gap, axPos(2), width, axPos(4)];
        case 'left'
            cb = colorbar(axCB);
            cb.Units = 'normalized';
            cb.Position = [axPos(1)-gap-width, axPos(2), width, axPos(4)];
        case 'bottom'
            cb = colorbar(axCB, 'Location','southoutside');
            cb.Units = 'normalized';
            cb.Position = [axPos(1), axPos(2)-gap-width, axPos(3), width];
        case 'top'
            cb = colorbar(axCB, 'Location','northoutside');
            cb.Units = 'normalized';
            cb.Position = [axPos(1), axPos(2)+axPos(4)+gap, axPos(3), width];
        otherwise
            error('drawThinColorbar:badSide', ...
                'side must be ''right''|''left''|''bottom''|''top'' (got ''%s'').', side);
    end
end
