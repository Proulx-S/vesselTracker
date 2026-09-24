function drawScaleBar(ax, lengthMm, color, location)
% DRAWSCALEBAR  Draws a physical-length reference bar + text label into a GIVEN axes already in
% center-voxel-relative mm units (see plotVesselPatch.m's own toX/toY convention) -- same "caller owns
% the axes" convention as plotVesselPatch.m. Reads the bar's own placement off the axes' CURRENT
% XLim/YLim (via a fixed margin fraction), so it stays correctly positioned regardless of the axes'
% own XLim/YLim/CLim/colormap -- draw it AFTER those are finalized (a later XLim/YLim change would
% require redrawing the bar).
%
% INPUT
%   ax       : target axes, already in mm units (e.g. via plotVesselPatch.m/drawModelImage.m).
%   lengthMm : bar length, mm (e.g. 1).
%   color    : RGB triplet for the bar/text (default [1 1 1], white).
%   location : 'southeast' (default) | 'southwest' | 'northeast' | 'northwest' -- which corner to anchor
%              to, margined in from the axes' own current XLim/YLim.
%
%   drawScaleBar(ax, 1, [1 1 1]);
%   drawScaleBar(ax, 1, [1 1 1], 'northwest');

    if nargin<4 || isempty(location); location = 'southeast'; end
    if nargin<3 || isempty(color);    color    = [1 1 1];     end

    xl = xlim(ax); yl = ylim(ax);
    dx = diff(xl); dy = diff(yl);
    marginX = 0.08*dx; marginY = 0.08*dy;

    switch location
        case 'southeast'
            x1 = xl(2)-marginX; x0 = x1-lengthMm; y0 = yl(1)+marginY; textVA = 'bottom';
        case 'southwest'
            x0 = xl(1)+marginX; x1 = x0+lengthMm; y0 = yl(1)+marginY; textVA = 'bottom';
        case 'northeast'
            x1 = xl(2)-marginX; x0 = x1-lengthMm; y0 = yl(2)-marginY; textVA = 'top';
        case 'northwest'
            x0 = xl(1)+marginX; x1 = x0+lengthMm; y0 = yl(2)-marginY; textVA = 'top';
        otherwise
            error('drawScaleBar:badLocation', ...
                'location must be ''southeast''|''southwest''|''northeast''|''northwest'' -- got ''%s''.', location);
    end

    hold(ax,'on');
    plot(ax, [x0 x1], [y0 y0], '-', 'Color',color, 'LineWidth',2);
    textOffsetY = 0.04*dy; if strcmp(textVA,'top'); textOffsetY = -textOffsetY; end
    text(ax, (x0+x1)/2, y0+textOffsetY, sprintf('%g mm', lengthMm), ...
        'Color',color, 'HorizontalAlignment','center', 'VerticalAlignment',textVA, 'FontSize',7);
end
