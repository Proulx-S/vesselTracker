function h = drawMeasuredVsPredictedScatter(ax, predicted, measured, color, unityLineColor, unityLineStyle, showGrid)
% DRAWMEASUREDVSPREDICTEDSCATTER  Draws a scatter of every fit voxel's predicted (x) vs. measured (y)
% intensity plus a unity reference line, into a GIVEN axes -- same "caller owns the axes" convention as
% plotVesselPatch.m. Sets DataAspectRatio/PlotBoxAspectRatio to [1 1 1] itself (same intensity units both
% axes, so the unity line reads as a true 45 degrees) -- unlike the radial-profile panels, this 1:1 data
% scale is part of what makes this panel meaningful, not just a caller-side styling choice, so it lives
% here rather than being left to the caller. Promoted out of fitVesselDiag.m's own panel-4 inline code
% (2026-08-18).
%
% INPUT
%   ax             : target axes.
%   predicted      : [n x 1] each fit voxel's own predicted intensity.
%   measured       : [n x 1] each fit voxel's own measured intensity (SAME order as predicted).
%   color          : RGB triplet for the scatter markers.
%   unityLineColor : RGB triplet (or []) for the reference line -- default [0.6 0.6 0.6] (dim gray, this
%                    file's own original behavior -- e.g. fitVesselDiag.m's own panel 4 leaves this at
%                    default so its neutral reference line stays independent of whatever method color
%                    the markers use).
%   unityLineStyle : LineStyle char for the reference line -- default '--' (this file's own original
%                    behavior).
%   showGrid       : true/false -- grid(ax,'on') if true -- default false (this file's own original
%                    behavior).
%
% OUTPUT
%   h : the scatter handle.
%
%   drawMeasuredVsPredictedScatter(ax, predicted, measured, [1 0.85 0.1]);                    % original behavior
%   drawMeasuredVsPredictedScatter(ax, predicted, measured, [1 1 1], [1 1 1], '-', true);      % solid white unity line + grid

    if nargin<7 || isempty(showGrid);       showGrid       = false;         end
    if nargin<6 || isempty(unityLineStyle); unityLineStyle = '--';          end
    if nargin<5 || isempty(unityLineColor); unityLineColor = [0.6 0.6 0.6]; end

    hold(ax,'on'); set(ax,'Color','k','XColor','w','YColor','w');
    h = scatter(ax, predicted, measured, 40, color, 'filled', 'MarkerEdgeColor','k');
    lims = [min([predicted(:);measured(:)]) max([predicted(:);measured(:)])];
    if diff(lims)==0; lims = lims + [-1 1]; end
    plot(ax, lims, lims, unityLineStyle, 'Color',unityLineColor);
    xlim(ax,lims); ylim(ax,lims);
    set(ax, 'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1]);
    if showGrid; grid(ax,'on'); end
end
