%Puget Sound routines
%04/02 MHA
%
%Routines for plotting contours of Puget Sound bathymetry, as well as designing cruise tracks
%and plotting along tracks.
%Functions:
%ComputeSWIMSSpacings: return the estimated distance and time between successive SWIMS profiles given tow speed and cycle depth
%
%LatLonColorPlotFcn: plot a variable along a track, on top of Puget Sound bathymetry contours
%MakeGridPattern: Make waypoints for a rectangular grid of specified dimensions, number of waypoints, location, and angle.
%MakeStarPattern: Make waypoints for a star pattern of specified radius, and position of starting location.
%OverplotProfileLocs: given a list of waypoint lat/lon pairs, overplot the spacing of SWIMS profiles using ComputeSWIMSSpacings (above).
%PSBathy: return a matrix of Puget Sound bathymetry
%PSBathyPlot: Specify plotting limits, the routine plots grayscale bathymetry
%PSShorePlot: Specify plotting limits, the routine plots shorelines in black
%SWIMSprofileLocs: Compute (using ComputeSWIMSSpacings.m) the locations of profiles given a list of waypoint lat/lon pairs

%Directory: PS_tides 
%Dave's routines for Lavelle's model
%Plus I have added a function that returns the tide for a given segment at the desired UTC ydays.