function fcs_corr_add(h, trace, channels)
% function fcs_corr_add(h, trace, channels)
% adds a piece of trace to the correlator. A correlator initialized with
% fcs_corr_init has to be passed to the routine. This routine can add data
% to i) autocorrelate a time trace, ii) crosscorrelate two time traces or
% iii) crosscorrelate a two-channel FIFO trace. If the entry 'channels' is
% empty, simple linear time traces are assumed otherwise FIFO arrival times
% (a single array of ascending times which describe the detection of an
% event, e.g. the arrival time of a photoelectron).
%
%       h               The correlator handle
%       trace           The piece(s) of trace to be added.
%                       i)   single trace for autocorrelation
%                            (double array)
%                       ii)  two traces for crosscorrelation
%                            (two-component double cell array)
%                       iii) macro arrival times
%                            (double array of arrival times)
%       channels        The channel indices if in arrival time mode
%                       (double array, 0 for channel 0, > 0 for channel 1)

if (nargin < 2 || nargin > 3)
    error('wrong number of parameters.')
elseif (nargin == 2)
    omexFcs(2, h, trace);
elseif (nargin == 3)
    omexFcs(2, h, trace, channels);
end
end