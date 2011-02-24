%MT_SPLINE_SUPPORT - Calculate clamped b-spline parameters
%
%  K = MT_SPLINE_SUPPORT(POINTS,POS,DEGREE,POS_RANGE)
%
% INPUT
%   POINTS			Number of b-spline knots
%   POS             Position along b-spline 
%   DEGREE          Degree of b-spline
%   POS_RANGE       Range of possible positions
%
% OUTPUT
%   K		        B-spline basis function values (1 * npoints)
%
% DESCRIPTION
% Determines clamped b-spline values. %

% SEE ALSO
% MT_CREATE_SEQ_MODEL
%
% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function k = mt_spline_support(points,pos,degree,pos_range)

n = points - 1;
m = n + degree + 1;
k = zeros(1,n + 1);
t = [zeros(1,degree) 0:1./(m - 2 * degree):1 ones(1,degree)];

tt = (pos - pos_range(1)) / (pos_range(2) - 1);
loc = find(t <= tt);
startloc = min(loc(end),points);

k(startloc) = 1;
for i = 1:degree
  k(startloc - i) = ((t(startloc + 1) - tt) ./ (t(startloc + 1) - t(startloc - i + 1))) * k(startloc - i + 1);
  for j = (startloc - i + 1):(startloc - 1)
      k(j) = ((tt - t(j)) ./ (t(j + i) - t(j))) * k(j) + ((t(j + i + 1) - tt) ./ (t(j + i + 1) - t(j + 1))) * k(j+1);
  end;
  k(startloc) = ((tt - t(startloc)) ./ (t(startloc + i) - t(startloc))) *  k(startloc);
end;
