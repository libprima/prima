function yw = year_week(time_zone)
% YEAR_WEEK returns the following number defined by the current year and week numbers for a given
% time zone:
% yw = 100*year_number + week_number,
% where year_number and week_number are both two-digit integers.

if nargin < 1
    time_zone = 'Asia/Shanghai';
end

dt = datetime('now', 'TimeZone', time_zone);
yw = 100*mod(year(dt), 100) + week(dt);
