function pinfo = GetPe2lyrInfo(pinfo_in,fname,routine)
%% GetPe2lyrInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetPe2lyrInfo(pinfo_in,fname,routine);
%
% pinfo_in  Name of existing pinfo struct, e.g. output from CheckModelCompatibility
% fname     Name of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo  = pinfo_in;
model  = ncreadatt(fname, '/', 'model');

if strcmpi(model,'pe2lyr') ~= 1
   error('Not so fast, this is not a pe2lyr model.')
end
if strcmpi(pinfo.model,'pe2lyr') ~= 1
   error('Not so fast, this is not a pe2lyr plotting structure.')
end

pinfo.levels = ncread(fname,'lev');
pinfo.lons   = ncread(fname,'lon');
pinfo.lats   = ncread(fname,'lat');

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(      pinfo.vars, pinfo.vars{1});
      [level, lvlind] = GetLevel(    pinfo.levels);
      [lat  , latind] = GetLatitude( pinfo.lats);
      [lon  , lonind] = GetLongitude(pinfo.lons);

      pinfo.fname      = fname;
      pinfo.var        = pgvar;
      pinfo.level      = level;
      pinfo.levelindex = lvlind;
      pinfo.longitude  = lon;
      pinfo.lonindex   = lonind;
      pinfo.latitude   = lat;
      pinfo.latindex   = latind;

   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(      pinfo.vars);
      [base_time, base_tmeind] = GetTime(     pinfo.time);
      [base_lvl,  base_lvlind] = GetLevel(    pinfo.levels);
      [base_lat,  base_latind] = GetLatitude( pinfo.lats);
      [base_lon,  base_lonind] = GetLongitude(pinfo.lons);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(  pinfo.vars,   base_var);
      [comp_lvl, comp_lvlind] = GetLevel(pinfo.levels, base_lvl);

      pinfo.fname       = fname;
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;

   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(      pinfo.vars);
      [base_time, base_tmeind] = GetTime(     pinfo.time);
      [base_lvl , base_lvlind] = GetLevel(    pinfo.levels);
      [base_lat , base_latind] = GetLatitude( pinfo.lats);
      [base_lon , base_lonind] = GetLongitude(pinfo.lons);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(      pinfo.vars,   base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    pinfo.levels, base_lvl);
      [comp_lat, comp_latind] = GetLatitude( pinfo.lats,   base_lat);
      [comp_lon, comp_lonind] = GetLongitude(pinfo.lons,   base_lon);

      pinfo.fname       = fname;
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;
      pinfo.comp_lat    = comp_lat;
      pinfo.comp_latind = comp_latind;
      pinfo.comp_lon    = comp_lon;
      pinfo.comp_lonind = comp_lonind;

   case 'plotsawtooth'

       pgvar          = GetVar(       pinfo.vars);
      [level, lvlind] = GetLevel(     pinfo.levels);
      [lat  , latind] = GetLatitude(  pinfo.lats);
      [lon  , lonind] = GetLongitude( pinfo.lons);
      copyindices     = SetCopyID(fname);
      copy            = length(copyindices);

      pinfo.var_names   = pgvar;
      pinfo.truth_file  = [];
      pinfo.level       = level;
      pinfo.levelindex  = lvlind;
      pinfo.latitude    = lat;
      pinfo.latindex    = latind;
      pinfo.longitude   = lon;
      pinfo.lonindex    = lonind;
      pinfo.copies      = copy;
      pinfo.copyindices = copyindices;

   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(       pinfo.vars);
      [var1_lvl, var1_lvlind] = GetLevel(     pinfo.levels);
      [var1_lat, var1_latind] = GetLatitude(  pinfo.lats);
      [var1_lon, var1_lonind] = GetLongitude( pinfo.lons);

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(       pinfo.vars,   var1    );
      [var2_lvl, var2_lvlind] = GetLevel(     pinfo.levels, var1_lvl);
      [var2_lat, var2_latind] = GetLatitude(  pinfo.lats,   var1_lat);
      [var2_lon, var2_lonind] = GetLongitude( pinfo.lons,   var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(       pinfo.vars,   var1    );
      [var3_lvl, var3_lvlind] = GetLevel(     pinfo.levels, var1_lvl);
      [var3_lat, var3_latind] = GetLatitude(  pinfo.lats,   var1_lat);
      [var3_lon, var3_lonind] = GetLongitude( pinfo.lons,   var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.fname       = fname;
      pinfo.var1name    = var1;
      pinfo.var2name    = var2;
      pinfo.var3name    = var3;
      pinfo.var1_lvl    = var1_lvl;
      pinfo.var1_lvlind = var1_lvlind;
      pinfo.var1_lat    = var1_lat;
      pinfo.var1_latind = var1_latind;
      pinfo.var1_lon    = var1_lon;
      pinfo.var1_lonind = var1_lonind;
      pinfo.var2_lvl    = var2_lvl;
      pinfo.var2_lvlind = var2_lvlind;
      pinfo.var2_lat    = var2_lat;
      pinfo.var2_latind = var2_latind;
      pinfo.var2_lon    = var2_lon;
      pinfo.var2_lonind = var2_lonind;
      pinfo.var3_lvl    = var3_lvl;
      pinfo.var3_lvlind = var3_lvlind;
      pinfo.var3_lat    = var3_lat;
      pinfo.var3_latind = var3_latind;
      pinfo.var3_lon    = var3_lon;
      pinfo.var3_lonind = var3_lonind;
      pinfo.ens_mem     = ens_mem;
      pinfo.ltype       = ltype;

   otherwise

end

function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = 'u'; end

str = sprintf(' %s ',prognostic_vars{1});
for i = 2:length(prognostic_vars),
   str = sprintf(' %s %s ',str,prognostic_vars{i});
end
fprintf('Default variable is ''%s'', if this is OK, <cr>;\n',pgvar)
fprintf('If not, please enter one of: %s\n',str)
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), pgvar = deblank(varstring); end


function [time, timeind] = GetTime(times, deftime)
%----------------------------------------------------------------------
if (nargin == 2), time = deftime; else time = times(1); end

fprintf('Default time is %f, if this is OK, <cr>;\n',time)
fprintf('If not, enter a time between %.4f and %.4f, we use the closest.\n', ...
                         min(times),max(times))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), time  = str2num(varstring); end

d       = abs(time - times);  % crude distance
ind     = find(min(d) == d);  % multiple minima possible
timeind = ind(1);             % use the first one
time    = times(timeind);


function [level, lvlind] = GetLevel(levels, deflevel)
%----------------------------------------------------------------------
% level and lvlind will not be equal for all models, (and probably
% shouldn't for the pe2lyr ... but for future expansion ...
if (nargin == 2), level = deflevel; else level = 1; end

fprintf('Default level is  %d, if this is OK, <cr>;\n',level)
fprintf('If not, enter a level between %d and %d, inclusive ...\n', ...
                         min(levels),max(levels))
varstring = input('we''ll use the closest (no syntax required)\n','s');

if ~isempty(varstring), level = str2num(varstring); end

d      = abs(level - levels);  % crude distance
ind    = find(min(d) == d);    % multiple minima possible
lvlind = ind(1);               % use the first one
level  = levels(lvlind);


function [lon, lonind] = GetLongitude(lons, deflon)
%----------------------------------------------------------------------
if (nargin == 2), lon = deflon; else lon = 255.0; end

fprintf('Default longitude is %f, if this is OK, <cr>;\n',lon)
fprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lons),max(lons))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2num(varstring); end

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
lonind = ind(1);             % use the first one
lon    = lons(lonind);


function [lat, latind] = GetLatitude(lats, deflat)
%----------------------------------------------------------------------
if (nargin == 2), lat = deflat; else lat = 40.0; end

fprintf('Default latitude is %f, if this is OK, <cr>;\n',lat)
fprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lats),max(lats))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lat = str2num(varstring); end

d      = abs(lat - lats);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
latind = ind(1);             % use the first one
lat    = lats(latind);


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
