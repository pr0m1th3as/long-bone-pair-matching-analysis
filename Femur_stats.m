## A script for analyzing the descriptive statistics
## of pair-matched Femur skeletal samples
pkg load csg-toolkit
Header = longbone_Measurements;

## Load matched and mismatched Femur data
load ("Femur-Data.mat");
F_pairs = F_match(:,1:2);
F_match(:,1:2) = [];
F_mismatch(:,1:2) = [];

## Get number of matched and mismatched cases
cases_m = size (F_match, 1);
cases_mm = size (F_mismatch, 1);

## Compute mean and standard deviation for each measurement diff
Femur_mu = mean (F_match, 1);
Femur_sd = std (F_match, 1);
Femur_lb = Femur_sd;
Femur_ub = Femur_sd;
Femur_id = ones (size (Femur_sd));

## Calculate bin centers to +- 7 STDs
bins = ([-70:70] * 0.1)' * Femur_sd + Femur_mu;
## Set percentage threshold below which a variable is ignored
threshold = 17;

## Handle Max Distance and Perimeter measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Max Distance", "Perimeter 20%", "Perimeter 35%", ...
            "Perimeter 50%", "Perimeter 65%", "Perimeter 80%"};
x_limits = [-1, 1] .* [30; 40; 15; 20; 15; 20];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [60000; 40000; 40000; 40000; 40000; 40000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.1f and %0.1f mm", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 6, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 6, v + 6);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig1-MaxD_Per.png");
close h1


## Handle Area measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Area 20%", "Area 35%", "Area 50%", "Area 65%", "Area 80%"};
x_limits = [-1, 1] .* [400; 150; 150; 200; 250];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [40000; 40000; 40000; 40000; 40000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.1f and %0.1f mm^2", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^2", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^2", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig2-Area.png");
close h1


## Handle ArPerIndex measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"ArPerIndex 20%", "ArPerIndex 35%", ...
            "ArPerIndex 50%", "ArPerIndex 65%", "ArPerIndex 80%"};
x_limits = [-1, 1] .* [0.35; 0.15; 0.2; 0.15; 0.1];
y_limit1 = [0, 1] .* [50; 50; 50; 50; 50];
y_limit2 = [0, 1] .* [7000; 7000; 7000; 7000; 7000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.2f and %0.2f", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig3-ArPerIndex.png");
close h1


## Handle Ix measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Ix 20%", "Ix 35%", "Ix 50%", "Ix 65%", "Ix 80%"};
x_limits = [-1, 1] .* [60000; 16000; 16000; 16000; 32000];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [30000; 30000; 30000; 30000; 30000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.0f and %0.0f mm^4", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  xt = get (gca, "xtick");
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig4-Ix.png");
close h1


## Handle Iy measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Iy 20%", "Iy 35%", "Iy 50%", "Iy 65%", "Iy 80%"};
x_limits = [-1, 1] .* [90000; 20000; 20000; 20000; 50000];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [30000; 30000; 30000; 30000; 30000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.0f and %0.0f mm^4", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  xt = get (gca, "xtick");
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig5-Iy.png");
close h1


## Handle Ixy measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Ixy 20%", "Ixy 35%", "Ixy 50%", "Ixy 65%", "Ixy 80%"};
x_limits = [-1, 1] .* [32000; 10000; 10000; 8000; 16000];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [8000; 8000; 8000; 8000; 8000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.0f and %0.0f mm^4", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  xt = get (gca, "xtick");
  set (gca, "xticklabel", xt * 0.001);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^3)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  set (gca, "xticklabel", xt * 0.001);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^3)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig6-Ixy.png");
close h1


## Handle Ix/Iy measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Ix/Iy 20%", "Ix/Iy 35%", "Ix/Iy 50%", "Ix/Iy 65%", "Ix/Iy 80%"};
x_limits = [-1, 1] .* [1; 1; 1; 1; 0.5];
y_limit1 = [0, 1] .* [30; 30; 30; 30; 30];
y_limit2 = [0, 1] .* [4000; 4000; 4000; 4000; 4000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.2f and %0.2f", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig7-Ix_Iy.png");
close h1


## Handle Imin measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Imin 20%", "Imin 35%", "Imin 50%", "Imin 65%", "Imin 80%"};
x_limits = [-1, 1] .* [60000; 16000; 16000; 16000; 32000];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [30000; 30000; 30000; 30000; 30000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.0f and %0.0f mm^4", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  xt = get (gca, "xtick");
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig8-Imin.png");
close h1


## Handle Imax measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Imax 20%", "Imax 35%", "Imax 50%", "Imax 65%", "Imax 80%"};
x_limits = [-1, 1] .* [90000; 20000; 20000; 20000; 50000];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [30000; 30000; 30000; 30000; 30000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.0f and %0.0f mm^4", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  xt = get (gca, "xtick");
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  set (gca, "xticklabel", xt * 0.0001);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in mm^4 (x10^4)", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig9-Imax.png");
close h1


## Handle Imax/Imin measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Imax/Imin 20%", "Imax/Imin 35%", ...
            "Imax/Imin 50%", "Imax/Imin 65%", "Imax/Imin 80%"};
x_limits = [-1, 1] .* [1.5; 1; 1; 1; 1];
y_limit1 = [0, 1] .* [30; 30; 30; 30; 30];
y_limit2 = [0, 1] .* [10000; 10000; 10000; 10000; 10000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.2f and %0.2f", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("unitless difference", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig10-Imax_Imin.png");
close h1


## Handle theta measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"theta 20%", "theta 35%", "theta 50%", "theta 65%", "theta 80%"};
x_limits = [-1, 1] .* [100; 100; 100; 100; 50];
y_limit1 = [0, 1] .* [40; 40; 40; 40; 40];
y_limit2 = [0, 1] .* [6000; 6000; 6000; 6000; 6000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.1f and %0.1f degrees", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in degrees", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (tname, "fontsize", 26);
  xlabel ("difference in degrees", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig11-theta.png");
close h1


## Handle Dihedral angle measurements for each cross-section
h1 = figure ("position", [400 0 3000 1800], "visible", "off");
varNames = {"Dihedral angle 20_35", "Dihedral angle 35_50", ...
            "Dihedral angle 50_65", "Dihedral angle 65_80", "Diaphyseal Bending"};
x_limits = [-1, 1] .* [3; 4; 4; 3; 10];
y_limit1 = [0, 1] .* [30; 30; 30; 30; 30];
y_limit2 = [0, 1] .* [5000; 5000; 5000; 5000; 5000];

for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (F_match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (F_mismatch(:,varID), bins(:,varID)');
  ## Find boundary bins of matched cases
  lb_var_ind = find (ne_var_m, 1, "first");
  ub_var_ind = find (ne_var_m, 1, "last");
  ## Leave one bin gap unless already at +- 7 STDs
  if (lb_var_ind > 1)
    lb_var_ind -= 1;
  endif
  if (ub_var_ind < 141)
    ub_var_ind += 1;
  endif
  lb_var = bins(lb_var_ind, varID);
  ub_var = bins(ub_var_ind, varID);
  ## Save boundary limits
  Femur_lb(varID) = lb_var;
  Femur_ub(varID) = ub_var;
  ## Calculate definite excluded mismatched cases
  excluded_cases = (sum (ne_var_mm([1:lb_var_ind])) + ...
                    sum (ne_var_mm([ub_var_ind:end]))) / cases_mm * 100;
  ## Ignore variable if excluded_cases drop below certain percentage
  if (excluded_cases < threshold)
    Femur_id(varID) = 0;
  endif
  ## Prepare annotation strings
  txt_var_exc = sprintf ("%0.2f%% of pairs excluded", excluded_cases);
  txt_var_bnd = sprintf ("bounded: %0.1f and %0.1f degrees", lb_var, ub_var);
  ## Plot results in histograms
  subplot (2, 5, v);
  bar (bc_var_m, ne_var_m, "facecolor", "r", "edgecolor", "b");
  xlim (x_limits(v,:));
  ylim (y_limit1(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- matched", varNames{v});
  title (strrep (tname, "_", "-"), "fontsize", 26);
  xlabel ("difference in degrees", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit1(v,2) / 5.5;
  text (0, -yt, txt_var_bnd, "horizontalalignment", "center", "fontsize", 22);
  subplot (2, 5, v + 5);
  bar (bc_var_mm, ne_var_mm, "facecolor", "g", "edgecolor", "y");
  xlim (x_limits(v,:));
  ylim (y_limit2(v,:));
  set (gca, "fontsize",18);
  tname = sprintf ("%s\n- mismatched", varNames{v});
  title (strrep (tname, "_", "-"), "fontsize", 26);
  xlabel ("difference in degrees", "fontsize", 20, "fontangle", "italic");
  ylabel ("number of cases", "fontsize", 20, "fontangle", "italic");
  yt = y_limit2(v,2) / 5.5;
  text (0, -yt, txt_var_exc, "horizontalalignment", "center", "fontsize", 22);
endfor
## save plot as an image
print (h1, "Femur_fig12-dihedral.png");
close h1


## Append descriptives for pair matching into 'descriptives_pair.mat'
filename = fullfile (pwd, 'private', 'descriptives_pair.mat');
Fvars = [Femur_mu; Femur_sd; Femur_lb; Femur_ub; Femur_id];
save ('-append', '-binary', filename, 'Fvars');
clear
