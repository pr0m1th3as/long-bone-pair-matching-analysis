## Process all measurement data and create matched-mismatched datasets
## for each bone and each collection, and assemble pair-matching
## descriptive statistics for interpopulation comparisons. Pool all
## populations together and create descriptives for the entire assemblage.

## Load required packages
pkg load csg-dataset
pkg load csg-toolkit

## Load CSG dataset
CSG = load ('csg-properties.mat').csg_properties;
nametoid = @(x) str2num (x([4:6]));

## Load measurement names
Header = longbone_Measurements;

## List variables in desired order for descriptive tables
varNames = {"Max Distance", "Perimeter 20%", "Perimeter 35%", ...
            "Perimeter 50%", "Perimeter 65%", "Perimeter 80%", ...
            "Area 20%", "Area 35%", "Area 50%", "Area 65%", "Area 80%", ...
            "ArPerIndex 20%", "ArPerIndex 35%", ...
            "ArPerIndex 50%", "ArPerIndex 65%", "ArPerIndex 80%", ...
            "Ix 20%", "Ix 35%", "Ix 50%", "Ix 65%", "Ix 80%", ...
            "Iy 20%", "Iy 35%", "Iy 50%", "Iy 65%", "Iy 80%", ...
            "Ixy 20%", "Ixy 35%", "Ixy 50%", "Ixy 65%", "Ixy 80%", ...
            "Ix/Iy 20%", "Ix/Iy 35%", "Ix/Iy 50%", "Ix/Iy 65%", "Ix/Iy 80%", ...
            "Imin 20%", "Imin 35%", "Imin 50%", "Imin 65%", "Imin 80%", ...
            "Imax 20%", "Imax 35%", "Imax 50%", "Imax 65%", "Imax 80%", ...
            "Imax/Imin 20%", "Imax/Imin 35%", ...
            "Imax/Imin 50%", "Imax/Imin 65%", "Imax/Imin 80%", ...
            "theta 20%", "theta 35%", "theta 50%", "theta 65%", "theta 80%", ...
            "Dihedral angle 20_35", "Dihedral angle 35_50", ...
            "Dihedral angle 50_65", "Dihedral angle 65_80", "Diaphyseal Bending"};


## Gather Femur samples from each collection into left and right sides
F_ALL = [];
for c = 1:3
  Femur = [];

  ## Left side
  data = CSG{2, (c - 1) * 8 + 1};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Femur = [Femur; sIdx, side, collection, measurements];

  ## Right side
  data = CSG{2, (c - 1) * 8 + 2};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Femur = [Femur; sIdx, side, collection, measurements];

  ## Append collection to pooled data
  F_ALL = [F_ALL; Femur];

  ## Split dataset into left and right side
  Femur_L = Femur(Femur(:,2) == 1, [6:66]);
  Femur_R = Femur(Femur(:,2) == 2, [6:66]);
  F_idx_L = Femur(Femur(:,2) == 1, [1,3:5]);
  F_idx_R = Femur(Femur(:,2) == 2, [1,3:5]);

  ## Compute absolute values
  Lvec = abs (Femur_L);
  Rvec = abs (Femur_R);

  ## Scan through all possible pairs and assemble matched/mismatched diffs
  ## betweeen left and right sides
  match = [];
  mismatch = [];
  for ii = 1:size (Femur_L, 1)
    ## Get absolute differences between left and right side
    Dvec = Lvec(ii,:) - Rvec;
    ## Get index of matching pair
    idx = all (F_idx_L(ii,1:3) == F_idx_R(:,1:3), 2);
    if (any (idx))
      match = [match; Dvec(idx,:)];
    endif
    ## Get indices of mismatched pairs
    idx = ! idx;
    mismatch = [mismatch; Dvec(idx,:)];
  endfor

  ## Get number of matched and mismatched cases
  F_cases_m(c) = size (match, 1);
  F_cases_mm(c) = size (mismatch, 1);

  ## Calculate bin centers to +- 7 STDs
  bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

  ## Compute descriptives for each population
  for v = 1:numel (varNames)
    ## Get column index in data
    varID = find (strcmpi (Header, varNames{v}));
    ## Compute mean and standard deviation
    FLmu(v,c) = mean (Femur_L(:,varID));
    FLsd(v,c) = std (Femur_L(:,varID));
    FRmu(v,c) = mean (Femur_R(:,varID));
    FRsd(v,c) = std (Femur_R(:,varID));
    F_mu(v,c) = mean (match(:,varID));
    F_sd(v,c) = std (match(:,varID));
    ## Compute histogram bins
    [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
    [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
    F_lb(v,c) = lb_var;
    F_ub(v,c) = ub_var;
    ## Calculate definite excluded mismatched cases
    F_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
                 sum (ne_var_mm([ub_var_ind:end]))) / F_cases_mm(c) * 100;
  endfor

  ## Clear unnecessary variables for next population
  clear -x CSG nametoid F_ALL F_mu F_sd F_lb F_ub F_er ...
           F_cases_m F_cases_mm Header varNames FLmu FLsd FRmu FRsd
endfor

## Process pooled collections
c = 4;

## Split dataset into left and right side
Femur_L = F_ALL(F_ALL(:,2) == 1, [6:66]);
Femur_R = F_ALL(F_ALL(:,2) == 2, [6:66]);
F_idx_L = F_ALL(F_ALL(:,2) == 1, [1,3:5]);
F_idx_R = F_ALL(F_ALL(:,2) == 2, [1,3:5]);

## Compute absolute values
Lvec = abs (Femur_L);
Rvec = abs (Femur_R);

## Scan through all possible pairs and assemble matched/mismatched diffs
## betweeen left and right sides
match = [];
mismatch = [];
for ii = 1:size (Femur_L, 1)
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (F_idx_L(ii,1:3) == F_idx_R(:,1:3), 2);
  if (any (idx))
    match = [match; Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  mismatch = [mismatch; Dvec(idx,:)];
endfor

## Get number of matched and mismatched cases
F_cases_m(c) = size (match, 1);
F_cases_mm(c) = size (mismatch, 1);

## Calculate bin centers to +- 7 STDs
bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

## Compute descriptives for pooled collections
for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute mean and standard deviation
  FLmu(v,c) = mean (Femur_L(:,varID));
  FLsd(v,c) = std (Femur_L(:,varID));
  FRmu(v,c) = mean (Femur_R(:,varID));
  FRsd(v,c) = std (Femur_R(:,varID));
  F_mu(v,c) = mean (match(:,varID));
  F_sd(v,c) = std (match(:,varID));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
  F_lb(v,c) = lb_var;
  F_ub(v,c) = ub_var;
  ## Calculate definite excluded mismatched cases
  F_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
               sum (ne_var_mm([ub_var_ind:end]))) / F_cases_mm(c) * 100;
endfor

## Assemble stats into tables for each collection and pooled sample
F1table = table (table (FLmu(:,1), FLsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (FRmu(:,1), FRsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (F_mu(:,1), F_sd(:,1), F_lb(:,1), F_ub(:,1), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (F_er(:,1), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames');
F2table = table (table (FLmu(:,2), FLsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (FRmu(:,2), FRsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (F_mu(:,2), F_sd(:,2), F_lb(:,2), F_ub(:,2), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (F_er(:,2), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
F3table = table (table (FLmu(:,3), FLsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (FRmu(:,3), FRsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (F_mu(:,3), F_sd(:,3), F_lb(:,3), F_ub(:,3), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (F_er(:,3), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
F4table = table (table (FLmu(:,4), FLsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (FRmu(:,4), FRsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (F_mu(:,4), F_sd(:,4), F_lb(:,4), F_ub(:,4), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (F_er(:,4), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);

## Save tables to CSV files
table2csv (F1table, "Femur-Stas-Athens.csv");
table2csv (F2table, "Femur-Stas-Crete.csv");
table2csv (F3table, "Femur-Stas-Granada.csv");
table2csv (F4table, "Femur-Stas-Pooled.csv");

## Clear unnecessary variables
clear -x CSG nametoid Header varNames F1table F2table F3table F4table



## Gather Humerus samples from each collection into left and right sides
H_ALL = [];
for c = 1:3
  Humerus = [];

  ## Left side
  data = CSG{2, (c - 1) * 8 + 1};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Humerus = [Humerus; sIdx, side, collection, measurements];

  ## Right side
  data = CSG{2, (c - 1) * 8 + 2};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Humerus = [Humerus; sIdx, side, collection, measurements];

  ## Append collection to pooled data
  H_ALL = [H_ALL; Humerus];

  ## Split dataset into left and right side
  Humerus_L = Humerus(Humerus(:,2) == 1, [6:66]);
  Humerus_R = Humerus(Humerus(:,2) == 2, [6:66]);
  H_idx_L = Humerus(Humerus(:,2) == 1, [1,3:5]);
  H_idx_R = Humerus(Humerus(:,2) == 2, [1,3:5]);

  ## Compute absolute values
  Lvec = abs (Humerus_L);
  Rvec = abs (Humerus_R);

  ## Scan through all possible pairs and assemble matched/mismatched diffs
  ## betweeen left and right sides
  match = [];
  mismatch = [];
  for ii = 1:size (Humerus_L, 1)
    ## Get absolute differences between left and right side
    Dvec = Lvec(ii,:) - Rvec;
    ## Get index of matching pair
    idx = all (H_idx_L(ii,1:3) == H_idx_R(:,1:3), 2);
    if (any (idx))
      match = [match; Dvec(idx,:)];
    endif
    ## Get indices of mismatched pairs
    idx = ! idx;
    mismatch = [mismatch; Dvec(idx,:)];
  endfor

  ## Get number of matched and mismatched cases
  H_cases_m(c) = size (match, 1);
  H_cases_mm(c) = size (mismatch, 1);

  ## Calculate bin centers to +- 7 STDs
  bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

  ## Compute descriptives for each population
  for v = 1:numel (varNames)
    ## Get column index in data
    varID = find (strcmpi (Header, varNames{v}));
    ## Compute mean and standard deviation
    HLmu(v,c) = mean (Humerus_L(:,varID));
    HLsd(v,c) = std (Humerus_L(:,varID));
    HRmu(v,c) = mean (Humerus_R(:,varID));
    HRsd(v,c) = std (Humerus_R(:,varID));
    H_mu(v,c) = mean (match(:,varID));
    H_sd(v,c) = std (match(:,varID));
    ## Compute histogram bins
    [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
    [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
    H_lb(v,c) = lb_var;
    H_ub(v,c) = ub_var;
    ## Calculate definite excluded mismatched cases
    H_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
                 sum (ne_var_mm([ub_var_ind:end]))) / H_cases_mm(c) * 100;
  endfor

  ## Clear unnecessary variables for next population
  clear -x CSG nametoid H_ALL H_mu H_sd H_lb H_ub H_er ...
           H_cases_m H_cases_mm Header varNames HLmu HLsd HRmu HRsd ...
           F1table F2table F3table F4table
endfor

## Process pooled collections
c = 4;

## Split dataset into left and right side
Humerus_L = H_ALL(H_ALL(:,2) == 1, [6:66]);
Humerus_R = H_ALL(H_ALL(:,2) == 2, [6:66]);
H_idx_L = H_ALL(H_ALL(:,2) == 1, [1,3:5]);
H_idx_R = H_ALL(H_ALL(:,2) == 2, [1,3:5]);

## Compute absolute values
Lvec = abs (Humerus_L);
Rvec = abs (Humerus_R);

## Scan through all possible pairs and assemble matched/mismatched diffs
## betweeen left and right sides
match = [];
mismatch = [];
for ii = 1:size (Humerus_L, 1)
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (H_idx_L(ii,1:3) == H_idx_R(:,1:3), 2);
  if (any (idx))
    match = [match; Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  mismatch = [mismatch; Dvec(idx,:)];
endfor

## Get number of matched and mismatched cases
H_cases_m(c) = size (match, 1);
H_cases_mm(c) = size (mismatch, 1);

## Calculate bin centers to +- 7 STDs
bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

## Compute descriptives for pooled collections
for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute mean and standard deviation
  HLmu(v,c) = mean (Humerus_L(:,varID));
  HLsd(v,c) = std (Humerus_L(:,varID));
  HRmu(v,c) = mean (Humerus_R(:,varID));
  HRsd(v,c) = std (Humerus_R(:,varID));
  H_mu(v,c) = mean (match(:,varID));
  H_sd(v,c) = std (match(:,varID));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
  H_lb(v,c) = lb_var;
  H_ub(v,c) = ub_var;
  ## Calculate definite excluded mismatched cases
  H_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
               sum (ne_var_mm([ub_var_ind:end]))) / H_cases_mm(c) * 100;
endfor

## Assemble stats into tables for each collection and pooled sample
H1table = table (table (HLmu(:,1), HLsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (HRmu(:,1), HRsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (H_mu(:,1), H_sd(:,1), H_lb(:,1), H_ub(:,1), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (H_er(:,1), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
H2table = table (table (HLmu(:,2), HLsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (HRmu(:,2), HRsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (H_mu(:,2), H_sd(:,2), H_lb(:,2), H_ub(:,2), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (H_er(:,2), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
H3table = table (table (HLmu(:,3), HLsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (HRmu(:,3), HRsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (H_mu(:,3), H_sd(:,3), H_lb(:,3), H_ub(:,3), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (H_er(:,3), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
H4table = table (table (HLmu(:,4), HLsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (HRmu(:,4), HRsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (H_mu(:,4), H_sd(:,4), H_lb(:,4), H_ub(:,4), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (H_er(:,4), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);

## Save tables to CSV files
table2csv (H1table, "Humerus-Stas-Athens.csv");
table2csv (H2table, "Humerus-Stas-Crete.csv");
table2csv (H3table, "Humerus-Stas-Granada.csv");
table2csv (H4table, "Humerus-Stas-Pooled.csv");

## Clear unnecessary variables
clear -x CSG nametoid Header varNames F1table F2table F3table F4table ...
         H1table H2table H3table H4table



## Gather Tibia samples from each collection into left and right sides
T_ALL = [];
for c = 1:3
  Tibia = [];

  ## Left side
  data = CSG{2, (c - 1) * 8 + 1};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Tibia = [Tibia; sIdx, side, collection, measurements];

  ## Right side
  data = CSG{2, (c - 1) * 8 + 2};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Tibia = [Tibia; sIdx, side, collection, measurements];

  ## Append collection to pooled data
  T_ALL = [T_ALL; Tibia];

  ## Split dataset into left and right side
  Tibia_L = Tibia(Tibia(:,2) == 1, [6:66]);
  Tibia_R = Tibia(Tibia(:,2) == 2, [6:66]);
  T_idx_L = Tibia(Tibia(:,2) == 1, [1,3:5]);
  T_idx_R = Tibia(Tibia(:,2) == 2, [1,3:5]);

  ## Compute absolute values
  Lvec = abs (Tibia_L);
  Rvec = abs (Tibia_R);

  ## Scan through all possible pairs and assemble matched/mismatched diffs
  ## betweeen left and right sides
  match = [];
  mismatch = [];
  for ii = 1:size (Tibia_L, 1)
    ## Get absolute differences between left and right side
    Dvec = Lvec(ii,:) - Rvec;
    ## Get index of matching pair
    idx = all (T_idx_L(ii,1:3) == T_idx_R(:,1:3), 2);
    if (any (idx))
      match = [match; Dvec(idx,:)];
    endif
    ## Get indices of mismatched pairs
    idx = ! idx;
    mismatch = [mismatch; Dvec(idx,:)];
  endfor

  ## Get number of matched and mismatched cases
  T_cases_m(c) = size (match, 1);
  T_cases_mm(c) = size (mismatch, 1);

  ## Calculate bin centers to +- 7 STDs
  bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

  ## Compute descriptives for each population
  for v = 1:numel (varNames)
    ## Get column index in data
    varID = find (strcmpi (Header, varNames{v}));
    ## Compute mean and standard deviation
    TLmu(v,c) = mean (Tibia_L(:,varID));
    TLsd(v,c) = std (Tibia_L(:,varID));
    TRmu(v,c) = mean (Tibia_R(:,varID));
    TRsd(v,c) = std (Tibia_R(:,varID));
    T_mu(v,c) = mean (match(:,varID));
    T_sd(v,c) = std (match(:,varID));
    ## Compute histogram bins
    [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
    [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
    T_lb(v,c) = lb_var;
    T_ub(v,c) = ub_var;
    ## Calculate definite excluded mismatched cases
    T_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
                 sum (ne_var_mm([ub_var_ind:end]))) / T_cases_mm(c) * 100;
  endfor

  ## Clear unnecessary variables for next population
  clear -x CSG nametoid T_ALL T_mu T_sd T_lb T_ub T_er ...
           T_cases_m T_cases_mm Header varNames TLmu TLsd TRmu TRsd ...
           F1table F2table F3table F4table H1table H2table H3table H4table
endfor

## Process pooled collections
c = 4;

## Split dataset into left and right side
Tibia_L = T_ALL(T_ALL(:,2) == 1, [6:66]);
Tibia_R = T_ALL(T_ALL(:,2) == 2, [6:66]);
T_idx_L = T_ALL(T_ALL(:,2) == 1, [1,3:5]);
T_idx_R = T_ALL(T_ALL(:,2) == 2, [1,3:5]);

## Compute absolute values
Lvec = abs (Tibia_L);
Rvec = abs (Tibia_R);

## Scan through all possible pairs and assemble matched/mismatched diffs
## betweeen left and right sides
match = [];
mismatch = [];
for ii = 1:size (Tibia_L, 1)
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (T_idx_L(ii,1:3) == T_idx_R(:,1:3), 2);
  if (any (idx))
    match = [match; Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  mismatch = [mismatch; Dvec(idx,:)];
endfor

## Get number of matched and mismatched cases
T_cases_m(c) = size (match, 1);
T_cases_mm(c) = size (mismatch, 1);

## Calculate bin centers to +- 7 STDs
bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

## Compute descriptives for pooled collections
for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute mean and standard deviation
  TLmu(v,c) = mean (Tibia_L(:,varID));
  TLsd(v,c) = std (Tibia_L(:,varID));
  TRmu(v,c) = mean (Tibia_R(:,varID));
  TRsd(v,c) = std (Tibia_R(:,varID));
  T_mu(v,c) = mean (match(:,varID));
  T_sd(v,c) = std (match(:,varID));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
  T_lb(v,c) = lb_var;
  T_ub(v,c) = ub_var;
  ## Calculate definite excluded mismatched cases
  T_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
               sum (ne_var_mm([ub_var_ind:end]))) / T_cases_mm(c) * 100;
endfor

## Assemble stats into tables for each collection and pooled sample
T1table = table (table (TLmu(:,1), TLsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (TRmu(:,1), TRsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (T_mu(:,1), T_sd(:,1), T_lb(:,1), T_ub(:,1), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (T_er(:,1), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
T2table = table (table (TLmu(:,2), TLsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (TRmu(:,2), TRsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (T_mu(:,2), T_sd(:,2), T_lb(:,2), T_ub(:,2), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (T_er(:,2), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
T3table = table (table (TLmu(:,3), TLsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (TRmu(:,3), TRsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (T_mu(:,3), T_sd(:,3), T_lb(:,3), T_ub(:,3), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (T_er(:,3), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
T4table = table (table (TLmu(:,4), TLsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (TRmu(:,4), TRsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (T_mu(:,4), T_sd(:,4), T_lb(:,4), T_ub(:,4), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (T_er(:,4), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);

## Save tables to CSV files
table2csv (T1table, "Tibia-Stas-Athens.csv");
table2csv (T2table, "Tibia-Stas-Crete.csv");
table2csv (T3table, "Tibia-Stas-Granada.csv");
table2csv (T4table, "Tibia-Stas-Pooled.csv");

## Clear unnecessary variables
clear -x CSG nametoid Header varNames F1table F2table F3table F4table ...
         H1table H2table H3table H4table T1table T2table T3table T4table



## Gather Ulna samples from each collection into left and right sides
U_ALL = [];
for c = 1:3
  Ulna = [];

  ## Left side
  data = CSG{2, (c - 1) * 8 + 1};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Ulna = [Ulna; sIdx, side, collection, measurements];

  ## Right side
  data = CSG{2, (c - 1) * 8 + 2};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Ulna = [Ulna; sIdx, side, collection, measurements];

  ## Append collection to pooled data
  U_ALL = [U_ALL; Ulna];

  ## Split dataset into left and right side
  Ulna_L = Ulna(Ulna(:,2) == 1, [6:66]);
  Ulna_R = Ulna(Ulna(:,2) == 2, [6:66]);
  U_idx_L = Ulna(Ulna(:,2) == 1, [1,3:5]);
  U_idx_R = Ulna(Ulna(:,2) == 2, [1,3:5]);

  ## Compute absolute values
  Lvec = abs (Ulna_L);
  Rvec = abs (Ulna_R);

  ## Scan through all possible pairs and assemble matched/mismatched diffs
  ## betweeen left and right sides
  match = [];
  mismatch = [];
  for ii = 1:size (Ulna_L, 1)
    ## Get absolute differences between left and right side
    Dvec = Lvec(ii,:) - Rvec;
    ## Get index of matching pair
    idx = all (U_idx_L(ii,1:3) == U_idx_R(:,1:3), 2);
    if (any (idx))
      match = [match; Dvec(idx,:)];
    endif
    ## Get indices of mismatched pairs
    idx = ! idx;
    mismatch = [mismatch; Dvec(idx,:)];
  endfor

  ## Get number of matched and mismatched cases
  U_cases_m(c) = size (match, 1);
  U_cases_mm(c) = size (mismatch, 1);

  ## Calculate bin centers to +- 7 STDs
  bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

  ## Compute descriptives for each population
  for v = 1:numel (varNames)
    ## Get column index in data
    varID = find (strcmpi (Header, varNames{v}));
    ## Compute mean and standard deviation
    ULmu(v,c) = mean (Ulna_L(:,varID));
    ULsd(v,c) = std (Ulna_L(:,varID));
    URmu(v,c) = mean (Ulna_R(:,varID));
    URsd(v,c) = std (Ulna_R(:,varID));
    U_mu(v,c) = mean (match(:,varID));
    U_sd(v,c) = std (match(:,varID));
    ## Compute histogram bins
    [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
    [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
    U_lb(v,c) = lb_var;
    U_ub(v,c) = ub_var;
    ## Calculate definite excluded mismatched cases
    U_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
                 sum (ne_var_mm([ub_var_ind:end]))) / U_cases_mm(c) * 100;
  endfor

  ## Clear variables for next population
  clear -x CSG nametoid U_ALL U_mu U_sd U_lb U_ub U_er ...
           U_cases_m U_cases_mm Header varNames ULmu ULsd URmu URsd ...
           F1table F2table F3table F4table H1table H2table H3table H4table ...
           T1table T2table T3table T4table
endfor

## Process pooled collections
c = 4;

## Split dataset into left and right side
Ulna_L = U_ALL(U_ALL(:,2) == 1, [6:66]);
Ulna_R = U_ALL(U_ALL(:,2) == 2, [6:66]);
U_idx_L = U_ALL(U_ALL(:,2) == 1, [1,3:5]);
U_idx_R = U_ALL(U_ALL(:,2) == 2, [1,3:5]);

## Compute absolute values
Lvec = abs (Ulna_L);
Rvec = abs (Ulna_R);

## Scan through all possible pairs and assemble matched/mismatched diffs
## betweeen left and right sides
match = [];
mismatch = [];
for ii = 1:size (Ulna_L, 1)
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (U_idx_L(ii,1:3) == U_idx_R(:,1:3), 2);
  if (any (idx))
    match = [match; Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  mismatch = [mismatch; Dvec(idx,:)];
endfor

## Get number of matched and mismatched cases
U_cases_m(c) = size (match, 1);
U_cases_mm(c) = size (mismatch, 1);

## Calculate bin centers to +- 7 STDs
bins = ([-70:70] * 0.1)' * std (match, 1) + mean (match, 1);

## Compute descriptives for pooled collections
for v = 1:numel (varNames)
  ## Get column index in data
  varID = find (strcmpi (Header, varNames{v}));
  ## Compute mean and standard deviation
  ULmu(v,c) = mean (Ulna_L(:,varID));
  ULsd(v,c) = std (Ulna_L(:,varID));
  URmu(v,c) = mean (Ulna_R(:,varID));
  URsd(v,c) = std (Ulna_R(:,varID));
  U_mu(v,c) = mean (match(:,varID));
  U_sd(v,c) = std (match(:,varID));
  ## Compute histogram bins
  [ne_var_m, bc_var_m] = hist (match(:,varID), bins(:,varID)');
  [ne_var_mm, bc_var_mm] = hist (mismatch(:,varID), bins(:,varID)');
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
  U_lb(v,c) = lb_var;
  U_ub(v,c) = ub_var;
  ## Calculate definite excluded mismatched cases
  U_er(v,c) = (sum (ne_var_mm([1:lb_var_ind])) + ...
               sum (ne_var_mm([ub_var_ind:end]))) / U_cases_mm(c) * 100;
endfor

## Assemble stats into tables for each collection and pooled sample
U1table = table (table (ULmu(:,1), ULsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (URmu(:,1), URsd(:,1), "VariableNames", {'mean', 'std'}), ...
                 table (U_mu(:,1), U_sd(:,1), U_lb(:,1), U_ub(:,1), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (U_er(:,1), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
U2table = table (table (ULmu(:,2), ULsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (URmu(:,2), URsd(:,2), "VariableNames", {'mean', 'std'}), ...
                 table (U_mu(:,2), U_sd(:,2), U_lb(:,2), U_ub(:,2), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (U_er(:,2), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
U3table = table (table (ULmu(:,3), ULsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (URmu(:,3), URsd(:,3), "VariableNames", {'mean', 'std'}), ...
                 table (U_mu(:,3), U_sd(:,3), U_lb(:,3), U_ub(:,3), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (U_er(:,3), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);
U4table = table (table (ULmu(:,4), ULsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (URmu(:,4), URsd(:,4), "VariableNames", {'mean', 'std'}), ...
                 table (U_mu(:,4), U_sd(:,4), U_lb(:,4), U_ub(:,4), "VariableNames", ...
                 {'mean', 'std', 'lower bound', 'upper bound'}), ...
                 table (U_er(:,4), "VariableNames", {'exclusion rate'}), ...
                 "VariableNames", {'Left Side', 'Right Side', ...
                                   'Matched Pairs', 'Mismatched Pairs'}, ...
                 "RowNames", varNames);

## Save tables to CSV files
table2csv (U1table, "Ulna-Stas-Athens.csv");
table2csv (U2table, "Ulna-Stas-Crete.csv");
table2csv (U3table, "Ulna-Stas-Granada.csv");
table2csv (U4table, "Ulna-Stas-Pooled.csv");

## Clear unnecessary variables
clear -x F1table F2table F3table F4table H1table H2table H3table H4table ...
         T1table T2table T3table T4table U1table U2table U3table U4table
