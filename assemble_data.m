## Process all measurement data and create matched-mismatched datasets
## for each bone across all collections

## Load CSG dataset
pkg load csg-dataset
CSG = load ('csg-properties.mat').csg_properties;
nametoid = @(x) str2num (x([4:6]));

## Gather Femur samples from all collections into left and right sides
Femur = [];
for c = 1:3
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
endfor

## Split dataset into left and right side
Femur_L = Femur(Femur(:,2) == 1, [6:66]);
Femur_R = Femur(Femur(:,2) == 2, [6:66]);
F_idx_L = Femur(Femur(:,2) == 1, [1,3:5]);
F_idx_R = Femur(Femur(:,2) == 2, [1,3:5]);

## Get number of available femur bones in each side
nsamples_L = size (Femur_L, 1);
jj = 1:size (Femur_R, 1);
jj = jj';

## Compute absolute values
Lvec = abs (Femur_L);
Rvec = abs (Femur_R);

## Scan through all possible ulna pairs and assemble matched/mismatched diffs
## betweeen left and right sides
F_match = [];
F_mismatch = [];
for ii = 1:nsamples_L
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (F_idx_L(ii,1:3) == F_idx_R(:,1:3), 2);
  if (any (idx))
    F_match = [F_match; ii, jj(idx), Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  F_mismatch = [F_mismatch; repmat(ii, sum (idx), 1), jj(idx), Dvec(idx,:)];
endfor

## Save to file
save ('-binary', 'Femur-Data.mat', 'F_match', 'F_mismatch', 'F_idx_L', 'F_idx_R')
clear -x CSG nametoid


## Gather Humerus samples from all collections into left and right sides
Humerus = [];
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 3};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Humerus = [Humerus; sIdx, side, collection, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 4};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Humerus = [Humerus; sIdx, side, collection, measurements];
endfor

## Split dataset into left and right side
Humerus_L = Humerus(Humerus(:,2) == 1, [6:66]);
Humerus_R = Humerus(Humerus(:,2) == 2, [6:66]);
H_idx_L = Humerus(Humerus(:,2) == 1, [1,3:5]);
H_idx_R = Humerus(Humerus(:,2) == 2, [1,3:5]);

## Get number of available humerus bones in each side
nsamples_L = size (Humerus_L, 1);
jj = 1:size (Humerus_R, 1);
jj = jj';

## Compute absolute values
Lvec = abs (Humerus_L);
Rvec = abs (Humerus_R);

## Scan through all possible ulna pairs and assemble matched/mismatched diffs
## betweeen left and right sides
H_match = [];
H_mismatch = [];
for ii = 1:nsamples_L
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (H_idx_L(ii,1:3) == H_idx_R(:,1:3), 2);
  if (any (idx))
    H_match = [H_match; ii, jj(idx), Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  H_mismatch = [H_mismatch; repmat(ii, sum (idx), 1), jj(idx), Dvec(idx,:)];
endfor

## Save to file
save ('-binary', 'Humerus-Data.mat', 'H_match', 'H_mismatch', 'H_idx_L', 'H_idx_R')
clear -x CSG nametoid


## Gather Tibia samples from all collections into left and right sides
Tibia = [];
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 5};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Tibia = [Tibia; sIdx, side, collection, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 6};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Tibia = [Tibia; sIdx, side, collection, measurements];
endfor

## Split dataset into left and right side
Tibia_L = Tibia(Tibia(:,2) == 1, [6:66]);
Tibia_R = Tibia(Tibia(:,2) == 2, [6:66]);
T_idx_L = Tibia(Tibia(:,2) == 1, [1,3:5]);
T_idx_R = Tibia(Tibia(:,2) == 2, [1,3:5]);

## Get number of available humerus bones in each side
nsamples_L = size (Tibia_L, 1);
jj = 1:size (Tibia_R, 1);
jj = jj';

## Compute absolute values
Lvec = abs (Tibia_L);
Rvec = abs (Tibia_R);

## Scan through all possible ulna pairs and assemble matched/mismatched diffs
## betweeen left and right sides
T_match = [];
T_mismatch = [];
for ii = 1:nsamples_L
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (T_idx_L(ii,1:3) == T_idx_R(:,1:3), 2);
  if (any (idx))
    T_match = [T_match; ii, jj(idx), Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  T_mismatch = [T_mismatch; repmat(ii, sum (idx), 1), jj(idx), Dvec(idx,:)];
endfor

## Save to file
save ('-binary', 'Tibia-Data.mat', 'T_match', 'T_mismatch', 'T_idx_L', 'T_idx_R')
clear -x CSG nametoid


## Gather Ulna samples from all collections into left and right sides
Ulna = [];
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 7};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Ulna = [Ulna; sIdx, side, collection, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 8};
  sIdx = cell2mat (cellfun (nametoid, data([2:end],1), "UniformOutput", false));
  side = 2 * ones (size (sIdx));
  collection = c * ones (size (sIdx));
  measurements = cell2mat (data([2:end], [2:64]));
  Ulna = [Ulna; sIdx, side, collection, measurements];
endfor

## Split dataset into left and right side
Ulna_L = Ulna(Ulna(:,2) == 1, [6:66]);
Ulna_R = Ulna(Ulna(:,2) == 2, [6:66]);
U_idx_L = Ulna(Ulna(:,2) == 1, [1,3:5]);
U_idx_R = Ulna(Ulna(:,2) == 2, [1,3:5]);

## Get number of available humerus bones in each side
nsamples_L = size (Ulna_L, 1);
jj = 1:size (Ulna_R, 1);
jj = jj';

## Compute absolute values
Lvec = abs (Ulna_L);
Rvec = abs (Ulna_R);

## Scan through all possible ulna pairs and assemble matched/mismatched diffs
## betweeen left and right sides
U_match = [];
U_mismatch = [];
for ii = 1:nsamples_L
  ## Get absolute differences between left and right side
  Dvec = Lvec(ii,:) - Rvec;
  ## Get index of matching pair
  idx = all (U_idx_L(ii,1:3) == U_idx_R(:,1:3), 2);
  if (any (idx))
    U_match = [U_match; ii, jj(idx), Dvec(idx,:)];
  endif
  ## Get indices of mismatched pairs
  idx = ! idx;
  U_mismatch = [U_mismatch; repmat(ii, sum (idx), 1), jj(idx), Dvec(idx,:)];
endfor

## Save to file
save ('-binary', 'Ulna-Data.mat', 'U_match', 'U_mismatch', 'U_idx_L', 'U_idx_R')
clear
