## Evaluate all datasets of randomized commingled cases
load ('datasets.mat');


## Evaluate datasets of 30 randomized commingled cases of Femur bones
for i = 1:numel (dataset)
  ## 30 individuals with both elements present (60 bones)
  testcase = dataset(i).F_i30p30;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i30p30 = [stats, Actual];
  else
    F_i30p30 = [F_i30p30; stats, Actual];
  endif

  ## 60 individuals with both elements present (120 bones)
  testcase = dataset(i).F_i60p60;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i60p60 = [stats, Actual];
  else
    F_i60p60 = [F_i60p60; stats, Actual];
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  testcase = dataset(i).F_i60p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i60p50 = [stats, Actual];
  else
    F_i60p50 = [F_i60p50; stats, Actual];
  endif

  ## 100 individuals with both elements present (200 bones)
  testcase = dataset(i).F_i100p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i100p100 = [stats, Actual];
  else
    F_i100p100 = [F_i100p100; stats, Actual];
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  testcase = dataset(i).F_i120p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i120p100 = [stats, Actual];
  else
    F_i120p100 = [F_i120p100; stats, Actual];
  endif

  ## 100 individuals with only one element per individual (100 bones)
  testcase = dataset(i).F_i100p0;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i100p0 = [stats, Actual];
  else
    F_i100p0 = [F_i100p0; stats, Actual];
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  testcase = dataset(i).F_i100p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i100p50 = [stats, Actual];
  else
    F_i100p50 = [F_i100p50; stats, Actual];
  endif

  ## 150 individuals with both elements present (300 bones)
  testcase = dataset(i).F_i150p150;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Femur');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    F_i150p150 = [stats, Actual];
  else
    F_i150p150 = [F_i150p150; stats, Actual];
  endif
endfor

## Create summary table for Femur
varnames = {'i30p30', 'i60p60', 'i60p50', 'i100p100', ...
            'i120p100', 'i100p0', 'i100p50', 'i150p150'};
rownames = {'Reportedly sorted elements', 'Correctly sorted elements', ...
            'TPR of sorted elements', 'Reportedly sorted pairs', ...
            'Correctly sorted pairs', 'TPR of sorted pairs', ...
            'Reportedly sorted singles', 'Correctly sorted singles', ...
            'TPR of sorted singles', 'Reportedly identified individuals', ...
            'Correctly identified individuals', 'TPR of identified individuals', ...
            'Unsorted elements', 'Plausible pairs'};
RSE = mean ([F_i30p30.Sorted, F_i60p60.Sorted,  F_i60p50.Sorted, F_i100p100.Sorted, ...
             F_i120p100.Sorted, F_i100p0.Sorted,  F_i100p50.Sorted, F_i150p150.Sorted]);
ASE = mean ([F_i30p30{:,8}*2+F_i30p30{:,9}, F_i60p60{:,8}*2+F_i60p60{:,9}, ...
             F_i60p50{:,8}*2+F_i60p50{:,9}, F_i100p100{:,8}*2+F_i100p100{:,9}, ...
             F_i120p100{:,8}*2+F_i120p100{:,9}, F_i100p0{:,8}*2+F_i100p0{:,9}, ...
             F_i100p50{:,8}*2+F_i100p50{:,9}, F_i150p150{:,8}*2+F_i150p150{:,9}]);
TPR_SE = ASE ./ RSE;
RSP = mean ([F_i30p30.Paired, F_i30p30.Paired,  F_i60p50.Paired, F_i100p100.Paired, ...
             F_i120p100.Paired, F_i100p0.Paired,  F_i100p50.Paired, F_i150p150.Paired]);
ASP = mean ([F_i30p30{:,8}, F_i30p30{:,8},  F_i60p50{:,8}, F_i100p100{:,8}, ...
             F_i120p100{:,8}, F_i100p0{:,8},  F_i100p50{:,8}, F_i150p150{:,8}]);
TPR_SP = ASP ./ RSP;
RSS = mean ([F_i30p30.Single, F_i30p30.Single,  F_i60p50.Single, F_i100p100.Single, ...
             F_i120p100.Single, F_i100p0.Single,  F_i100p50.Single, F_i150p150.Single]);
ASS = mean ([F_i30p30{:,9}, F_i30p30{:,9},  F_i60p50{:,9}, F_i100p100{:,9}, ...
             F_i120p100{:,9}, F_i100p0{:,9},  F_i100p50{:,9}, F_i150p150{:,9}]);
TPR_SS = ASS ./ RSS;
RII = mean ([F_i30p30.Individuals, F_i30p30.Individuals,  F_i60p50.Individuals, ...
             F_i100p100.Individuals, F_i120p100.Individuals, F_i100p0.Individuals, ...
             F_i100p50.Individuals, F_i150p150.Individuals]);
AII = mean ([F_i30p30{:,9}, F_i30p30{:,9},  F_i60p50{:,9}, F_i100p100{:,9}, ...
             F_i120p100{:,9}, F_i100p0{:,9},  F_i100p50{:,9}, F_i150p150{:,9}]);
TPR_II = AII ./ RII;
UE  = mean ([F_i30p30.Unsorted, F_i60p60.Unsorted,  F_i60p50.Unsorted, ...
             F_i100p100.Unsorted, F_i120p100.Unsorted, F_i100p0.Unsorted, ...
             F_i100p50.Unsorted, F_i150p150.Unsorted]);
PP  = mean ([F_i30p30.Plausible, F_i60p60.Plausible,  F_i60p50.Plausible, ...
             F_i100p100.Plausible, F_i120p100.Plausible, F_i100p0.Plausible, ...
             F_i100p50.Plausible, F_i150p150.Plausible]);
M = [RSE; ASE; TPR_SE; RSP; ASP; TPR_SP; RSS; ASS; TPR_SS; RII; AII; TPR_II; UE; PP];
values = round (M * 100) / 100;
Femur = array2table (values, 'VariableNames',  varnames, 'RowNames', rownames);




## Evaluate datasets of randomized commingled cases of Humerus bones
for i = 1:numel (dataset)
  ## 30 individuals with both elements present (60 bones)
  testcase = dataset(i).H_i30p30;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i30p30 = [stats, Actual];
  else
    H_i30p30 = [H_i30p30; stats, Actual];
  endif

  ## 60 individuals with both elements present (120 bones)
  testcase = dataset(i).H_i60p60;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i60p60 = [stats, Actual];
  else
    H_i60p60 = [H_i60p60; stats, Actual];
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  testcase = dataset(i).H_i60p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i60p50 = [stats, Actual];
  else
    H_i60p50 = [H_i60p50; stats, Actual];
  endif

  ## 100 individuals with both elements present (200 bones)
  testcase = dataset(i).H_i100p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i100p100 = [stats, Actual];
  else
    H_i100p100 = [H_i100p100; stats, Actual];
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  testcase = dataset(i).H_i120p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i120p100 = [stats, Actual];
  else
    H_i120p100 = [H_i120p100; stats, Actual];
  endif

  ## 100 individuals with only one element per individual (100 bones)
  testcase = dataset(i).H_i100p0;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i100p0 = [stats, Actual];
  else
    H_i100p0 = [H_i100p0; stats, Actual];
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  testcase = dataset(i).H_i100p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i100p50 = [stats, Actual];
  else
    H_i100p50 = [H_i100p50; stats, Actual];
  endif

  ## 150 individuals with both elements present (300 bones)
  testcase = dataset(i).H_i150p150;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Humerus');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    H_i150p150 = [stats, Actual];
  else
    H_i150p150 = [H_i150p150; stats, Actual];
  endif
endfor

## Create summary table for Humerus
varnames = {'i30p30', 'i60p60', 'i60p50', 'i100p100', ...
            'i120p100', 'i100p0', 'i100p50', 'i150p150'};
rownames = {'Reportedly sorted elements', 'Correctly sorted elements', ...
            'TPR of sorted elements', 'Reportedly sorted pairs', ...
            'Correctly sorted pairs', 'TPR of sorted pairs', ...
            'Reportedly sorted singles', 'Correctly sorted singles', ...
            'TPR of sorted singles', 'Reportedly identified individuals', ...
            'Correctly identified individuals', 'TPR of identified individuals', ...
            'Unsorted elements', 'Plausible pairs'};
RSE = mean ([H_i30p30.Sorted, H_i60p60.Sorted,  H_i60p50.Sorted, H_i100p100.Sorted, ...
             H_i120p100.Sorted, H_i100p0.Sorted,  H_i100p50.Sorted, H_i150p150.Sorted]);
ASE = mean ([H_i30p30{:,8}*2+H_i30p30{:,9}, H_i60p60{:,8}*2+H_i60p60{:,9}, ...
             H_i60p50{:,8}*2+H_i60p50{:,9}, H_i100p100{:,8}*2+H_i100p100{:,9}, ...
             H_i120p100{:,8}*2+H_i120p100{:,9}, H_i100p0{:,8}*2+H_i100p0{:,9}, ...
             H_i100p50{:,8}*2+H_i100p50{:,9}, H_i150p150{:,8}*2+H_i150p150{:,9}]);
TPR_SE = ASE ./ RSE;
RSP = mean ([H_i30p30.Paired, H_i30p30.Paired,  H_i60p50.Paired, H_i100p100.Paired, ...
             H_i120p100.Paired, H_i100p0.Paired,  H_i100p50.Paired, H_i150p150.Paired]);
ASP = mean ([H_i30p30{:,8}, H_i30p30{:,8},  H_i60p50{:,8}, H_i100p100{:,8}, ...
             H_i120p100{:,8}, H_i100p0{:,8},  H_i100p50{:,8}, H_i150p150{:,8}]);
TPR_SP = ASP ./ RSP;
RSS = mean ([H_i30p30.Single, H_i30p30.Single,  H_i60p50.Single, H_i100p100.Single, ...
             H_i120p100.Single, H_i100p0.Single,  H_i100p50.Single, H_i150p150.Single]);
ASS = mean ([H_i30p30{:,9}, H_i30p30{:,9},  H_i60p50{:,9}, H_i100p100{:,9}, ...
             H_i120p100{:,9}, H_i100p0{:,9},  H_i100p50{:,9}, H_i150p150{:,9}]);
TPR_SS = ASS ./ RSS;
RII = mean ([H_i30p30.Individuals, H_i30p30.Individuals,  H_i60p50.Individuals, ...
             H_i100p100.Individuals, H_i120p100.Individuals, H_i100p0.Individuals, ...
             H_i100p50.Individuals, H_i150p150.Individuals]);
AII = mean ([H_i30p30{:,9}, H_i30p30{:,9},  H_i60p50{:,9}, H_i100p100{:,9}, ...
             H_i120p100{:,9}, H_i100p0{:,9},  H_i100p50{:,9}, H_i150p150{:,9}]);
TPR_II = AII ./ RII;
UE  = mean ([H_i30p30.Unsorted, H_i60p60.Unsorted,  H_i60p50.Unsorted, ...
             H_i100p100.Unsorted, H_i120p100.Unsorted, H_i100p0.Unsorted, ...
             H_i100p50.Unsorted, H_i150p150.Unsorted]);
PP  = mean ([H_i30p30.Plausible, H_i60p60.Plausible,  H_i60p50.Plausible, ...
             H_i100p100.Plausible, H_i120p100.Plausible, H_i100p0.Plausible, ...
             H_i100p50.Plausible, H_i150p150.Plausible]);
M = [RSE; ASE; TPR_SE; RSP; ASP; TPR_SP; RSS; ASS; TPR_SS; RII; AII; TPR_II; UE; PP];
values = round (M * 100) / 100;
Humerus = array2table (values, 'VariableNames',  varnames, 'RowNames', rownames);




## Evaluate datasets of randomized commingled cases of Tibia bones
for i = 1:numel (dataset)
  ## 30 individuals with both elements present (60 bones)
  testcase = dataset(i).T_i30p30;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i30p30 = [stats, Actual];
  else
    T_i30p30 = [T_i30p30; stats, Actual];
  endif

  ## 60 individuals with both elements present (120 bones)
  testcase = dataset(i).T_i60p60;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i60p60 = [stats, Actual];
  else
    T_i60p60 = [T_i60p60; stats, Actual];
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  testcase = dataset(i).T_i60p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i60p50 = [stats, Actual];
  else
    T_i60p50 = [T_i60p50; stats, Actual];
  endif

  ## 100 individuals with both elements present (200 bones)
  testcase = dataset(i).T_i100p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i100p100 = [stats, Actual];
  else
    T_i100p100 = [T_i100p100; stats, Actual];
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  testcase = dataset(i).T_i120p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i120p100 = [stats, Actual];
  else
    T_i120p100 = [T_i120p100; stats, Actual];
  endif

  ## 100 individuals with only one element per individual (100 bones)
  testcase = dataset(i).T_i100p0;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i100p0 = [stats, Actual];
  else
    T_i100p0 = [T_i100p0; stats, Actual];
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  testcase = dataset(i).T_i100p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i100p50 = [stats, Actual];
  else
    T_i100p50 = [T_i100p50; stats, Actual];
  endif

  ## 150 individuals with both elements present (300 bones)
  testcase = dataset(i).T_i150p150;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Tibia');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    T_i150p150 = [stats, Actual];
  else
    T_i150p150 = [T_i150p150; stats, Actual];
  endif
endfor

## Create summary table for Tibia
varnames = {'i30p30', 'i60p60', 'i60p50', 'i100p100', ...
            'i120p100', 'i100p0', 'i100p50', 'i150p150'};
rownames = {'Reportedly sorted elements', 'Correctly sorted elements', ...
            'TPR of sorted elements', 'Reportedly sorted pairs', ...
            'Correctly sorted pairs', 'TPR of sorted pairs', ...
            'Reportedly sorted singles', 'Correctly sorted singles', ...
            'TPR of sorted singles', 'Reportedly identified individuals', ...
            'Correctly identified individuals', 'TPR of identified individuals', ...
            'Unsorted elements', 'Plausible pairs'};
RSE = mean ([T_i30p30.Sorted, T_i60p60.Sorted,  T_i60p50.Sorted, T_i100p100.Sorted, ...
             T_i120p100.Sorted, T_i100p0.Sorted,  T_i100p50.Sorted, T_i150p150.Sorted]);
ASE = mean ([T_i30p30{:,8}*2+T_i30p30{:,9}, T_i60p60{:,8}*2+T_i60p60{:,9}, ...
             T_i60p50{:,8}*2+T_i60p50{:,9}, T_i100p100{:,8}*2+T_i100p100{:,9}, ...
             T_i120p100{:,8}*2+T_i120p100{:,9}, T_i100p0{:,8}*2+T_i100p0{:,9}, ...
             T_i100p50{:,8}*2+T_i100p50{:,9}, T_i150p150{:,8}*2+T_i150p150{:,9}]);
TPR_SE = ASE ./ RSE;
RSP = mean ([T_i30p30.Paired, T_i30p30.Paired,  T_i60p50.Paired, T_i100p100.Paired, ...
             T_i120p100.Paired, T_i100p0.Paired,  T_i100p50.Paired, T_i150p150.Paired]);
ASP = mean ([T_i30p30{:,8}, T_i30p30{:,8},  T_i60p50{:,8}, T_i100p100{:,8}, ...
             T_i120p100{:,8}, T_i100p0{:,8},  T_i100p50{:,8}, T_i150p150{:,8}]);
TPR_SP = ASP ./ RSP;
RSS = mean ([T_i30p30.Single, T_i30p30.Single,  T_i60p50.Single, T_i100p100.Single, ...
             T_i120p100.Single, T_i100p0.Single,  T_i100p50.Single, T_i150p150.Single]);
ASS = mean ([T_i30p30{:,9}, T_i30p30{:,9},  T_i60p50{:,9}, T_i100p100{:,9}, ...
             T_i120p100{:,9}, T_i100p0{:,9},  T_i100p50{:,9}, T_i150p150{:,9}]);
TPR_SS = ASS ./ RSS;
RII = mean ([T_i30p30.Individuals, T_i30p30.Individuals,  T_i60p50.Individuals, ...
             T_i100p100.Individuals, T_i120p100.Individuals, T_i100p0.Individuals, ...
             T_i100p50.Individuals, T_i150p150.Individuals]);
AII = mean ([T_i30p30{:,9}, T_i30p30{:,9},  T_i60p50{:,9}, T_i100p100{:,9}, ...
             T_i120p100{:,9}, T_i100p0{:,9},  T_i100p50{:,9}, T_i150p150{:,9}]);
TPR_II = AII ./ RII;
UE  = mean ([T_i30p30.Unsorted, T_i60p60.Unsorted,  T_i60p50.Unsorted, ...
             T_i100p100.Unsorted, T_i120p100.Unsorted, T_i100p0.Unsorted, ...
             T_i100p50.Unsorted, T_i150p150.Unsorted]);
PP  = mean ([T_i30p30.Plausible, T_i60p60.Plausible,  T_i60p50.Plausible, ...
             T_i100p100.Plausible, T_i120p100.Plausible, T_i100p0.Plausible, ...
             T_i100p50.Plausible, T_i150p150.Plausible]);
M = [RSE; ASE; TPR_SE; RSP; ASP; TPR_SP; RSS; ASS; TPR_SS; RII; AII; TPR_II; UE; PP];
values = round (M * 100) / 100;
Tibia = array2table (values, 'VariableNames',  varnames, 'RowNames', rownames);




## Evaluate datasets of randomized commingled cases of Ulna bones
for i = 1:numel (dataset)
  ## 30 individuals with both elements present (60 bones)
  testcase = dataset(i).U_i30p30;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i30p30 = [stats, Actual];
  else
    U_i30p30 = [U_i30p30; stats, Actual];
  endif

  ## 60 individuals with both elements present (120 bones)
  testcase = dataset(i).U_i60p60;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i60p60 = [stats, Actual];
  else
    U_i60p60 = [U_i60p60; stats, Actual];
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  testcase = dataset(i).U_i60p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i60p50 = [stats, Actual];
  else
    U_i60p50 = [U_i60p50; stats, Actual];
  endif

  ## 100 individuals with both elements present (200 bones)
  testcase = dataset(i).U_i100p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i100p100 = [stats, Actual];
  else
    U_i100p100 = [U_i100p100; stats, Actual];
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  testcase = dataset(i).U_i120p100;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i120p100 = [stats, Actual];
  else
    U_i120p100 = [U_i120p100; stats, Actual];
  endif

  ## 100 individuals with only one element per individual (100 bones)
  testcase = dataset(i).U_i100p0;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i100p0 = [stats, Actual];
  else
    U_i100p0 = [U_i100p0; stats, Actual];
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  testcase = dataset(i).U_i100p50;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  s = isnan (sorted.Score);
  singles = [sorted{s,1}; sorted{s,2}];
  singles = singles(! ismissing (singles));
  fcn = @(x) sum (strncmpi (testcase(:,1), x, 6));
  n = sum (cellfun (fcn, cellstr (singles)) == 1);
  p = sum (strncmp (sorted{! s,1}, sorted{! s,2}, 6));
  Actual = table (p, n, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i100p50 = [stats, Actual];
  else
    U_i100p50 = [U_i100p50; stats, Actual];
  endif

  ## 150 individuals with both elements present (300 bones)
  testcase = dataset(i).U_i150p150;
  [sorted, stats, unsorted] = longbone_Pair (testcase, 'Ulna');
  ## Check single and paired bones
  p = sum (strncmp (sorted{:,1}, sorted{:,2}, 6));
  Actual = table (p, 0, 'VariableNames', {'True Pairs', 'True Singles'});
  if (i == 1)
    U_i150p150 = [stats, Actual];
  else
    U_i150p150 = [U_i150p150; stats, Actual];
  endif
endfor

## Create summary table for Ulna
varnames = {'i30p30', 'i60p60', 'i60p50', 'i100p100', ...
            'i120p100', 'i100p0', 'i100p50', 'i150p150'};
rownames = {'Reportedly sorted elements', 'Correctly sorted elements', ...
            'TPR of sorted elements', 'Reportedly sorted pairs', ...
            'Correctly sorted pairs', 'TPR of sorted pairs', ...
            'Reportedly sorted singles', 'Correctly sorted singles', ...
            'TPR of sorted singles', 'Reportedly identified individuals', ...
            'Correctly identified individuals', 'TPR of identified individuals', ...
            'Unsorted elements', 'Plausible pairs'};
RSE = mean ([U_i30p30.Sorted, U_i60p60.Sorted,  U_i60p50.Sorted, U_i100p100.Sorted, ...
             U_i120p100.Sorted, U_i100p0.Sorted,  U_i100p50.Sorted, U_i150p150.Sorted]);
ASE = mean ([U_i30p30{:,8}*2+U_i30p30{:,9}, U_i60p60{:,8}*2+U_i60p60{:,9}, ...
             U_i60p50{:,8}*2+U_i60p50{:,9}, U_i100p100{:,8}*2+U_i100p100{:,9}, ...
             U_i120p100{:,8}*2+U_i120p100{:,9}, U_i100p0{:,8}*2+U_i100p0{:,9}, ...
             U_i100p50{:,8}*2+U_i100p50{:,9}, U_i150p150{:,8}*2+U_i150p150{:,9}]);
TPR_SE = ASE ./ RSE;
RSP = mean ([U_i30p30.Paired, U_i30p30.Paired,  U_i60p50.Paired, U_i100p100.Paired, ...
             U_i120p100.Paired, U_i100p0.Paired,  U_i100p50.Paired, U_i150p150.Paired]);
ASP = mean ([U_i30p30{:,8}, U_i30p30{:,8},  U_i60p50{:,8}, U_i100p100{:,8}, ...
             U_i120p100{:,8}, U_i100p0{:,8},  U_i100p50{:,8}, U_i150p150{:,8}]);
TPR_SP = ASP ./ RSP;
RSS = mean ([U_i30p30.Single, U_i30p30.Single,  U_i60p50.Single, U_i100p100.Single, ...
             U_i120p100.Single, U_i100p0.Single,  U_i100p50.Single, U_i150p150.Single]);
ASS = mean ([U_i30p30{:,9}, U_i30p30{:,9},  U_i60p50{:,9}, U_i100p100{:,9}, ...
             U_i120p100{:,9}, U_i100p0{:,9},  U_i100p50{:,9}, U_i150p150{:,9}]);
TPR_SS = ASS ./ RSS;
RII = mean ([U_i30p30.Individuals, U_i30p30.Individuals,  U_i60p50.Individuals, ...
             U_i100p100.Individuals, U_i120p100.Individuals, U_i100p0.Individuals, ...
             U_i100p50.Individuals, U_i150p150.Individuals]);
AII = mean ([U_i30p30{:,9}, U_i30p30{:,9},  U_i60p50{:,9}, U_i100p100{:,9}, ...
             U_i120p100{:,9}, U_i100p0{:,9},  U_i100p50{:,9}, U_i150p150{:,9}]);
TPR_II = AII ./ RII;
UE  = mean ([U_i30p30.Unsorted, U_i60p60.Unsorted,  U_i60p50.Unsorted, ...
             U_i100p100.Unsorted, U_i120p100.Unsorted, U_i100p0.Unsorted, ...
             U_i100p50.Unsorted, U_i150p150.Unsorted]);
PP  = mean ([U_i30p30.Plausible, U_i60p60.Plausible,  U_i60p50.Plausible, ...
             U_i100p100.Plausible, U_i120p100.Plausible, U_i100p0.Plausible, ...
             U_i100p50.Plausible, U_i150p150.Plausible]);
M = [RSE; ASE; TPR_SE; RSP; ASP; TPR_SP; RSS; ASS; TPR_SS; RII; AII; TPR_II; UE; PP];
values = round (M * 100) / 100;
Ulna = array2table (values, 'VariableNames',  varnames, 'RowNames', rownames);
