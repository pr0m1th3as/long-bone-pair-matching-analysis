## Generate randomized datasets for pair-matching analysis
nsets = 50;

## Load CSG dataset
pkg load csg-dataset
CSG = load ('csg-properties.mat').csg_properties;

## Gather Femur samples from all collections into left and right sides
Ldata = {};
Rdata = {};
Lname = {};
Rname = {};
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 1};
  sIdx = data([2:end],1);
  Lname = [Lname; sIdx];
  side = num2cell (ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Ldata = [Ldata; sIdx, side, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 2};
  sIdx = data([2:end],1);
  Rname = [Rname; sIdx];
  side = num2cell (2 * ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Rdata = [Rdata; sIdx, side, measurements];
endfor
Femur = [Ldata; Rdata];

## Get paired and single bone lists
Lname = cellfun (@(x) x(1:6), Lname, "UniformOutput", false);
Rname = cellfun (@(x) x(1:6), Rname, "UniformOutput", false);
L_TF = ismember (Lname, Rname);
paired_samples = sum (L_TF);

## Create datasets of randomized commingled cases of Femur bones
for i = 1:nsets
  ## 30 individuals with both elements present (60 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 30));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).F_i30p30 = Femur(idx,:);
  ## Check parity
  s = dataset(i).F_i30p30(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 30)
    error ("F_i30p30 set %d does not have correct parity.", i);
  elseif (rows (s) != 60)
    error ("F_i30p30 set %d does not have correct samples.", i);
  endif

  ## 60 individuals with both elements present (120 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 60));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).F_i60p60 = Femur(idx,:);
  ## Check parity
  s = dataset(i).F_i60p60(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 60)
    error ("F_i60p60 set %d does not have correct parity.", i);
  elseif (rows (s) != 120)
    error ("F_i60p60 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 5));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 5));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).F_i60p50 = Femur([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).F_i60p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("F_i60p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 110)
    error ("F_i60p50 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements present (200 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
  dataset(i).F_i100p100 = Femur(idx,:);
  ## Check parity
  s = dataset(i).F_i100p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("F_i100p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 200)
    error ("F_i100p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 10));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 10));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).F_i120p100 = Femur([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).F_i120p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("F_i120p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 220)
    error ("F_i120p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with only one element per individual (100 bones)
  single_list_L = Lname(randperm (numel (Lname), 50));
  TF = ismember (Rname, single_list_L);
  single_list_R = Rname(! TF)(randperm (sum (! TF), 50));
  ids = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).F_i100p0 = Femur(ids,:);
  ## Check parity
  s = dataset(i).F_i100p0(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 0)
    error ("F_i100p0 set %d does not have correct parity.", i);
  elseif (rows (s) != 100)
    error ("F_i100p0 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 25));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 25));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).F_i100p50 = Femur([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).F_i100p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("F_i100p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 150)
    error ("F_i100p50 set %d does not have correct samples.", i);
  endif

  ## 150 individuals with both elements present (300 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 150));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).F_i150p150 = Femur(idx,:);
  ## Check parity
  s = dataset(i).F_i150p150(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 150)
    error ("F_i150p150 set %d does not have correct parity.", i);
  elseif (rows (s) != 300)
    error ("F_i150p150 set %d does not have correct samples.", i);
  endif
endfor


## Gather Humerus samples from all collections into left and right sides
Ldata = {};
Rdata = {};
Lname = {};
Rname = {};
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 3};
  sIdx = data([2:end],1);
  Lname = [Lname; sIdx];
  side = num2cell (ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Ldata = [Ldata; sIdx, side, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 4};
  sIdx = data([2:end],1);
  Rname = [Rname; sIdx];
  side = num2cell (2 * ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Rdata = [Rdata; sIdx, side, measurements];
endfor
Humerus = [Ldata; Rdata];

## Get paired and single bone lists
Lname = cellfun (@(x) x(1:6), Lname, "UniformOutput", false);
Rname = cellfun (@(x) x(1:6), Rname, "UniformOutput", false);
L_TF = ismember (Lname, Rname);
paired_samples = sum (L_TF);

## Create datasets of randomized commingled cases of Humerus bones
for i = 1:nsets
  ## 30 individuals with both elements present (60 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 30));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).H_i30p30 = Humerus(idx,:);
  ## Check parity
  s = dataset(i).H_i30p30(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 30)
    error ("H_i30p30 set %d does not have correct parity.", i);
  elseif (rows (s) != 60)
    error ("H_i30p30 set %d does not have correct samples.", i);
  endif

  ## 60 individuals with both elements present (120 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 60));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).H_i60p60 = Humerus(idx,:);
  ## Check parity
  s = dataset(i).H_i60p60(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 60)
    error ("H_i60p60 set %d does not have correct parity.", i);
  elseif (rows (s) != 120)
    error ("H_i60p60 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 5));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 5));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).H_i60p50 = Humerus([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).H_i60p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("H_i60p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 110)
    error ("H_i60p50 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements present (200 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
  dataset(i).H_i100p100 = Humerus(idx,:);
  ## Check parity
  s = dataset(i).H_i100p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("H_i100p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 200)
    error ("H_i100p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 10));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 10));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).H_i120p100 = Humerus([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).H_i120p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("H_i120p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 220)
    error ("H_i120p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with only one element per individual (100 bones)
  single_list_L = Lname(randperm (numel (Lname), 50));
  TF = ismember (Rname, single_list_L);
  single_list_R = Rname(! TF)(randperm (sum (! TF), 50));
  ids = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).H_i100p0 = Humerus(ids,:);
  ## Check parity
  s = dataset(i).H_i100p0(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 0)
    error ("H_i100p0 set %d does not have correct parity.", i);
  elseif (rows (s) != 100)
    error ("H_i100p0 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 25));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 25));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).H_i100p50 = Humerus([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).H_i100p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("H_i100p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 150)
    error ("H_i100p50 set %d does not have correct samples.", i);
  endif

  ## 150 individuals with both elements present (300 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 150));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).H_i150p150 = Humerus(idx,:);
  ## Check parity
  s = dataset(i).H_i150p150(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 150)
    error ("H_i150p150 set %d does not have correct parity.", i);
  elseif (rows (s) != 300)
    error ("H_i150p150 set %d does not have correct samples.", i);
  endif
endfor


## Gather Tibia samples from all collections into left and right sides
Ldata = {};
Rdata = {};
Lname = {};
Rname = {};
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 5};
  sIdx = data([2:end],1);
  Lname = [Lname; sIdx];
  side = num2cell (ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Ldata = [Ldata; sIdx, side, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 6};
  sIdx = data([2:end],1);
  Rname = [Rname; sIdx];
  side = num2cell (2 * ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Rdata = [Rdata; sIdx, side, measurements];
endfor
Tibia = [Ldata; Rdata];

## Get paired and single bone lists
Lname = cellfun (@(x) x(1:6), Lname, "UniformOutput", false);
Rname = cellfun (@(x) x(1:6), Rname, "UniformOutput", false);
L_TF = ismember (Lname, Rname);
paired_samples = sum (L_TF);

## Create datasets of randomized commingled cases of Tibia bones
for i = 1:nsets
  ## 30 individuals with both elements present (60 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 30));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).T_i30p30 = Tibia(idx,:);
  ## Check parity
  s = dataset(i).T_i30p30(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 30)
    error ("T_i30p30 set %d does not have correct parity.", i);
  elseif (rows (s) != 60)
    error ("T_i30p30 set %d does not have correct samples.", i);
  endif

  ## 60 individuals with both elements present (120 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 60));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).T_i60p60 = Tibia(idx,:);
  ## Check parity
  s = dataset(i).T_i60p60(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 60)
    error ("T_i60p60 set %d does not have correct parity.", i);
  elseif (rows (s) != 120)
    error ("T_i60p60 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 5));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 5));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).T_i60p50 = Tibia([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).T_i60p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("T_i60p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 110)
    error ("T_i60p50 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements present (200 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
  dataset(i).T_i100p100 = Tibia(idx,:);
  ## Check parity
  s = dataset(i).T_i100p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("T_i100p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 200)
    error ("T_i100p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 10));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 10));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).T_i120p100 = Tibia([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).T_i120p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("T_i120p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 220)
    error ("T_i120p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with only one element per individual (100 bones)
  single_list_L = Lname(randperm (numel (Lname), 50));
  TF = ismember (Rname, single_list_L);
  single_list_R = Rname(! TF)(randperm (sum (! TF), 50));
  ids = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).T_i100p0 = Tibia(ids,:);
  ## Check parity
  s = dataset(i).T_i100p0(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 0)
    error ("T_i100p0 set %d does not have correct parity.", i);
  elseif (rows (s) != 100)
    error ("T_i100p0 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 25));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 25));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).T_i100p50 = Tibia([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).T_i100p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("T_i100p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 150)
    error ("T_i100p50 set %d does not have correct samples.", i);
  endif

  ## 150 individuals with both elements present (300 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 150));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).T_i150p150 = Tibia(idx,:);
  ## Check parity
  s = dataset(i).T_i150p150(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 150)
    error ("T_i150p150 set %d does not have correct parity.", i);
  elseif (rows (s) != 300)
    error ("T_i150p150 set %d does not have correct samples.", i);
  endif
endfor


## Gather Ulna samples from all collections into left and right sides
Ldata = {};
Rdata = {};
Lname = {};
Rname = {};
for c = 1:3
  ## Left side
  data = CSG{2, (c - 1) * 8 + 7};
  sIdx = data([2:end],1);
  Lname = [Lname; sIdx];
  side = num2cell (ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Ldata = [Ldata; sIdx, side, measurements];
  ## Right side
  data = CSG{2, (c - 1) * 8 + 8};
  sIdx = data([2:end],1);
  Rname = [Rname; sIdx];
  side = num2cell (2 * ones (size (sIdx)));
  measurements = data([2:end], [4:64]);
  Rdata = [Rdata; sIdx, side, measurements];
endfor
Ulna = [Ldata; Rdata];

## Get paired and single bone lists
Lname = cellfun (@(x) x(1:6), Lname, "UniformOutput", false);
Rname = cellfun (@(x) x(1:6), Rname, "UniformOutput", false);
L_TF = ismember (Lname, Rname);
paired_samples = sum (L_TF);

## Create datasets of randomized commingled cases of Ulna bones
for i = 1:nsets
  ## 30 individuals with both elements present (60 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 30));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).U_i30p30 = Ulna(idx,:);
  ## Check parity
  s = dataset(i).U_i30p30(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 30)
    error ("U_i30p30 set %d does not have correct parity.", i);
  elseif (rows (s) != 60)
    error ("U_i30p30 set %d does not have correct samples.", i);
  endif

  ## 60 individuals with both elements present (120 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 60));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).U_i60p60 = Ulna(idx,:);
  ## Check parity
  s = dataset(i).U_i60p60(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 60)
    error ("U_i60p60 set %d does not have correct parity.", i);
  elseif (rows (s) != 120)
    error ("U_i60p60 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 10 with single elements (110 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 5));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 5));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).U_i60p50 = Ulna([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).U_i60p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("U_i60p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 110)
    error ("U_i60p50 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements present (200 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
  dataset(i).U_i100p100 = Ulna(idx,:);
  ## Check parity
  s = dataset(i).U_i100p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("U_i100p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 200)
    error ("U_i100p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with both elements and 20 with single elements (220 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 100));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 10));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 10));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).U_i120p100 = Ulna([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).U_i120p100(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 100)
    error ("U_i120p100 set %d does not have correct parity.", i);
  elseif (rows (s) != 220)
    error ("U_i120p100 set %d does not have correct samples.", i);
  endif

  ## 100 individuals with only one element per individual (100 bones)
  single_list_L = Lname(randperm (numel (Lname), 50));
  TF = ismember (Rname, single_list_L);
  single_list_R = Rname(! TF)(randperm (sum (! TF), 50));
  ids = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).U_i100p0 = Ulna(ids,:);
  ## Check parity
  s = dataset(i).U_i100p0(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 0)
    error ("U_i100p0 set %d does not have correct parity.", i);
  elseif (rows (s) != 100)
    error ("U_i100p0 set %d does not have correct samples.", i);
  endif

  ## 50 individuals with both elements and 50 with single elements (150 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 50));
  LpTF = ismember (Lname, paired_list);
  RpTF = ismember (Rname, paired_list);
  tmpL = Lname(! LpTF);
  single_list_L = tmpL(randperm (numel (tmpL), 25));
  RsTF = ismember (Rname, single_list_L);
  tmpR = Rname(! (RpTF | RsTF));
  single_list_R = tmpR(randperm (numel (tmpR), 25));
  sidx = [ismember(Lname, single_list_L); ismember(Rname, single_list_R)];
	dataset(i).U_i100p50 = Ulna([LpTF; RpTF] | sidx,:);
  ## Check parity
  s = dataset(i).U_i100p50(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 50)
    error ("U_i100p50 set %d does not have correct parity.", i);
  elseif (rows (s) != 150)
    error ("U_i100p50 set %d does not have correct samples.", i);
  endif

  ## 150 individuals with both elements present (300 bones)
  paired_list = Lname(L_TF)(randperm (paired_samples, 150));
  idx = [ismember(Lname, paired_list); ismember(Rname, paired_list)];
	dataset(i).U_i150p150 = Ulna(idx,:);
  ## Check parity
  s = dataset(i).U_i150p150(:,1:2);
  p = ismember (cellfun (@(x) x(1:6), s(cellfun (@(x) x == 1, s(:,2)), 1), ...
                         "UniformOutput", false), ...
                cellfun (@(x) x(1:6), s(cellfun (@(x) x == 2, s(:,2)), 1), ...
                         "UniformOutput", false));
  if (sum (p) != 150)
    error ("U_i150p150 set %d does not have correct parity.", i);
  elseif (rows (s) != 300)
    error ("U_i150p150 set %d does not have correct samples.", i);
  endif
endfor

## Save dataset into file
save ('-binary', 'datasets.mat', 'dataset')
clear -x dataset
