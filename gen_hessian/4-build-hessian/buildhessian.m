% Build the Hessian for an arbitrary structure from the RNL

%%% Read structure and rnl
struc = dlmread('SNEW');
rnl = dlmread('RNL');
nl = dlmread('NLNEW');
fh_ijlist = fopen('IJLIST', 'w');
fh_iilist = fopen('IILIST', 'w');
[N blah] = size(struc);

fprintf(fh_iilist, 'ii\n');
fprintf(fh_ijlist, 'ij\n');
%%% Read LAT file
fh_lat = fopen('LAT', 'r');
l = fgetl(fh_lat);
A = [];
for a = 1:3
  l = fgetl(fh_lat);
  A = [A; str2num(l)];
end
l = fgetl(fh_lat);
nb = str2num(l(1));
basisd = [];
for a = 1:nb
  l = fgetl(fh_lat); 
  basisd = [basisd; str2num(l)];
end
basisc = basisd*A;

strucd = struc*inv(A);

% Loop through all self interactions and write FCs
% Self interactions correspond to entries in the RNL where bj = 0
bjindx = find(rnl(:,5)==0);
iifcs = rnl(bjindx, :);
%for a = 1:N
%  % Determine basis atom type of ii
%  id = strucd(a,:);
%  for b = 1:nb
%    test = id - basisd(b,:);
%    if round(test.*10^4)./(10^4) - round(test) == 0
%      typei = b;
%    end
%  end
%  % Find where typei matches bi in iifcs
%  matchindx = find(iifcs(:,1)==typei);
%  fcs = iifcs(matchindx,:);
%  fcs = fcs(6:end);
%  % Write these fcs to the Hessian file
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 1 a 1 fcs(1)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 1 a 2 fcs(2)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 1 a 3 fcs(3)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 2 a 1 fcs(4)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 2 a 2 fcs(5)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 2 a 3 fcs(6)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 3 a 1 fcs(7)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 3 a 2 fcs(8)]);
%  fprintf(fh_iilist, '%i %i %i %i %.15f\n', [a 3 a 3 fcs(9)]);
%end

% Loop through all pairs and write the FCs
[l w] = size(nl);
for a = 1:l
  pair = nl(a,:);
  itag = pair(1);
  jtag = pair(2);
  if itag == jtag
    % Determine basis atom type of ii
    id = strucd(itag,:);
    for b = 1:nb
      test = id - basisd(b,:);
      if round(test.*10^4)./(10^4) - round(test) == 0
        typei = b;
      end
    end
    % Find where typei matches bi in iifcs
    matchindx = find(iifcs(:,1)==typei);
    fcs = iifcs(matchindx,:);
    fcs = fcs(6:end);
    % Write these fcs to the Hessian file
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 1 itag 1 fcs(1)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 1 itag 2 fcs(2)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 1 itag 3 fcs(3)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 2 itag 1 fcs(4)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 2 itag 2 fcs(5)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 2 itag 3 fcs(6)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 3 itag 1 fcs(7)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 3 itag 2 fcs(8)]);
    fprintf(fh_iilist, '%i %i %i %i %.15f\n', [itag 3 itag 3 fcs(9)]);
  else
    disp = pair(3:5);
    dispd = disp*inv(A);
    % Determine basis atom type of i and j
    id = strucd(itag,:);
    jd = strucd(jtag,:);
    for b = 1:nb
      itest = id - basisd(b,:);
      jtest = jd - basisd(b,:);
      % Check if whole numbers
      if round(itest.*1e4)/1e4 == round(itest)
        bi = b;
      end
      if round(jtest.*1e4)/1e4 == round(jtest)
        bj = b;
      end
    end
    % Get n values for displacement
    n = dispd + basisd(bi, :) - basisd(bj,:);
    % Form a [bi n1 n2 n3 bj] vector for comparison to RNL
    testvec = [bi n bj];
    testvec = round(testvec.*1e8)./1e8;
    % Find where this combo exists in the RNL
    rnleqtestvec = rnl(:,1:5)==testvec;
    rnleqtestvect = rnleqtestvec';
    allrnleqtestvect = all(rnleqtestvect);
    rnlpair = find(allrnleqtestvect);
    fcpair = rnl(rnlpair,:);
    fcs = fcpair(6:end);
    % Write these fcs to the Hessian file
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 1 jtag 1 fcs(1)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 1 jtag 2 fcs(2)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 1 jtag 3 fcs(3)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 2 jtag 1 fcs(4)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 2 jtag 2 fcs(5)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 2 jtag 3 fcs(6)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 3 jtag 1 fcs(7)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 3 jtag 2 fcs(8)]);
    fprintf(fh_ijlist, '%i %i %i %i %.15f\n', [itag 3 jtag 3 fcs(9)]);
  end
end
fclose(fh_ijlist);
fclose(fh_iilist);
