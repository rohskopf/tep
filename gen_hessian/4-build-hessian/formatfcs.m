%%% Format .fcs file for my purpse

fh_debug = fopen('DEBUG', 'w');
%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%
fh_fcs = fopen('si222.fcs', 'r');
nl = dlmread('NLORIG');
%N = 70; % Number of nearest neighbors
struc = dlmread('SORIG'); % Structure in cartesian coordinates

% Lattice vector matrix
A = 5.472.*[0.0 0.5 0.5;
            0.5 0.0 0.5;
            0.5 0.5 0.0];

% Basis vector matrix
basisd = [0.0 0.0 0.0;
          0.25 0.25 0.25];

% IDs of basis atoms
batom = [1 33]; % You must know the atom ids that the basis consists of (unfortunately)

% Minimum cutoff
mincut = 3; % Force constants below this distance will be set to zero


%%%%%%%%%%%%% CODE (don't edit) %%%%%%%%%%%%%%%
[nb blah] = size(basisd);

% Declare matrix to store rnl in
%rnl = zeros(nb*(N+1), 14);
rnl=[];

a1 = A(1,:);
a2 = A(2,:);
a3 = A(3,:);
strucd = struc*inv(A);

% Build basis 
%basisd = [];
%for a = 1:length(batoms)
%  basisd = [basisd; strucd(batoms(a),:)];
%end

% Write LAT file
fh_lat = fopen('LAT','w');
fprintf(fh_lat,'LATTICE VECTORS\n');
for a = 1:3
  fprintf(fh_lat, '%.8f %.8f %.8f\n', A(a,:));
end
fprintf(fh_lat,'%i BASIS VECTORS\n', nb);
for a = 1:nb
  fprintf(fh_lat, '%.8f %.8f %.8f\n', basisd(a,:));
end
fclose(fh_lat);


l = fgetl(fh_fcs);
while (~strcmp(l,' ------------------------ All FCs below ------------------------'))
  l = fgetl(fh_fcs);

end
l = fgetl(fh_fcs);

counter = 0;
totalfcs = 0;
fclist = [];
while(ischar(l))
%while(~strcmp(l,' **FC3')) % USE THIS LOOP WHEN THERE ARE 3RD ORDER IFCS IN FILE
  l = fgetl(fh_fcs);
  if any(l == '#') % Then we have a FC header line
    blah = l;
    % Find the 4th space indx
    spaceindx = find(l == ' ');
    space = spaceindx(4);
    goodstuff = l(space:end);
    goodstuff = str2num(goodstuff);
    num = goodstuff(1);
    totalfcs = totalfcs + num;
    fc = goodstuff(2);
    %if fc < 0
    %  fc = -1*fc;
    %end
    % Now loop through these FCs
    for a = 1:num
      l = fgetl(fh_fcs);
      % Find indx of the period
      periodindx = find(l == '.');
      % Find first space after period (get 2nd half of line)
      lcut = l(periodindx:end);
      spaceindx = find(lcut == ' ');
      spaceindx = spaceindx(1);
      l2nd = lcut(spaceindx:end);
      % Get first half of line
      % The rest of the multiple should be 5 chars after the period
      l1st = l(1:periodindx+6);

      % Alright... Now we have the multiple (in l1st) and the IFC component (in l2nd)
      l1stnum = str2num(l1st);
      mult = l1stnum(2);

      % How can we interpret the IFC component (i alpha, j beta)?
      % The two ints will be i and j
      % The two chars will be alpha and beta
      % Start with omitting all x,y,z to get only i and j
      xyzindx = [strfind(l2nd, 'x') strfind(l2nd, 'y') strfind(l2nd, 'z')];
      ij = l2nd;
      ij(xyzindx) = [];

      % Now omit all 0,1,2,3,4,... to get only alpha and beta 
      numindx = [strfind(l2nd, '0') strfind(l2nd, '1') strfind(l2nd, '2') strfind(l2nd, '3') strfind(l2nd, '4') strfind(l2nd, '5') strfind(l2nd, '6') strfind(l2nd, '7') strfind(l2nd, '8')  strfind(l2nd, '9')];
      xyz = l2nd;
      xyz(numindx) = [];


      % Now extract i and j from struc
      ij = str2num(ij);
      ii = struc(ij(1),:)*inv(A);
      jj = struc(ij(2),:)*inv(A);

      % Determine basis type of ii
      for b = 1:nb
        test = ii - basisd(b,:);
        if round(test.*10^4)./(10^4) - round(test) == 0
          typeii = b;
        end
      end

      % Determine basis type of jj
      for b = 1:nb
        test = jj - basisd(b,:);
        if round(test.*10^4)./(10^4) - round(test) == 0
          typejj = b;
        end
      end

      % Now get the nvector
      % Need to look up the displacement in the neighbor list
      if ij(1) == ij(2)
        n = [0 0 0];
        if mincut > 0
          fc = 0; %%% i-i interactions: Set fc=fc if unchanged. Set fc=0.0 if you want 0.
        end
      else
        itagindx = find(nl(:,1)==ij(1));
        iinl = nl(itagindx,:);
        jtagindx = find(iinl(:,2) == ij(2));
        pairs = iinl(jtagindx,:);
        disp = pairs(:,3:5);
        if norm(disp) < mincut
          fc = 0; %%% i-j interactions: Set fc=fc if unchanged. Set fc=0.0 if you want 0.
        end
        dispd = disp*inv(A);
        n = dispd + basisd(typeii, :) - basisd(typejj, :);
        [n_length n_width] = size(n);
        %for t = 1:n_length
        %  ntemp = n(t,:);
        %  fprintf(fh_debug, '%i %i %i %i %i %i %i \n', ij(1), ij(2), typeii, round(ntemp), typejj);
        %end
      end

      % Need FC component
      spaceindices = strfind(xyz,' ');
      xyz(spaceindices)=[]; % Now xyz has 2 components (no spaces)
      for c = 1:9
        if xyz =='xx'
          fcab = 1;
        elseif xyz =='xy'
          fcab = 2;
        elseif xyz =='xz'
          fcab = 3;
        elseif xyz =='yx'
          fcab = 4;
        elseif xyz =='yy'
          fcab = 5;
        elseif xyz =='yz'
          fcab = 6;
        elseif xyz =='zx'
          fcab = 7;
        elseif xyz =='zy'
          fcab = 8;
        elseif xyz =='zz'
          fcab = 9;
        end
      end

      % Now store i,j, ii type, n, jj type, FC, and FC component
      [n_length n_width] = size(n);
      for t = 1:n_length
        ntemp = n(t,:);
        fclist = [fclist; ij(1) ij(2) typeii ntemp typejj fc*mult fcab]; 
      end

      
    end
  end
  counter = counter + 1;
end
fclose(fh_fcs);


% Now we have the fclist, just build the rnl
% Loop over all basis atoms
rnlindx = 0;
for b = 1:nb
  ii = batom(b);
  % Do the self interaction first
  fcindx = find(fclist(:,1)==ii);
  fcnl=fclist(fcindx,:);
  fcpairindx = find(fcnl(:,2)==ii);
  fcpair = fcnl(fcpairindx,:);
  nvals = fcpair(1,3:7);
  fcs = fcpair(:,end-1:end);
  [numfcs blah] = size(fcs);
  % Put this pair into the rnl matrix
  % This is a self interaction so set basis atom 2 to be 0
  nvals(5) = 0;
  rnl(1+rnlindx,1:5) = nvals;
  for f = 1:numfcs
    % n values
   rnl(1+rnlindx,fcs(f,2)+5) = fcs(f,1);
  end

  % Now do the ij interactions
  nlindx = find(nl(:,1)==ii);
  iinl = nl(nlindx,:);
  %dlmwrite(fh_debug, iinl, ' ');
  [N blah] = size(iinl); % Get number of neighbors N
  fcindx = find(fclist(:,1)==ii);
  fcnl = fclist(fcindx,:);
  for p = 1:N % Loop over all pair (neighbors)
    fprintf(fh_debug, 'NEW PAIR---------------\n');
    pair = iinl(p,:); %xtract a single pair
    fprintf(fh_debug, '%i %i %f %f %f\n', pair);
    disp = pair(3:5);
    dispd = disp*inv(A);
    % Determine basis atom types in the pair
    itag = pair(1);
    jtag = pair(2);
    id = strucd(itag,:);
    jd = strucd(jtag,:);
    fprintf(fh_debug, 'id: %f %f %f\n', id);
    fprintf(fh_debug, 'jd: %f %f %f\n', jd);
    % Determine basis type of ii
    for ba = 1:nb
      test = id - basisd(ba,:);
      if round(test.*10^4)./(10^4) - round(test) == 0
        bi = ba;
      end
    end

    % Determine basis type of jj
    for ba = 1:nb
      test = jd - basisd(ba,:);
      if round(test.*10^4)./(10^4) - round(test) == 0
        bj = ba;
      end
    end
    nvals = dispd +basisd(bi,:) - basisd(bj,:);
    nvals = [bi round(nvals.*1e4)./1e4 bj];
    fcpairindx = find(fcnl(:,2) == pair(2)); 
    fcpair = fcnl(fcpairindx,:);
    %nvals = fcpair(1,3:7);
    fprintf(fh_debug, '%i %i %i %i %i\n', nvals);
    fcs = fcpair(:,end-1:end);
    [numfcs blah] = size(fcs);
    % Put this pair into the rnl matrix
    rnl(p+1+rnlindx,1:5) = nvals;
    for f = 1:numfcs
      % n values
      rnl(p+1+rnlindx,fcs(f,2)+5) = fcs(f,1);
    end
  end
  rnlindx = rnlindx + N+1;
end

% Now write this RNL
fh_rnl = fopen('RNL','w');
rnl = round(rnl.*1e8)./(1e8);
[numlines blah] = size(rnl);
for a = 1:numlines
  fprintf(fh_rnl, '%i %i %i %i %i %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', rnl(a,:) );
end
fclose(fh_rnl);
fclose(fh_debug);
%dlmwrite('RNL', rnl,' ','precision','%i %i %i %i %i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n');
