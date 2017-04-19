%%% Create a lattice based on lattice vectors

%%% INPUTS %%%
a = 5.472; % Lattice constant
A = (a/2).*[0 1 1; 1 0 1; 1 1 0]; % Lattice vector matrix
basisc = (a/4).*[0 0 0; 1 1 1]; % Basis in cartesian coordinates
nvals = 10; % Scope of n values to check
supercell = [2; 2; 2]; % 2x2x2 supercell
bx = [a 0 0; 
      0 a 0; 
      0 0 a]; % Box shape

%%% CODE (don't edit) %%%
bx = bx.*supercell; 


a1 = A(1,:);
a2 = A(2,:);
a3 = A(3,:);

[batoms blah] = size(basisc);

posarr = [];
hold on
for b = 1:batoms
  for n1 = -nvals:nvals
    for n2 = -nvals:nvals
      for n3 = -nvals:nvals
        r = n1.*a1 + n2.*a2 + n3.*a3 + basisc(b,:);
        if 0 <= r(1) && r(1) < bx(1,1) && 0 <= r(2) && r(2) < bx(2,2) && 0 <= r(3) && r(3) < bx(3,3)
          scatter3(r(1), r(2), r(3));
          posarr = [posarr; r];
        end
      end
    end
  end
end
hold off
