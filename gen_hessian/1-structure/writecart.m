%%% Use after "strucgen.m"
%%% Creates a file with direct coordinate values

% Inputs
% None

%posarr = posarr/boxvectors;

[N blah] = size(posarr);

fh = fopen('cartcoor', 'w');
fprintf(fh, '%i\n', N);
fprintf(fh, '%.10f\n', 1);
fprintf(fh, '%.10f %.10f %.10f\n', bx(1,:));
fprintf(fh, '%.10f %.10f %.10f\n', bx(2,:));
fprintf(fh, '%.10f %.10f %.10f\n', bx(3,:));

[N blah] = size(posarr);

ids = 1:N;
ids = ids';
types = [ones(N,1)];
append = [];%[ids types];
arr = [ append posarr];
for l = 1:N
  fprintf(fh, '%.10f %.10f %.10f\n', arr(l,:));

end

fclose(fh);
