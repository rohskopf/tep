%%% Use after "strucgen.m"
%%% Creates a file with direct coordinate values

% Inputs
% None

posarr = posarr/bx;

[N blah] = size(posarr);

fh = fopen('directcoor', 'w');
fprintf(fh, '%i\n', N);
fprintf(fh, '%.10f\n', 1);
fprintf(fh, '%.10f %.10f %.10f\n', bx(1,:));
fprintf(fh, '%.10f %.10f %.10f\n', bx(2,:));
fprintf(fh, '%.10f %.10f %.10f\n', bx(3,:));

[N blah] = size(posarr);

one = ones(64,1);
append = [one]; %[one; two];
arr = [append posarr];
for l = 1:N
  fprintf(fh, '%i %.10f %.10f %.10f\n', arr(l,:));

end

fclose(fh);
