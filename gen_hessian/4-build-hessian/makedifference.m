%%% Take 2 IILIST and IJLIST files (DFT and potential), subtract them, and then make a file that fills in the difference
iilist_dft = dlmread('IILIST_dft');
iilist_tbc1 = dlmread('IILIST_tbc1');
ijlist_dft = dlmread('IJLIST_dft');
ijlist_tbc1 = dlmread('IJLIST_tbc1');

iilist_diff = iilist_dft(:,5) - iilist_tbc1(:,5);
ijlist_diff = ijlist_dft(:,5) - ijlist_tbc1(:,5);

iilist_diff = [iilist_dft(:,1:4) iilist_diff];
ijlist_diff = [ijlist_dft(:,1:4) ijlist_diff];

dlmwrite('IILIST', iilist_diff, ' ');
dlmwrite('IJLIST', ijlist_diff, ' ');

