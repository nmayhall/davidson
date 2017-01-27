clc;
clear all;
close all;

# ******************************************
# Davidson Test (davidsonTest.m)
#  Last Edited: Jan 26, 2017 (Shannon Houck)
#
# Tests the Davidson algorithm (davidson.m)
#  and compares results to expected values.
#
# ******************************************
rand('seed',2);
# define dimensions of test matrix
dim = 4000;

# define number of eigenvalues to find
eigVals = 3;

# generate random test matrix
ATest = rand(dim,dim);
ATest = ATest'+ATest; # symmetrize
ATest = ATest + diag(eye(dim,1)) - 5000*eye(dim);

# form initial search subspace
vTest = orth(rand(dim,eigVals));

# define convergence criteria
cutTest = 1e-6; # cutoff
maxTest = 300; # max iteration

printf("Mine:\n");
davidson(ATest, vTest, cutTest, maxTest);

printf("\n\n Expected:\n")
[u,l] = eigs(ATest, eigVals);
for s=diag(l)
    printf(" %16.12f\n",s);
end
printf("\n");
