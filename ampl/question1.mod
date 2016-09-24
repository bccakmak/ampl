# This is the AMPL program that solves the first question posed on the final project.
# It reads in a training set of Y's and N's to compute a separating line. Assuming this
# line of separation is perfect, that is to say there are no exceptions, all test data
# located under this line is classified as N and all above is classified as Y.

reset;

param xY{1..50, 1..64};
param xN{1..50, 1..64};
param xTest{1..5, 1..64};

var alpha{1..64};
var b;
var e{1..100} >= 0;

param weight{1..100} := 1/50;

maximize distance:
sum{i in 1..100} weight[i] * e[i];

s.t. yDist {i in 1..50}:
e[i] = -b + sum{j in 1..64} alpha[j] * xY[i,j];

s.t. nDist {i in 1..50}:
e[i+50] = b - sum{j in 1..64} alpha[j] * xN[i,j];

s.t. normal_vector:
sum{i in 1..64} alpha[i]^2 <= 1;

data;

read {i in 1..50, j in 1..64} xY[i,j] < Y.txt;
read {i in 1..50, j in 1..64} xN[i,j] < N.txt;
read {i in 1..5, j in 1..64} xTest[i,j] < YesNo.txt;

option solver minos;

solve;

for {i in 1..5} {
	if -b + sum{j in 1..64} alpha[j] * xTest[i,j] >= 0 then {
		printf "Y\n";
	} else {
		printf "N\n";
	}
}