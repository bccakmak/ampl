# This is the AMPL program that solves the second question posed on the final project.
# It reads in a training set of characters from A to Z (and space) to compute a
# separating hyperplane for each. And then, the code assigns a character to each 8x8 image
# vector in the test files "Slogan3", "Slogan4", "Slogan 5", depending on their distance
# to the respective hyperplane.

reset;

set alphabet := {"ABCDEFGHIJKLMNOPQRSTUVWXYZ "};
param m; # indicates which character is used from above to solve the optimization problem
param likelihood; # a conservative measure that allows us to take the mathematically
# closest data point in case a hyperplane isn't perfect and for one data point, there are
# no hyperplanes which is located below the data point. Also lets us avoid printing
# multiple letters for one data point by just taking the best (mathematically maximum
# distance) option.

param xData{1..50*27,1..64};

param xTest3{1..23,1..64};
param xTest4{1..39,1..64};
param xTest5{1..43,1..64};

var alpha{1..64, 1..27};
var b{1..27};
var e{1..50*27, 1..27} >= 0;

param weightCurrent := 1/50;
param weightOthers := 1/(50*26);

maximize distance:
sum{i in 1..50} weightCurrent * e[i,m] + sum{j in 51..50*27} weightOthers * e[j,m];

s.t. otherDist1 {i in 1..50*(m-1)}:
e[i,m] = b[m] - sum{j in 1..64} alpha[j,m] * xData[i,j];

s.t. currDist {i in 50*(m-1)+1..50*m}:
e[i,m] = -b[m] + sum{j in 1..64} alpha[j,m] * xData[i,j];

s.t. otherDist2 {i in 50*m+1..50*27}:
e[i,m] = b[m] - sum{j in 1..64} alpha[j,m] * xData[i,j];

s.t. normal_vector:
sum{i in 1..64} alpha[i,m]^2 <= 1;

data;

read {i in 1..50, j in 1..64} xData[i,j] < A.txt;
read {i in 51..100, j in 1..64} xData[i,j] < B.txt;
read {i in 101..150, j in 1..64} xData[i,j] < C.txt;
read {i in 151..200, j in 1..64} xData[i,j] < D.txt;
read {i in 201..250, j in 1..64} xData[i,j] < E.txt;
read {i in 251..300, j in 1..64} xData[i,j] < F.txt;
read {i in 301..350, j in 1..64} xData[i,j] < G.txt;
read {i in 351..400, j in 1..64} xData[i,j] < H.txt;
read {i in 401..450, j in 1..64} xData[i,j] < I.txt;
read {i in 451..500, j in 1..64} xData[i,j] < J.txt;
read {i in 501..550, j in 1..64} xData[i,j] < K.txt;
read {i in 551..600, j in 1..64} xData[i,j] < L.txt;
read {i in 601..650, j in 1..64} xData[i,j] < M.txt;
read {i in 651..700, j in 1..64} xData[i,j] < N.txt;
read {i in 701..750, j in 1..64} xData[i,j] < O.txt;
read {i in 751..800, j in 1..64} xData[i,j] < P.txt;
read {i in 801..850, j in 1..64} xData[i,j] < Q.txt;
read {i in 851..900, j in 1..64} xData[i,j] < R.txt;
read {i in 901..950, j in 1..64} xData[i,j] < S.txt;
read {i in 951..1000, j in 1..64} xData[i,j] < T.txt;
read {i in 1001..1050, j in 1..64} xData[i,j] < U.txt;
read {i in 1051..1100, j in 1..64} xData[i,j] < V.txt;
read {i in 1101..1150, j in 1..64} xData[i,j] < W.txt;
read {i in 1151..1200, j in 1..64} xData[i,j] < X.txt;
read {i in 1201..1250, j in 1..64} xData[i,j] < Y.txt;
read {i in 1251..1300, j in 1..64} xData[i,j] < Z.txt;
read {i in 1301..1350, j in 1..64} xData[i,j] < Space.txt;

read {i in 1..23, j in 1..64} xTest3[i,j] < Slogan3.txt;
read {i in 1..39, j in 1..64} xTest4[i,j] < Slogan4.txt;
read {i in 1..43, j in 1..64} xTest5[i,j] < Slogan5.txt;

for {o in 1..27} {

let m := o;

option solver knitro;
solve;

}

for {i in 1..23} {
	let likelihood := max({j in 1..27} -b[j] + sum{k in 1..64} alpha[k,j] * xTest3[i,k]);
	for {j in 1..27} {
		if -b[j] + sum{k in 1..64} alpha[k,j] * xTest3[i,k] = likelihood then {
			printf "%s", {t in alphabet} substr(t,j,1);
		}
	}
}

printf "\n";

for {i in 1..39} {
	let likelihood := max({j in 1..27} -b[j] + sum{k in 1..64} alpha[k,j] * xTest4[i,k]);
	for {j in 1..27} {
		if -b[j] + sum{k in 1..64} alpha[k,j] * xTest4[i,k] = likelihood then {
			printf "%s", {t in alphabet} substr(t,j,1);
		}
	}
}

printf "\n";

for {i in 1..43} {
	let likelihood := max({j in 1..27} -b[j] + sum{k in 1..64} alpha[k,j] * xTest5[i,k]);
	for {j in 1..27} {
		if -b[j] + sum{k in 1..64} alpha[k,j] * xTest5[i,k] = likelihood then {
			printf "%s", {t in alphabet} substr(t,j,1);
		}
	}
}

printf "\n";