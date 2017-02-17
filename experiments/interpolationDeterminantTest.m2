-- find determinants
-- via interpolation


restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

prime = 11
n = 3 -- number of variables
m = 4 -- size of matrix
d = 1 -- degree of entries

-- the field
K = ZZ/prime
-- the ring of entries
R = K[x_1..x_n]
apply(20,m->(
	  print m;
	  -- a random m x m matrix with degree d entries in n variables
	  M = random(R^{m:0},R^{m:-d});
	  time det M;
     )
)

-- used 1.55608 seconds (m=16,n=3)
-- used 2.28799 seconds (m=17,n=3)
-- used 3.3404 seconds (m=18,n=3)
-- used 4.68447 seconds (m=19,n=3)

m=10
M = random(R^{m:0},R^{m:-d})

-- evaluate the determinant at a point
detM = point -> matrix{{det sub(M,point)}}
-- !!!result must be a matrix!!!

S=K[y_{0,0}..y_{m-1,m-1}]
RS = R**S

point = matrix{{1,5,1_K}}

-- jacobianAt for Matrices
quickCoRankJacobianAt = (point) -> (
     r = numcols M - rank sub(M,point);
     MY = ideal (sub(M,RS)*sub(genericMatrix(S,m,r),RS));
     pointY = point|flatten syz sub(M,point)|matrix{{m*(m-r):0_K}};
     assert (0==sub(MY,pointY));
     print"kjhgkjhgkj";
     m-rank sub((jacobian(MY))^{0..m-1},pointY)
     )


-- make a black box from the evaluation function
bb = blackBoxIdealFromEvaluation(R,detM);
pointProperties bb
-- register quick Jacobian
bb = bb.upp("rankJacobianAt",quickCoRankJacobianAt);
pointProperties bb


-- find points on the determinant
e = new Experiment from bb;
e.watchedProperties()

time e.run 10 
-- !!! manchmal sehr langsam !!!
-- used 48.6921 seconds (for 10)

point = first points e
bb.rankJacobianAt point

e.watchedProperties()
e.counts()

time e.createAllInterpolatedIdeals(4,2)

E = K[ee]   --/(ee^2)

point = sum apply(40,i->random(E^{0},E^{3:-i}));
time det sub(M,point);
-- used 0.183348 seconds (10, m=10, n=3)
-- used 0.831813 seconds (20)
-- used 1.94263 seconds (30)
-- used 3.4279 seconds (40)
time det M;
-- used 0.094525 seconds

point1 = (e.pointsByKey({1}))#0
M1 = sub(M,point1)
syz M1
M1*(
     (M1^{0..8}_{0..8})^-1
     |syz M1)

time detM(sub(point,K))
M1 = sub(M,point)
--M1*M1^-1
time 
det sub(M,point)
S = ring eps


e.printInterpolatedIdeals()
ideal det M
keys e

ed = e.experimentData()
ht = new HashTable from ed.interpolatedIdeals
detInterpol = (ht#"ideal_0").ideal
-- !!! this must work automatically !!!
assert (detInterpol == ideal det M)

