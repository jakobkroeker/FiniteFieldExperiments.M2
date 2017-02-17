-- look at the hilbert scheme of singular
-- elliptic sextics

quit -- F11 F11 F12

path = append(path,"/Users/bothmer/Desktop/projekte/strudel/Jakob2010/GitHub/padicLiftM2/")

uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
viewHelp BlackBoxIdeals
check BlackBoxIdeals



uninstallPackage"M2Logging"
installPackage"M2Logging"
check M2Logging

uninstallPackage"IntervalPkg"
installPackage"IntervalPkg"
check IntervalPkg


uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"
check FiniteFieldExperiments

viewHelp BlackBoxIdeals
viewHelp FiniteFieldExperiments


-- here the test case start
restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/3
R = K[x_0..x_5]

-- make a blackbox describing the parameter space of cubicc
bbC = blackBoxParameterSpace(6*9,K);
-- this Black Box is still empty:
bbC.knownPointProperties()

-- point for testing purposes
testPoint = random(K^1,K^54)

-- from the coefficients make a symmetric 3x3 matrix
-- with linear entries whose 2x2 minors define the sextic
-- elliptic curve. (this is Miles Reids Tom-format)
tomMatrixAt = (point) -> (
     l := 0;
     matrix apply(3,i->apply(3,j->sum apply(6,k->(
			 coeff := point_l_0;
			 l = l+1;
			 coeff*(symbol x_k)_R
			 ))))
     )

tomMatrixAt(random(K^1,K^54))

-- normalize Matrix
M = tomMatrixAt(testPoint)

tomNormalMatrixAt = (point) -> (
     M := tomMatrixAt(point);
     J := jacobian (matrix{flatten entries M}|vars R);
     B := J_{};
     Bplus := B;
     apply(rank source J,i->(
	  Bplus = B|J_{i};
	  if rank Bplus > rank B then B = Bplus
	  ));
     --print rank B;
     sub(M,vars R*B^-1)
     )


tomNormalMatrixAt(testPoint)
nullPoint = matrix{{0,0,0,0,0,0}}|random(K^1,K^48)
tomNormalMatrixAt(nullPoint)


tomIdealAt = (point) -> minors(2,tomNormalMatrixAt(point))
-- test
betti res tomIdealAt(random(K^1,K^54))
-- 0: 1 .  . . .
-- 1: . 9 16 9 .
-- 2: . .  . . 1

-- calculate singular locus of cubic
singularLocusAt = (point) -> singularLocus(tomIdealAt(point))

-- invariants singular locus
codimDegSingularLocusAt = (point) -> (
     s := singularLocusAt(point);
     (codim s,degree s)
     )     

-- register this function in the BlackBox
bbC = bbC.registerPointProperty codimDegSingularLocusAt;
pointProperties bbC
-- {codimDegSingularLocusAt}

e = new Experiment from bbC;
-- so far nothing is observed:
watchedProperties e
-- {}

--e.watchProperty("degreeSingularLocusAt")
e.watchProperty "codimDegSingularLocusAt"
-- now we watch for the degree of the singular locus
watchedProperties e

-- lets look at 100 curves
time e.run 1
time e.run 10 -- 30 sec on T42
time e.run 100
-- {(4,  1)} => 10   
-- {(4,  2)} => 4
-- {(5,  1)} => 493
-- {(5,  2)} => 655
-- {(5,  3)} => 270
-- {(5,  4)} => 37
-- {(5,  5)} => 3
-- {(5, 18)} => 6
-- {(6, 37)} => 3422
degComponents = (point) -> sort apply(decompose tomIdealAt(point),degree)
bbC = bbC.registerPointProperty("degComponents",degComponents);

e.tryProperty "degComponents" 
point = (e.pointsByKey({(5,2)}))#0

-- find example with non expected codimension
codimDegTomAt = (point) -> (
     codim tomIdealAt(point),
     degree tomIdealAt(point)
     )

-- test
codimDegTomAt(testPoint)

-- look at singular examples found
bbC = bbC.registerPointProperty codimDegTomAt;
e.tryProperty("codimDegTomAt")
-- none of unexpected dimension

eCodim = new Experiment from bbC;
-- so far nothing is observed:
eCodim.watchedProperties()
-- {}


eCodim.watchProperty("codimDegTomAt")
-- now we watch for the codim of ideal
eCodim.watchedProperties()

time eCodim.run(100)
time eCodim.run(1000)
-- used 62.7552 seconds
-- {(3, 1)} => 129
-- {(3, 2)} => 33
-- {(3, 3)} => 38
-- {(3, 4)} => 4
-- {(4, 6)} => 71796 

point = (eCodim.pointsByKey({(3,4)}))#1
betti res tomIdealAt(point)
-- 0: 1 . . .
-- 1: . 6 8 3


dI = primaryDecomposition(tomIdealAt(point))
tally apply(dI,i->(codim i,degree i, degree radical i,betti res i))
tomNormalMatrixAt(point)

-- direct matrix construction
dI = decompose tomIdealAt(point)
betti res dI#0
codim singularLocus dI#0
M = tomNormalMatrixAt(point)
sub(M,{x_0=>0,x_1=>0,x_2=>0,x_3=>x_4})

M = matrix{
     {0,x_1,x_2},
     {x_0,random(1,R),random(1,R)},     
     {x_3,random(1,R),random(1,R)}
     }

betti (res (I = minors(2,M)))
dI = primaryDecomposition I
apply(dI,i->(degree radical i,degree i))
dI#0
dI#2

-- save experiment data

f = openOut("tomExperiment1char5.newIO.m2")
f << toExternalString e.experimentData()
f << flush
f << close

g = openIn("tomExperiment1char5.newIO.m2")
experimentDataString = read g
experimentDataString

f = openOut("tomExperiment1char5.m2")
f << "watchedProperties = "|toString(e.watchedProperties()) << endl
L = e.pointLists();
f << "pointList = new MutableHashTable" << endl
apply(keys L, k->(
	f << ("pointList#"|toString k|" = "|toString(L#k)) << endl
	))  
f << flush
f << close

