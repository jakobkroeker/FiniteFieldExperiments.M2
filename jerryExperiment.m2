-- look at the hilbert scheme of singular
-- elliptic sextics

restart

uninstallPackage"M2Logging"
installPackage"M2Logging"
check M2Logging

uninstallPackage"IntervalPkg"
installPackage"IntervalPkg"
check IntervalPkg

uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
check BlackBoxIdeals

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
bbC = blackBoxParameterSpace(6*8,K);
-- this Black Box is still empty:
bbC.knownPointProperties()

-- point for testing purposes
testPoint = random(K^1,K^(6*8))

-- from the coefficients make a symmetric 3x3 matrix
-- with linear entries whose 2x2 minors define the sextic
-- elliptic curve. (this is Miles Reids Tom-format)
jerryCubeAt = (point) -> (
     l := 0;
     apply(2,i1->apply(2,i2->apply(2,i3->sum apply(6,k->(
			 coeff := point_l_0;
			 l = l+1;
			 coeff*(symbol x_k)_R
			 ))))
     ))

jerryCubeAt(testPoint)

-- normalize Cube
M = tomMatrixAt(testPoint)

jerryNormalCubeAt = (point) -> (
     C := jerryCubeAt(point);
     J := jacobian (matrix{flatten flatten C}|vars R);
     B := J_{};
     Bplus := B;
     apply(rank source J,i->(
	  Bplus = B|J_{i};
	  if rank Bplus > rank B then B = Bplus
	  ));
     --print rank B;
     varsB = vars R*B^-1;
     apply(C,i->apply(i,j->apply(j,L->sub(L,varsB))))
     )

jerryNormalCubeAt(testPoint)

nullPoint = matrix{{0,0,0,0,0,0}}|random(K^1,K^42)
C = jerryNormalCubeAt(nullPoint)


jerryIdealAt = (point) -> (
     C := jerryNormalCubeAt(point);
     ideal (
	  apply(2,j->det matrix apply(2,i1->apply(2,i2->C#i1#i2#j)))|
	  apply(2,j->det matrix apply(2,i1->apply(2,i2->C#i1#j#i2)))|
	  apply(2,j->det matrix apply(2,i1->apply(2,i2->C#j#i1#i2)))|
	  {
	       det matrix apply(2,i1->apply(2,i2->C#i1#i1#i2)),
     	       det matrix apply(2,i1->apply(2,i2->C#i1#i2#i1)),
     	       det matrix apply(2,i1->apply(2,i2->C#i2#i1#i1))
	  }
     ))

-- test
betti res jerryIdealAt(testPoint)
-- 0: 1 .  . . .
-- 1: . 9 16 9 .
-- 2: . .  . . 1

-- calculate singular locus of cubic
singularLocusAt = (point) -> singularLocus(jerryIdealAt(point))

-- invariants singular locus
codimDegSingularLocusAt = (point) -> (
     s := singularLocusAt(point);
     (codim s,degree s)
     )     

-- register this function in the BlackBox
bbC = bbC.registerPointProperty("codimDegSingularLocusAt",codimDegSingularLocusAt);
bbC.knownPointProperties()
-- {codimDegSingularLocusAt}

e = new Experiment from bbC;
-- so far nothing is observed:
e.recordedProperties()
-- {}

-- NICE TO HAVE: e.wachedProperties()
--e.watchProperty("degreeSingularLocusAt")
e.recordProperty("codimDegSingularLocusAt")
-- now we watch for the degree of the singular locus
e.watchedProperties()

-- lets look at 100 curves
time e.run(100)
-- {(4, 1)} => 4 
-- {(4, 2)} => 2
-- {(4, 3)} => 1
-- {(5, 1)} => 8
-- {(5, 2)} => 26
-- {(5, 3)} => 14
-- {(5, 4)} => 11
-- {(5, 5)} => 2
-- {(5, 13)} => 1
-- {(5, 18)} => 1
-- {(6, 37)} => 40

degComponents = (point) -> sort apply(decompose jerryIdealAt(point),degree)
bbC = bbC.registerPointProperty("degComponents",degComponents);

-- black-box can not be extracted from experiment
tryProperty = (ex,bb,tProp) -> (
     L = ex.pointLists();
     tally flatten apply(keys L, k->apply(L#k,point->(k,(bb.pointProperty(tProp))(point))))
     )

tally apply(primaryDecomposition jerryIdealAt(point),i->(codim i,degree i,degree radical i,betti res i))

-- find example with non expected codimension
codimDegJerryAt = (point) -> (
     codim jerryIdealAt(point),
     degree jerryIdealAt(point)
     )

-- test
codimDegJerryAt(testPoint)

-- look at singular examples found
bbC = bbC.registerPointProperty("codimDegJerryAt",codimDegJerryAt);
tryProperty(e,bbC,"codimDegJerryAt")
-- none of unexpected dimension

eCodim = new Experiment from bbC;
-- so far nothing is observed:
eCodim.recordedProperties()
-- {}


eCodim.recordProperty("codimDegJerryAt")
-- now we watch for the codim of ideal
eCodim.watchedProperties()


time eCodim.run(100000)
-- used 62.7552 seconds
-- {(3, 1)} => 129
-- {(3, 2)} => 33
-- {(3, 3)} => 38
-- {(3, 4)} => 4
-- {(4, 6)} => 71796 

point = (eCodim.pointsByKey({(3,4)}))#0
betti (F = res jerryIdealAt(point))
-- 0: 1 . . .
-- 1: . 6 8 3
F.dd_2

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
f = openOut("tomExperiment1char5.m2")
f << "watchedProperties = "|toString(e.recordedProperties()) << endl
L = e.pointLists();
f << "pointList = new MutableHashTable" << endl
apply(keys L, k->(
	f << ("pointList#"|toString k|" = "|toString(L#k)) << endl
	))  
f << flush
f << close








