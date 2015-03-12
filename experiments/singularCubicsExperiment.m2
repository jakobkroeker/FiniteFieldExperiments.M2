-- look at the hilbert scheme of singular
-- cubics

restart

path = append(path,"/Users/bothmer/Desktop/projekte/strudel/Jakob2010/GitHub/padicLiftM2/")

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

K = ZZ/7
R = K[x,y,z,w]

-- make a blackbox describing the parameter space of cubicc
bbC = blackBoxParameterSpace(20,K);
-- this Black Box is still empty:
bbC.knownPointProperties()

-- there are 20 monomials in degree 3
mons3 = matrix entries transpose super basis(3,R)

-- coefficients of x^3+y^3+z^3 for testing purposes
cubicCone = matrix{{x^3+y^3+z^3}}
coeffCubicCone = contract(transpose mons3,cubicCone)
-- it has one singular point at (0:0:0:1)

-- coefficients of the smooth fermat cubic for testing purposes
cubicFermat = matrix{{x^3+y^3+z^3+w^3}}
coeffCubicFermat = contract(transpose mons3,cubicFermat)

-- from the coefficients make a cubic
cubicAt = (point) -> matrix entries (point*mons3)
-- test: do we get the cone-equations from the cone coefficients?
assert (cubicAt(coeffCubicCone) == cubicCone)
assert (cubicAt(coeffCubicFermat) == cubicFermat)

-- register this function in the BlackBox
bbC = bbC.registerPointProperty("cubicAt",cubicAt);
bbC.knownPointProperties()
-- {cubicAt}

-- calculate singular locus of cubic
singularLocusAt = (bb,point) -> ideal jacobian bb.cubicAt(point)
bbC = bbC.registerPointProperty("singularLocusAt",singularLocusAt);

-- Test: is the multiplicity of the singular locus of the cubic cone correct?
assert (8==degree bbC.singularLocusAt(coeffCubicCone))
-- Test: is the singular locus of the cubic cone in codim 3 (i.e. points)?
assert (3==codim bbC.singularLocusAt(coeffCubicCone))
-- Test: is the fermat cubic smooth?
assert (4==codim bbC.singularLocusAt(coeffCubicFermat))

-- calculate number of singular points
-- special treatment for cubics with 0 or infinitely many singular points
degreeSingularLocusAt = (bb,point) -> (
     s := bb.singularLocusAt(point);
     if dim s == 0 then return 0;
     if dim s == 1 then return degree s;
     if dim s >= 2 then return infinity
     -- these are affine dimensions
     )
-- check if cubic cone is treated correctly
assert (8==degreeSingularLocusAt(bbC,coeffCubicCone))
-- check if a smooth cubic is treated correctly
assert (0==degreeSingularLocusAt(bbC,coeffCubicFermat))

-- register this function
bbC = bbC.registerPointProperty("degreeSingularLocusAt",degreeSingularLocusAt);
-- MANDATORY: It must be possible to register updated versions of properties
-- since when developing registered properties often contain bugs.
-- SOLVED: Use updatePointProperty

-- now we can study the distribution of singular cubics in the parameter
-- space of all cubics:

e = new Experiment from bbC;
-- so far nothing is observed:
e.watchedProperties()
-- {}

e.watchProperty("degreeSingularLocusAt")
-- now we watch for the degree of the singular locus
e.watchedProperties()
-- {degreeSingularLocusAt}


-- lets look at 100 cubics
e.run(100)
-- {0} => 76}
-- {1} => 14
-- {2} => 9
-- {3} => 1

-- lets look at some more
e.run(100)
-- {0} => 160
-- {1} => 24
-- {2} => 15
-- {3} => 1

-- what does this mean? Lets try to estimate the stratification of 
-- the parameter space of cubics by this property:
e.estimateStratification()
--
-- 2.7 <= {3}
-- 1.3 <= {2}
-- 1.1 <= {1}
-- .1 <= {0}
--
-- it seems that:
--     1) a dense open subset (of codim 0) contains smooth cubics
--     2) the cubics with degree 1 singularity form a divisor
--     3) the remaining estimates are strange

-- does it help to look at more points?
e.run(200)
e.estimateStratification()
--
-- 2.7 <= {3}
-- 1.5 <= {2}
-- 1 <= {1}
-- .1 <= {0}
--

-- no, the degree 2 statum still has a codimension between 1 and 2
-- possibly this locus has several components?

-- lets look at the multiplicity structure of the singular locus
multiplicitiesSingularLocusAt = (bb,point) -> (
     sing := bb.singularLocusAt(point);
     d := bb.degreeSingularLocusAt(point);
     if d == 0 then return {};
     if d == infinity then return {infinity};
     pSing := primaryDecomposition(sing);
     sort flatten apply(pSing,i->(
	       if dim i == 0 then return {};
	       r := degree radical i;
	       d := degree i;
	       apply(r,j->d//r)
	       ))
     )

-- Test: The cubic cone has one singular point of mulitplicity 8
assert ({8}==multiplicitiesSingularLocusAt(bbC,coeffCubicCone))
-- Test: The fermat cubic has no singular points
assert ({}==multiplicitiesSingularLocusAt(bbC,coeffCubicFermat))

-- register this property
bbC = bbC.registerPointProperty("multiplicitiesSingularLocusAt",multiplicitiesSingularLocusAt);

-- try this out on the collected points
e.collectedCount()

e.tryProperty("multiplicitiesSingularLocusAt")
-- ({0}, {}) => 10      }
-- ({1}, {1}) => 10
-- ({2}, {1, 1}) => 4
-- ({2}, {2}) => 6
-- ({3}, {1, 1, 1}) => 1
-- ({3}, {2, 1}) => 1
--
-- is seems that this property divides the degree 2 locus in two of similar size
-- lets run the experiment again and watch this property
e.clear() 
-- this erases statistics
e.collectedCount()
e.countData()
-- NICE TO HAVE: e.clear does not need to erase collected points. They only
-- need to be reclassified if wached properties are added or removed
e.watchedProperties()
-- {degreeSingularLocusAt}
e.watchProperty("multiplicitiesSingularLocusAt")
e.watchedProperties()
-- {degreeSingularLocusAt, multiplicitiesSingularLocusAt}

time e.run(1000)
-- used 9.04075 seconds
--
-- {0, {}} => 817        
-- {1, {1}} => 131
-- {2, {1, 1}} => 21
-- {2, {2}} => 20
-- {3, {1, 1, 1}} => 2
-- {3, {1, 2}} => 3
-- {3, {3}} => 3
-- {4, {1, 1, 1, 1}} => 1
-- {4, {1, 3}} => 1
-- {4, {4}} => 1

e.estimateStratification()
--
-- .1 <= {0, {}}
-- 1 <= {1, {1}}
-- 2 <= {2, {1, 1}}
-- 2 <= {2, {2}}
-- 3 <= {3, {3}}
-- 3 <= {3, {1, 2}}
-- 3.2 <= {3, {1, 1, 1}}
-- 3.5 <= {4, {4}}
-- 3.5 <= {4, {1, 3}}
-- 3.5 <= {4, {1, 1, 1, 1}}
--
-- indeed it seems that for every partition of 2 and 3 we get an
-- irreducible component of codimension 2 and 3 respectively
-- 
-- lets do more points
time e.run(10000)    
-- used 94.8239 seconds
e.countData()
e.estimateStratification()     

-- notice that we only count point on the open stratum corresponding
-- to each partition. The number of points on this stratum reflects
-- the number of boundary strata, i.e. strata of the same dimension with 
-- more boundary strata have less points.

-- Notice also that for the partition {4} there seem to be two components.
-- Indeed this is the case: one correspoinding to the Dynkin/
-- diagramm A_4 and one to D_4.

-- is singular locus for multiplicity {1,1,1} indeed reduced?
tally apply(e.pointsByKey({3, {1,1,1}}),coeff->(
	  sing = singularLocusAt(bbC,coeff);
	  (
	       degree sing,
	       degree radical sing
	  )
     ))
-- yes

-- look at decomposition over finite field (not over the algebraic closure)
tally apply(e.pointsByKey({3, {1,1,1}}),coeff->(
	  sing = singularLocusAt(bbC,coeff);
	  apply(decompose sing,i->degree i)
	  ))
-- {1, 1, 1} => 1}
-- {1, 2} => 5
-- {3} => 3

-- indeed we have approximately the relation 1:3:2 corresponding
-- to the distribution of cycle lengths in S_3

