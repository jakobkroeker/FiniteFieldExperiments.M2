-- look at the hilbert scheme of singular
-- cubics

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

K = ZZ/7
R = K[x,y,z,w]
-- there are 20 monomials in degree 3
mons3 = matrix entries transpose super basis(3,R)

-- coefficients of x^3+y^3+z^3 for testing purposes
cubicCone = matrix{{x^3+y^3+z^3}}
coeffCubicCone = contract(transpose mons3,cubicCone)
-- it has one singular point at (0:0:0:1)

-- from the coefficients make a cubic
cubicAt = (point) -> point*mons3
cubicAt(coeffCubicCone) == cubicCone
-- calculate singular locus of cubic
singularLocusAt = (point) -> ideal jacobian cubicAt(point)
betti res singularLocusAt(coeffCubicCone)
-- calculate number of singular points
degreeSingularLocusAt = (point) -> (
     s := singularLocusAt(point);
     if dim s == 0 then return 0;
     if dim s == 1 then return degree s;
     if dim s > 2 then return infinity
     -- these are affine dimensions
     )
degreeSingularLocusAt(coeffCubicCone)

bb = blackBoxParameterSpace(20,K)
bb= bb.registerPointProperty("degreeSingularLocusAt", degreeSingularLocusAt);


e = new Experiment from bb;
-- e.recordedProperties()
e.recordProperty("degreeSingularLocusAt")
--tally apply(100,i->degreeSingularLocusAt(random(K^1,K^20)))
e.run(100)
-- e.countData()

multiplicitiesSingularLocusAt = (point) -> (
     sing := singularLocusAt(point);
     d := degreeSingularLocusAt(point);
     if d == 0 then return {};
     if d == infinity then return {infinity};
     pSing := primaryDecomposition(sing);
     (d,flatten apply(pSing,i->(
	       if dim i == 0 then return {};
	       r := degree radical i;
	       d := degree i;
	       apply(d//r,j->r)
	       )))
     )

bb= bb.registerPointProperty("multiplicitiesSingularLocusAt", multiplicitiesSingularLocusAt);
e.clear() --should it also clear the list of recorded properties?
e.clearRecordList()
e.recordProperty( "multiplicitiesSingularLocusAt" )
-- tally apply(1250,i->multiplicitiesSingularLocusAt(random(K^1,K^20)))
time e.run(1250)-- should it always return countData?
--e.countData()

    
     

