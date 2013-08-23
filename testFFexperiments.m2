-- a extremely simple test case
-- for finite field experiments

-- this has to be run after any change in FiniteFieldExperiments
restart

uninstallPackage"M2Logging"
installPackage"M2Logging"
check M2Logging

uninstallPackage"IntervalPkg"
installPackage"IntervalPkg"
check IntervalPkg

uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
<<<<<<< HEAD
check BlackBoxIdeals

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"
check FiniteFieldExperiments

=======
viewHelp BlackBoxIdeals
viewHelp FiniteFieldExperiments
>>>>>>> blackbox

-- here the test case start
restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/5
R = K[x,y,z]

-- ideal of a plane and a line that intersect at the origin
I = ideal (x*z,y*z)

-- make a black box from the ideal
bbI = blackBoxIdeal I
bbI = new BlackBoxIdeal from I  -- same as above


bbI.knownPointProperties()
bbI.knownPointPropertiesAsSymbols()
<<<<<<< HEAD
--  {rankJacobianAt, rankJacobianAtDup, valuesAt, bareJacobianAt, isZeroAt jacobianAt}

bbI.knownMethods()
-- { knownMethods, knownPointProperties, knownPointPropertiesAsSymbols,  hasPointProperty, pointProperty,
--   registerPointProperty,  updatePointProperty   }

bbI.knownAttributes()
--  {ideal, numVariables, jacobian, numGenerators, ring, type, unknowns, coefficientRing, equations}

=======
>>>>>>> blackbox

assert (2== bbI.rankJacobianAt(matrix{{0,0,1_K}}))
assert (1== bbI.rankJacobianAt(matrix{{1,2,0_K}}))
assert (0== bbI.rankJacobianAt(matrix{{0,0,0_K}}))
-- this is a point where the ideal does not vanish.
<<<<<<< HEAD
assert (2==bbI.rankJacobianAt(matrix{{1,1,1_K}}))


-- make an experiment without rankJacobian at
e = new Experiment from bbI
-- test: here rankJacobianAt is not watched
e.watchedProperties()
-- niceToHave: make a test like this for a blackbox not from an ideal

-- register new property
 bbI.registerPointProperty("rankJacobianAtDup",(bb,point)->(rank bb.jacobianAt(point)))

propSymb:= bbI.knownPointPropertiesAsSymbols()

-- better: nil when not isZero(point) (for rankJacobianAt?)


bbI = rebuildBlackBox bbI   
bbI.knownPointProperties()
propSymb := bbI.knownPointPropertiesAsSymbols() -- not ok any more...
bbI#(propSymb#5) --# ok
bbI.rankJacobianAtDup --ok

-- make an experiment from the black box
-- not necessary!
--e = new Experiment from bbI
-- test: here "rankJacobianAt" is watched
=======
-- the rank here has no meaning
assert (2==bbI.rankJacobianAt(matrix{{1,1,1_K}}))

-- make an experiment without rankJacobian at
e = new Experiment from bbI;
-- test: here rankJacobianAt is automatically watched
>>>>>>> blackbox
e.watchedProperties()

-- look at 1000 random points
e. setMinPointsPerComponent(20)
time e.run(1000)

-- how many times was each wached property oberved?
e.countData()
-- {0} => 16 
-- {1} => 385
-- {2} => 57

-- some points with these properties
e.pointLists()
-- only about 10 points for each component are collected:
apply(keys e.pointLists(),k->print (k => #((e.pointLists())#k))); --jk: not necessary to program; see collectedCount
e.collectedCount() -- jk, Q: rename this method?
-- {0} => 7
-- {1} => 12
-- {2} => 13
-- 
-- here the estimated number of components are calculated with rankJacobianAt
-- I do not know what happens when rankJacobianAt does not exist.
--
-- the collected number is somewhat higher than 10 since the maximum 
-- of the confidence interval is used at every step.
--
-- MANDATORY: Use estimate components to calculate wanted number of points.
-- At the moment this does not work, since estimateNumberOfComponents was rewritten
-- It now takes an experiment and a countKey, at the point where
-- this is coded the experiment does not jet exist, but just "experimentData".
-- If estimateNumberOfComponents is changed to take an experimentData, then
-- it can not be used in estimateDecomposition.
-- WORK AROUND: numberOfComponents now takes a MutableHashTable
-- fixed: the class of experiment is set to 'Experiment' using 'newClass'.


e.estimateStratification()
-----
-- estimated codim <= {wachtched properties}
-----
-- 3.4 <= {0}
-- 2.2 <= {2}
-- 1.1 <= {1}
-----
-- estimated codimension <= {watched properites}
-- niceToHave: see trailing zeros

e.estimateDecomposition()
-- (estimated codim, estimated number of components [confidence interval] <= {watched Properties})
--
-- {0} =>  2.23 [0.27, 4.2]
-- {1} =>  0.96 [0.82, 1.1]
-- {2} =>  1.03 [0.68, 1.39]
--
-- {watched properties} => estimated number of components [confidence interval of estimation]
-- not good

-- estimate Decomposition is only useful if jacobianRankAt was watched. 

-- this does not work since at the point where estimateNumberOfComponents is
-- called experiment is not an Experiment but a MutableHashTable
-- WORK AROUND: numberOfComponents now takes a MutableHashTable -- fixed

-- NiceToHave: trailing zeros to correct length
-- NiceToHave: maybe without showing the center of the estimation

estimateDecomposition(e)
-- I have no Ideal why this does not work -
--
-- (jk) because estimateDecomposition expected a mutable hash table,
-- but the experiment is a nonmutable hash table.  fixed 

-- interpolation using jets
interpolate = (mons,jet) -> if jetP#"succeeded" then (
	       s = syz sub(last coefficients sub(mons,jetP#"jet"),K);
	       I = ideal mingens ideal(mons*s);
	       I
	       )


maxDeg = 2
mons = matrix {flatten apply(maxDeg+1,i->flatten entries super(basis(i,R)))}

decomposeResult := apply(keys e.pointLists(),k->(
	  L = e.pointsByKey(k);
	  -- nice to have: e.pointList(key) 
          --   (jk): same method name with 0 and 1 parameters not easily possible 
          --         opt todo: ask Dan about syntax
     	  unique flatten apply(unique L,P->(
	  	    rank bbI.jacobianAt(P);
	  	    time jetP = JetAt(bbI,P,20,1);
	  	    {interpolate(mons,jetP)}
	  	    ))
     	  ))

assert(decomposeResult==={{null}, {ideal(z)}, {ideal (y, x)}});

-- better: 
--    0) empty ideal list
--    1) for each point make a short jet (eg. length 10) and check if it is contained
-- 	 in any of the ideals in the list 
--    2) if not, make a very long jet, interpolate, add ideal to list
--
-- even better: make several long, but not very long, jets by some heuristic
--
-- possibly: test wether interpolated ideal + I == interpolated ideal. This
--           only makes sense, if the ideal is not given as a blackbox
--
-- possibly: ideal list could be stored in the experiment. 
-- possibly: a watchable property could be the ideal which contains a jet of
--           a found point. 



-- register new property
bbI.registerPointProperty("rankJacobianAtDup",(bb,point)->(rank bb.jacobianAt(point)))

propSymb:= bbI.knownPointPropertiesAsSymbols()

-- better: nil when not isZero(point) (for rankJacobianAt?)


bbI = rebuildBlackBox bbI   
propSymb := bbI.knownPointPropertiesAsSymbols()
bbI#(propSymb#5) --# ok
bbI.rankJacobianAtDup --ok


-- stdio:20:4:(3): error: key not found in hash table
-- does not work. Not exported?
--fixed

-- niceToHave: make a test like this for a blackbox not from an ideal
