-- a extremely simple test case
-- for finite field experiments

-- this has to be run after any change in FiniteFieldExperiments
restart
uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"
uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"

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
bbI.knownPointProperties()
bbI.knownPointPropertiesAsSymbols()
-- {isZeroAt, jacobianAt,, valuesAt, valuesAt, jacobianAt, bareJacobianAt, bareJacobianAt}
-- ERROR: Why to properties appear several times???
---- no, they do not appear several times , Macaulay2 does print key strings without quotes, so in fact we can
--access the properties by .jacobianAt (fragile, because of how M2 deals with symbols) and #"jacobianAt" 
-- fixed, by returning only strings.

-- MANDATORY: rankJacobianAt must exist in BlackBoxIdeals
-- Fixed
bbI.rankJacobianAt(matrix{{0,0,1_K}})
bbI.rankJacobianAt(matrix{{1,2,0_K}})
bbI.rankJacobianAt(matrix{{0,0,0_K}})

-- register new property
 bbI.registerPointProperty("rankJacobianAtDup",(bb,point)->(rank bb.jacobianAt(point)))

propSymb:= bbI.knownPointPropertiesAsSymbols()

-- better: nil when not isZero(point)
assert (2==(bbI.pointProperty("rankJacobianAt"))(matrix{{0,0,1_K}}))
assert (1==(bbI.pointProperty("rankJacobianAt"))(matrix{{1_K,2_K,0_K}}))
assert (0==(bbI.pointProperty("rankJacobianAt"))(matrix{{0,0,0_K}}))
-- OK

bbI = rebuildBlackBox bbI   
propSymb := bbI.knownPointPropertiesAsSymbols()
bbI#(propSymb#5) --# ok
bbI.rankJacobianAtDup --ok


-- stdio:20:4:(3): error: key not found in hash table
-- does not work. Not exported?
--fixed

-- make an experiment from the black box
e = new Experiment from bbI

-- watch rank of jacobi matrix
-- MAYOR: This should be set automatically when making a new experiment from an ideal
e.clear()
e.watchProperties({"isZeroAt"})
e.watchedProperties()
e.watchProperties({"rankJacobianAt"})
e.watchedProperties()

-- look at 1000 random points
time e.run(1000)

-- how many times was each wached property oberved?
e.countData()
-- {0} => 16 
-- {1} => 385
-- {2} => 57

-- some points with these properties
e.pointLists()
-- only about 20 points for each component are collected:
apply(keys e.pointLists(),k->print (k => #((e.pointLists())#k)));
-- {0} => 6
-- {1} => 17
-- {2} => 18
-- 
-- here the estimated number of components are calculated without
-- using jacobiRankAt. If jacobiRankAt is watched, it should also be used 
-- (see estimateDecomposition below)
--
-- the estimated number of components for {1} and {2} is somewhat lower than 1
-- since there is the singular origin that lies on both components, but is
-- counted at {0}

e.estimateStratification2()
-----
-- 3.4 <= {0}
-- 2.2 <= {2}
-- 1.1 <= {1}
-----
-- estimated codimension <= {watched properites}

e.estimateStratification()
-- {0} =>  4.35 [3, 5.7]}
-- {1} =>  1 [0.9, 1.1]
-- {2} =>  2.2 [2, 2.4]

-- why are estimates different?
-- ERROR: problem with new interval. Center of interval is calculated not given. 
--        In a poisson Interval the center is not the average of the max and min

-- In applications Stratification (with intervals) was to much information. 
-- Stratification2 turned out to be much more useful 

e.estimateDecomposition()
-- {0} =>  2.23 [0.27, 4.2]
-- {1} =>  0.96 [0.82, 1.1]
-- {2} =>  1.03 [0.68, 1.39]
--
-- {watched properties} => estimated number of components [confidence interval of estimation]
-- not good

-- estimate Decomposition is only useful if jacobianRankAt was watched. 
-- here comes a rudementary implementation
assert ((e.watchedProperties())#0 == "rankJacobianAt")
-- later: find position of "rankJacobianAt"
c = e.countData()
print "(estimated codim, estimated number of Components <= {watched Properties})"
apply(sort apply(keys c,k->(net(k#0,estimateNumberOfComponents(e.trials(),k#0,c#k,char K)))|" <= "|net k),print)
-- (0,  0.01 [0., 0.01]) <= {0}
-- (1,  0.94 [0.84, 1.03]) <= {1}
-- (2,  0.74 [0.55, 0.93]) <= {2}
-- 
-- (est codim, est number of Components <= {watched Properties})
--
-- NiceToHave: trailing zeros to correct length
-- NiceToHave: maybe without showing the center of the estimation


-- interpolation using jets
interpolate = (mons,jet) -> if jetP#"succeeded" then (
	       s = syz sub(last coefficients sub(mons,jetP#"jet"),K);
	       I = ideal mingens ideal(mons*s);
	       I
	       )


maxDeg = 2
mons = matrix {flatten apply(maxDeg+1,i->flatten entries super(basis(i,R)))}

apply(keys e.pointLists(),k->(
	  L = (e.pointLists())#k;
	  -- nice to have: e.pointList(key)
     	  unique flatten apply(unique L,P->(
	  	    rank bbI.jacobianAt(P);
	  	    time jetP = JetAt(bbI,P,20,1);
	  	    {interpolate(mons,jetP)}
	  	    ))
     	  ))
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
