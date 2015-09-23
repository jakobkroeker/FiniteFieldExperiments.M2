-- a extremely simple test case
-- for finite field experiments

-- this has to be run after any change in FiniteFieldExperiments
restart

uninstallPackage"M2Logging"
installPackage("M2Logging",UserMode =>true)
check (M2Logging,UserMode =>true)

uninstallPackage"IntervalPkg"
installPackage("IntervalPkg",UserMode =>true)
check (IntervalPkg,UserMode =>true)

uninstallPackage"BlackBoxIdeals"
installPackage("BlackBoxIdeals",UserMode =>true)
check (BlackBoxIdeals,UserMode =>true)

uninstallPackage"FiniteFieldExperiments"
installPackage("FiniteFieldExperiments",UserMode =>true)
check (FiniteFieldExperiments,UserMode =>true)

--viewHelp BlackBoxIdeals
--viewHelp FiniteFieldExperiments


-- here the test case start
restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/5
R = K[x,y,z]

-- ideal of a plane and a line that intersect at the origin
I = ideal (x*z,y*z)

-- make a black box from the ideal
bbI = blackBoxIdeal I;
bbI = new BlackBoxIdeal from I;  -- same as above


print bbI.knownPointProperties()
bbI.knownPointPropertiesAsSymbols()
-- {jacobianAt, rankJacobianAt, isZeroAt, valuesAt, bareJacobianAt}


print bbI.knownMethods()
-- { knownMethods, knownAttributes, knownPointProperties, 
--   knownPointPropertiesAsSymbols, hasPointProperty, pointProperty, 
--   registerPointProperty, updatePointProperty }

print bbI.knownAttributes()
-- {ideal, jacobian, numVariables, type, unknowns, 
--  ring, equations, coefficientRing}
bbI.type       	    -- BlackBoxIdeal
bbI.unknowns 	    -- {x, y, z}
bbI.equations	    -- | xz yz |
bbI.ideal           -- ideal (x*z, y*z)
bbI.ring            -- R
bbI.coefficientRing -- K
bbI.numVariables    -- 3
bbI.numGenerators() -- 2 -- ?why brackets?
bbI.jacobian
-- | z 0 |
-- | 0 z |
-- | x y |

apply(bbI.knownAttributes(),attribute->print (bbI#attribute))


assert (2== bbI.rankJacobianAt(matrix{{0,0,1_K}}))
assert (1== bbI.rankJacobianAt(matrix{{1,2,0_K}}))
assert (0== bbI.rankJacobianAt(matrix{{0,0,0_K}}))
-- this is a point where the ideal does not vanish.
try (bbI.rankJacobianAt(matrix{{1,1,1_K}})) then (error "should fail") else();

assert (bbI.isZeroAt(matrix{{0,0,1_K}}))
assert (bbI.isZeroAt(matrix{{1,2,0_K}}))
assert (bbI.isZeroAt(matrix{{0,0,0_K}}))
-- this is a point where the ideal does not vanish.
assert (not bbI.isZeroAt(matrix{{1,1,1_K}}))



-- make an experiment with rankJacobian
e = new Experiment from bbI
-- test: here rankJacobianAt is automatically watched
e.watchedProperties()
-- niceToHave: make a test like this for a blackbox not from an ideal

-- register new property
bbI.registerPointProperty("rankJacobianAtDup",(bb,point)->(rank bb.jacobianAt(point)))
bbI.knownPointProperties()
-- rankJacobianAtDup is registered
-- bbI.rankJacobianAtDup(matrix{{0,0,1_K}})
-- stdio:73:4:(3): error: key not found in hash table
-- but not accessible

-- the blackbox must be updated
bbI = rebuildBlackBox bbI   

bbI.rankJacobianAtDup(matrix{{0,0,1_K}})
-- now accessible

-- faster:
bbI = bbI.registerPointProperty("rankJacobianAtTrip",(bb,point)->(rank bb.jacobianAt(point)));
bbI.rankJacobianAtTrip(matrix{{0,0,1_K}})
-- immediately accessible

-- shorter
bbI = bbI.rpp("rankJacobianAtQuad",(bb,point)->(4000+rank bb.jacobianAt(point)));


-- if an property is already registerd it can be updated
bbI.updatePointProperty("rankJacobianAtDup",(bb,point)->(1000+rank bb.jacobianAt(point)))
-- no rebuild is necessary
bbI.rankJacobianAtDup(matrix{{0,0,1_K}})
-- 1002


-- nice to have
-- bbI = bbI.rpp("rankJacobianAtTrip",(bb,point)->(rank bb.jacobianAt(point)));

-- Experiment tries to collect a fixed number of points per component
-- this is useful since in an random experiment many more points are found
-- on big components. Not all of them need to be saved.
e.setPointsPerComponent(20)


-- look at 1000 random points
time e.run(1000)
-- {0} => 16 
-- {1} => 385
-- {2} => 57
-- the same information can be accessed via
e.count()
-- these are the number of points found most of these are immedeately forgotten

-- some points with these properties
e.pointLists()
-- only about 20 points for each component are collected:
e.collectedCount() 
-- {0} => 11
-- {1} => 28
-- {2} => 21
-- these are the number of points collected

-- 
-- here the estimated number of components are calculated with rankJacobianAt
-- I do not know what happens when rankJacobianAt does not exist.
--
-- the collected number is somewhat higher than 20 since the maximum 
-- of the confidence interval is used at every step.

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
-- (0,  [0., 0.02]) <= {0}
-- (1,  [0.82, 1.09]) <= {1}
-- (2,  [0.50, 1.05]) <= {2}

-- estimate Decomposition is only useful if jacobianRankAt was watched. 

-- NiceToHave: trailing zeros to correct length

-- second usage of estimate Decomposition
-- functions of this type of usage can be implemented without changing "FiniteFieldExperiments"
estimateDecomposition(e)



maxDeg = 2
mons = matrix {flatten apply(maxDeg+1,i->flatten entries basis(i,R))}

decomposeResult := apply(keys e.pointLists(),k->(
	  L = e.pointsByKey(k);
     	  unique flatten apply(unique L,P->(
	  	    rank bbI.jacobianAt(P);
	  	    time jetP = jetAt(bbI,P,20,1);
		    if not isCertainlySingularAt(bbI,P,20,1) then
	  	    {interpolate(mons,{jetP})}
	  	    ))
     	  ))
decomposeResult
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

