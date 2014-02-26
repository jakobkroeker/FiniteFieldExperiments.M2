-- find polynomials containing a component
-- via interpolation


restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/5
R = K[x,y,z]

-- ideal of a plane and a line that intersect at the origin
I = ideal (x*z,y*z)

-- make a black box from the ideal
bbI = blackBoxIdeal I;

point2 = matrix{{0,0,1_K}}
point1 = matrix{{1,2,0_K}}
point0 = matrix{{0,0,0_K}}
pointNothing = matrix{{1,1,1_K}}

assert (2== bbI.rankJacobianAt(point2))
assert (1== bbI.rankJacobianAt(point1))
assert (0== bbI.rankJacobianAt(point0))
-- this is a point where the ideal does not vanish.
assert (2==bbI.rankJacobianAt(pointNothing))


jet2 = jetAt(bbI,point2,20,2)
jet1a = jetAt(bbI,point1,20,2)
jet1b = jetAt(bbI,point1,20,2)
jet0 = jetAt(bbI,point0,20,2)
-- HashTable{lift => null      }
--           succeeded => false
--
-- richtig waere:
--
-- { failedJetLength => 2}
--               jet => | 0 eps -2eps |
--         succeeded => false

jet0 = new HashTable from {
     "failedJetLength" => 2,
     --"jet" => matrix{{0,eps,-2eps}},
     "succeeded" => false     
     }

-- test: are the rings of different jets of the same precision the same?
assert ((ring jet1a#"jet") === (ring jet2#"jet"))



                  jet => | 0 eps -2eps |
--         succeeded => false
class jet1

jetNothing = jetAt(bbI,pointNothing,20,2)
-- error:  function is not zero at the point

-- richtig waere:
--
-- { failedJetLength => 0}
--               jet => null
--         succeeded => false

-- es fehlt ein testbeispiel fuer failedJetLength => 1

viewHelp jetAt



---------------------

-- interpolation using jets
interpolate = (mons,jetP) -> if jetP#"succeeded" then (
	       s = syz sub(last coefficients sub(mons,jetP#"jet"),K);
	       I = ideal mingens ideal(mons*s);
	       I
	       )

mons1 = basis(1,R)

interpolate(mons1,jet2)
interpolate(mons1,jet1a)
interpolate(mons1,jet0)
-- ueberpruefen was hier passiert / passieren soll
-- es gibt kein sinnvolles ergebnis. Fehlermeldung.

interpolate = (mons,jetList) -> (
     jetsSucceeded = select(jetList,jetP->jetP#"succeeded");
     if #jetsSucceeded != #jetList then error "the point is not smooth";
     -- substitute jets into the monomials and take coefficients
     coeffList = apply(jetList,jetP ->  sub(last coefficients sub(mons,jetP#"jet"),K));
     -- find interpolation solution
     s = syz fold((a,b)->(a||b),coeffList);
     -- make polynomials from the solution
     I = ideal mingens ideal(mons*s);
     --I = ideal(mons*s);
     I
     )

interpolate(mons1,{jet1a,jet1b})
interpolate(mons1,{jet2})

interpolateMaxDeg = (maxDegree,R,jetList) -> (
     -- monomials in R of degree at most maxDegree
     mons := matrix {flatten apply(maxDegree+1,i->flatten entries basis(i,R))};
     interpolate(mons,jetList)
     )

interpolateMaxDeg(1,R,{jet1a,jet1b})
-- ideal(z)
interpolateMaxDeg(2,R,{jet1a,jet1b})
-- ideal(z)

interpolateBB = (maxDegree,BB,point) -> (
     R := BB.ring;
     mons := matrix {flatten apply(maxDegree+1,i->flatten entries basis(i,R))};
     -- find one jet with precision 10 more then number of monomials
     jetP := jetAt(BB,point,rank source mons+10,2);
     -- !!!this heuristic must be tested!!!
     -- Test: see if interpolated polynomials are in at least one
     -- irreducible component of the BlackBoxIdeal.
     interpolate(mons,{jetP})
     )

component1 = interpolateBB(2,bbI,point1)     
component2 = interpolateBB(2,bbI,point2)     

-- are we finished?
component1 + bbI.ideal == component1
-- yes
component2 + bbI.ideal == component2
-- yes     


--------------------
-- second example --
--------------------

bb2 = blackBoxIdeal (intersect(ideal(x,y),ideal(x^2+y^2+z^2)))

e = new Experiment from bb2;
e.run(100)

point1 = (e.pointsByKey({1}))#0
point2 = (e.pointsByKey({2}))#0

interpolateBB(1,bb2,point1)
-- ideal()
component1 = interpolateBB(2,bb2,point1)
-- Ideal(x^2+y^2+z^2)
component2 = interpolateBB(2,bb2,point2)

--
-- test if a point lies on a component defined by an iterpolation ideal
--

isOnComponent = (BB,interpolationIdeal,point,prec) -> (
     jetP = jetAt(BB,point,prec,1);
     if jetP#"succeeded" then (
     	  0==sub(interpolationIdeal,jetP#"jet")
	  ) else error "point not smooth"
     )
-- !! Possibly there is a usefull answer even if the point is not smooth !!


isOnComponent(bb2,component1,point1,10)
isOnComponent(bb2,component1,point2,10)
isOnComponent(bb2,component2,point2,10)

-- make a more complicated example where we do not have all 
-- generators in the interpolation ideal


-------------------
-- third example --
-------------------

bb3 = blackBoxIdeal (intersect(ideal(x,y^2+z^3+1),ideal(x^2+y^2+z^2)));

e = new Experiment from bb3;
e.run(100)

point1 = (e.pointsByKey({1}))#0
point2 = (e.pointsByKey({2}))#0

interpolateBB(1,bb3,point1)
-- ideal()
component1 = interpolateBB(2,bb3,point1)
-- Ideal(x^2+y^2+z^2)
component2 = interpolateBB(2,bb3,point2)
-- ideal(x)
-- this is not the complete ideal defining the component containing point2
interpolateBB(3,bb3,point2)
-- ideal(x,z^3+y^2+1)
-- this would be the complete ideal defining the component

-- for demonstration purposes we work with the incomplete ideal
component2
-- ideal(x)

-- there are four types of points:
--  a) those on the curve defined by ideal(x,z^3+y^2+1)
--  b) those on the surface defined by ideal(x^2+y^2+z^2)
--  c) those on the interesection of curve and surface
--  d) those on the surface which are not on the curve
--     but nevertheless satisfying x=0
--
-- the points of type d) should be classified to lie on the surface
-- even though they also lie in the partial ideal that contains the curve

point3 = (select(e.pointsByKey({1}),P->0==sub(component2,P)))#0

-- point3 lies on both partial ideals
assert (0==sub(component1,point3))
assert (0==sub(component2,point3))

-- but jets starting at point3 usually only lie on component1
assert isOnComponent(bb3,component1,point3,10)
assert not isOnComponent(bb3,component2,point3,10)
-- so using the Blackbox to generate jets we can check
-- if a point is on a component even if we have not interpolated
-- the full ideal. Generically this works as soon as we have a
-- single Hypersurface that contains only one component.
