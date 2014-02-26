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

bb2 = blackBoxIdeal (intersect(ideal(x,y),ideal(x^2+y^2+z^2)));

e = new Experiment from bb2;
e.run(100)

point1 = (e.pointsByKey({1}))#0
point2 = (e.pointsByKey({2}))#0

interpolateBB(1,bb2,point1)
-- ideal()
component1 = interpolateBB(2,bb2,point1)
-- ideal(x^2+y^2+z^2)
component2 = interpolateBB(2,bb2,point2)
-- ideal (y, x)

assert isOnComponent(bb2,component1,point1,10)
assert not isOnComponent(bb2,component1,point2,10)
assert isOnComponent(bb2,component2,point2,10)

-- are we finished?
component1 + bb2.ideal == component1
-- yes
component2 + bb2.ideal == component2
-- yes     

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
-- ideal(x^2+y^2+z^2)
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
--  a) those on the curve defined by ideal(x,z^3+y^2+1)   -- {2}
--  b) those on the surface defined by ideal(x^2+y^2+z^2) -- {1}
--  c) those on the interesection of curve and surface    -- {0}
--  d) those on the surface which are not on the curve    -- {1}
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
-- the full ideal. Generically this works as soon as we have found a
-- single Hypersurface that contains only one component.

-- are we finished?
component1 + bb2.ideal == component1
-- yes
component2 + bb2.ideal == component2
-- no

