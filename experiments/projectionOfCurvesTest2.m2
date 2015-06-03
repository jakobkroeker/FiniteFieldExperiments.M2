-- test projection of a rational normal curve

restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

prime = 11
K = ZZ/prime
R = K[x,y,z,w]
betti res (I = minors(2,random(R^{2:0},R^{3:-1})))
projI = eliminate(w,I)

bb = blackBoxIdeal I;
e = new Experiment from bb;
time e.run(300)
-- used 0.150134 seconds
e.count()

-- project into IP^2, forgetting the forth variable
S = K[a,b,c]
mapS = matrix{{x,y,z}} 
f = map(R,S,mapS)  -- this is a map of rings !!
-- we need a map of spaces. This is defined by mapS:
sub(mapS,point)

-- a point on the curve
point = (e.pointsByKey({2}))#0
sub(mapS,point)
jet = jetAt(bb,point,20,1)    

imageJet = new HashTable from {
     "failedLength" => null,
     "jet" => sub(mapS,jet#"jet"),
     "succeeded" => true
     }
     
     
mons3 = basis(3,S)
interpolProjI = interpolate(mons3,{imageJet})
assert (interpolProjI == sub(projI,matrix{{a,b,c,0}}))

-- ToDo
-- define maps
-- apply maps to jets & points
-- interpolation with map != identity
-- ring of image either given or created



-- two rational normal curves
restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

prime = 11
K = ZZ/prime
R = K[x,y,z,w]
betti res (rationalNormalCurve = minors(2,random(R^{2:0},R^{3:-1})))
betti res (line = ideal random(R^{0},R^{2:-1}))
betti res (I = intersect(rationalNormalCurve,line))

projI = eliminate(w,I)
decompose projI

bb = blackBoxIdeal I;
e = new Experiment from bb;
time e.run(300)
-- used 0.150134 seconds
e.count()

-- a point on the curve
point = (e.pointsByKey({2}))#0

-- project into IP^2, forgetting the forth variable
S = K[a,b,c]
mapS = matrix{{x,y,z}} 
f = map(R,S,mapS)  -- this is a map of rings !!
-- we need a map of spaces. This is defined by mapS:
sub(mapS,point)

-- interpolate for all points
-- monomials up to degree 3
mons3 = basis(0,S)|basis(1,S)|basis(2,S)|basis(3,S)

apply(e.pointsByKey({2}),point->(
	  jet = jetAt(bb,point,20,1);
	  -- calculate the image of this jet
	  imageJet = new HashTable from {
     	       "failedLength" => null,
     	       "jet" => sub(mapS,jet#"jet"),
     	       "succeeded" => true
     	       };
	  ideal mingens interpolate(mons3,{imageJet})
     ))


assert (interpolProjI == sub(projI,matrix{{a,b,c,0}}))

