-- analyse the modulispace of
-- degree d foliations in AA^3
-- using finite field experiments
--
-- try only closed differential forms as domega

restart

needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/101
R = K[x,y,z,w]
--R = K[x,y,z]
D = K[dx,dy,dz,SkewCommutative=>true]
RD = R**D

d = 3
indA = flatten entries super basis({d,1},RD); #indA
A = K[apply(indA,i->a_i)]
indB = flatten entries super basis({d-1,2},RD); #indB
B = K[apply(indB,i->b_i)]
AB = A**B
ABRD = AB**RD

-- differentiate
differentialD = (w) -> (
     localRing = ring w;
     lx = (symbol x)_localRing;
     ly = (symbol y)_localRing;
     lz = (symbol z)_localRing;
     ldx = (symbol dx)_localRing;
     ldy = (symbol dy)_localRing;
     ldz = (symbol dz)_localRing;
     diff(lx,w)*ldx+diff(ly,w)*ldy+diff(lz,w)*ldz
     )

omegaA = sum apply(indA,i->a_i*sub(i,ABRD));
omegaB = sum apply(indB,i->b_i*sub(i,ABRD));

-- integrability condition w \wedge dw = 0
mons3 = super basis({2*d-1,3},RD)
betti (Iint = ideal sub(contract(sub(mons3,ABRD),omegaA*omegaB),AB))

-- diffential condition d(w_A) = w_B
betti (Idiff = sub(ideal mingens minors(2,contract( matrix{apply(indB,i->sub(i,ABRD))},
	  matrix {{differentialD(omegaA)},{omegaB}})),AB))

-- total conditions
betti (I = ideal mingens (Iint+Idiff))
assert isHomogeneous I

bb = blackBoxIdeal(I);

-- smart point search
-- choose coefficients B first, then A if possible

-- eliminate A
betti (M = matrix entries sub(contract(sub(transpose vars A,AB),mingens I),B))
-- 0: 18 46
assert isHomogeneous M

coeffDomegaA = contract(sub(super basis({d-1,2},RD),ABRD),matrix{{differentialD(omegaA)}})
MclosedA = contract(transpose vars A,sub(coeffDomegaA,A))

Mat = (pointB) -> sub(M,pointB)
rankMat = (pointB) -> rank Mat(pointB)

end
---

restart
uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

---

restart
load"experiments/codim1foliationsSmartSearch.m2"

-- split point in A part and B part
pointAat = (point) -> point_{0..(#gens A-1)}
pointBat = (point) -> point_{#gens A..#gens A + #gens B -1}
assert (pointAat(matrix{bb.unknowns})|pointBat(matrix{bb.unknowns}) == matrix{bb.unknowns})

-- differential Forms
omegaAat = (point) -> sub(omegaA,point|vars RD)
omegaBat = (point) -> sub(omegaB,point|vars RD)


-- consider only closed 2-forms for B
-- condition for closedness
Iclosed = ideal sub(contract(super basis({0,0,d-2,3},ABRD),differentialD(omegaB)),B)
-- basis for closed subspace
basisClosed = sub(transpose syz transpose jacobian Iclosed,K)


randomPoint = () -> (
     pointA := null;
     -- random closed 2 form for pointB
     pointB := random(K^1,K^(rank target basisClosed))*basisClosed;
     -- points A that satisfy omegaA * omegaB = 0 and D(omegaA) = l * omegaB
     allPointsA := transpose syz transpose Mat(pointB);
     -- if there are no such points A return null
     if rank target allPointsA == 0 then return null;
     -- if there is a one dimensional space of such points return a basis vector
     if rank target allPointsA == 1 then (
     	  pointA = allPointsA
     	  ) else (
	  -- otherwise a random linear combination of basis vectors
	  pointA = (random(K^1,K^(rank target allPointsA))*allPointsA);
	  );
     pointA|pointB
     )

-- stratification by betti tableau
betti (coeffB = matrix entries contract(super basis({0,0,0,2},ABRD),omegaB))
coeffBat = (point) -> sub(sub(coeffB,point|sub(vars RD,ABRD)),R)
bettiAt = (point) -> betti res ideal coeffBat(point)
bb.rpp("bettiAt",bettiAt);

-- check wether forms are closed
isAclosedAt = (point) -> (
     wA := omegaAat(point);
     if wA===null then null else 0==differentialD wA
     )
bb.rpp("isAclosedAt",isAclosedAt);

isBclosedAt = (point) -> (
     wB := omegaBat(point);
     if wB===null then null else 0==differentialD wB
     )
bb.rpp("isBclosedAt",isBclosedAt);

isAclosedAt2 = (point) -> (
     rank(Mat(point)|MclosedA) == rank Mat(point)
     )
-- more precise but slower


--- the experiment
e = new Experiment from bb;
-- on second trial error:
-- method interpolatedIdealKeys seems already registered, please use 'updatePointProperty' for updating.

-- find points by the method defined above
e.setPointGenerator(randomPoint)

-- ignore points for which A is closed
e.setIsInteresting ( (point) -> not isAclosedAt(point) )

e.watchProperty("bettiAt")
time e.run(10)
-- used 25. seconds (for 4)
e.trials()
e.watchedProperties()
e.tryProperty("isAclosedAt")
-- no A is closed (because of isInteresting)
e.tryProperty("isBclosedAt")
-- all B are closed
e.tryProperty("bettiAt")
e.tryProperty("isCertainlySingularAt")

-- sysygies of three cubics
restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/3
R = K[x,y,z,w]
d = 3
mons = super basis(d,R)
indB = flatten entries mons
B = K[flatten apply(3,i->apply(indB,j->b_{i,j}))]

betti (M = transpose matrix apply(3,i->apply(indB,j->b_{i,j})))
(M^{0})
M
Mat = (point) -> sub(M,point)
-- Blackbox
Iat = (point) -> ideal(mons*Mat(point))



numVariables = #(gens B)
coefficientRing B
--bb= blackBoxIdealFromEvaluation(numVariables,K,Iat) --fails.
bbBetti = blackBoxParameterSpace(#(gens B),K)

bbBetti.registerPointProperty("Mat" , (point) -> sub(M,point))
bbBetti = bbBetti.registerPointProperty("Iat" , (bb,point) -> ideal(mons*bb.Mat(point) ))

bbBetti = bbBetti.registerPointProperty("bettiAt" , (bb,point) -> betti res bb.Iat(point) )

--


--bbBetti = blackBoxParameterSpace(#(gens B),K)
--bbBetti.registerPointProperty("bettiAt",bettiAt)
eBetti = new Experiment from bbBetti
--eBetti.numVariables = #(gens B)
--eBetti.coefficientRing = K
--eBetti.property = bb.bettiAt
eBetti.watchProperty("bettiAt")
eBetti.setIsInteresting ( (i)->true)  -- default

time eBetti.run(1000)
-- used 9.11658 seconds on my old ThinkPad
time eBetti.run(10000)

time eBetti.run(100000)
--eBetti.estimateStratification()
estimateStratification2(eBetti)

-- Darboux examples
restart
load"experiments/codim1foliationsSmartSearch.m2"

lambdaPartition = {2,1,1,1}
idealFromPartition = (lambdaPartition) -> (
     lambdaPolynomials = apply(lambdaPartition,i->random({i,0},RD));
     w = sum apply(#lambdaPartition,i->
     	  random(K)*product apply(#lambdaPartition,j->(
	       	    if i==j then 
	       	    differentialD(lambdaPolynomials#j)
	       	    else
	       	    lambdaPolynomials#j
	       	    ))
     	  );
     dw = differentialD(w);
     assert (0==w*differentialD(w));
     (I = sub(ideal contract(super basis({0,2},RD),dw),R))
     )

tally apply(partitions 5,i->(
	  P = toList i;
	  I = idealFromPartition(P);
	  if 0==I then dI = {} else time dI = primaryDecomposition I;
	  curveComponents = select(dI,i->codim i==2);
	  if 0==#curveComponents then curve = ideal{1_R} else curve = intersect curveComponents;
	  pointComponents = select(dI,i->codim i==3);
	  if 0==#pointComponents then points = ideal{1_R} else points = intersect pointComponents;
	  irrelevant = select(dI,i->codim i==4);
	  assert (irrelevant=={});
	  print (P,degree I,betti res curve, #curveComponents,degree points,betti res I)
	  ))
