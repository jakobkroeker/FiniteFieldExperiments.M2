-- analyse the modulispace of
-- degree d foliations in AA^3
-- using finite field experiments

restart

needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/3
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

-- eliminate A
betti (M = matrix entries sub(contract(sub(transpose vars A,AB),mingens I),B))
-- 0: 18 46
assert isHomogeneous M

Mat = (point) -> sub(M,point)
rankMat = (point) -> rank Mat(point)
rankMatEx = (bb,point)-> rank Mat(point)
end
---

restart
uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

---

restart
load"experiments/codim1foliations.m2"

bbRankM = blackBoxParameterSpace(#(gens B),K)
bbRankM = bbRankM.registerPointProperty("rankJacobianAt",rankMat)
keys bbRankM
pointProperties bbRankM
-- {rankJacobianAt}

bbRankM.rankJacobianAt --ok

e = new Experiment from bbRankM
keys e
e.setPointsPerComponent(20);
e.pointsPerComponent() --20
e.watchProperties {"rankJacobianAt"};
watchedProperties e -- {rankJAcobianAt}

time e.run 400
-- 1.63 seconds for 400
pointPropertyVals = realizedPointProperties e
--points = e.pointsByKey({54}) 

points = e.pointsByKey(first pointPropertyVals)
point = points#0
-- howto apply properties:
bbRankM.rankJacobianAt point 
bbRankM#"rankJacobianAt" point -- safer

counts e
sort counts e

-- ?
estimateStratification e
estimateDecomposition e
-- e.stratificationIntervalView() -- (jk) does not work any more, because stratification() method was modified... 
e.collectedCount()
e.jacobianAtKey() --null
-- ??
keys e

--tally apply( e.pointKeys(), key ->( #((e.pointLists())#key)=>key))
-- lieber e.points().
   
estimateStratification e



omegaBat = (point) -> sub(sub(omegaB,sub(vars A,ABRD)|point|sub(vars RD,ABRD)),RD)
randomPointAat = (point) -> ( 
     allPointsA = transpose syz transpose Mat(point);
     if rank target allPointsA == 0 then return null;
     if rank target allPointsA == 1 then (
     	  randomPointA = allPointsA
     	  ) else (
	  randomPointA = (random(K^1,K^(rank target allPointsA))*allPointsA);
	  );
     randomPointA
     )

omegaAat = (point) -> (     
     sub(sub(omegaA,randomPointAat(point)|sub(vars B,ABRD)|sub(vars RD,ABRD)),RD)
     )


wA = omegaAat(point)
wB = omegaBat(point)
assert (wA*wB==0);
assert (1==rank source mingens ideal(differentialD(wA),wB));
assert (0==differentialD(wB));
-- no
assert (0==differentialD(wA))
-- yes

-- sometimes A is closed

-- matrix drops rank generically only by one

-- stratification by betti tableau
betti (coeffB = matrix entries contract(
	  super basis({0,0,0,2},ABRD)
	  ,omegaB))
isHomogeneous coeffB
coeffBat = (point) -> sub(sub(coeffB,sub(vars A,ABRD)|point|sub(vars RD,ABRD)),R)
bettiAt = (point) -> betti res ideal coeffBat(point)
--bettiAtEx = (bb,point)->bettiAt(point)
--bbBetti = blackBoxIdealFromProperties(#(gens B),K,i->(rankMat i, bettiAt i))
bbBetti = blackBoxParameterSpace(#(gens B),K)
bbBetti.registerPointProperty bettiAt
bettiInteresting = (point) -> (rankMat(point) < #(gens A))

eBetti = new Experiment from bbBetti
eBetti.setIsInteresting(bettiInteresting)
eBetti.watchProperties {"bettiAt"}
time eBetti.run 10000
-- used 19.2529 seconds (3 variables)
-- used 40.5503 seconds (4 variables)
estimateStratification2(eBetti)
-- sollte nach intervallen sortiert sein
-- eventuell nur keys mit mindestens 10 punkten?
-- eventuell nur keys mit codim < vorgabe?
-- uebersichtlicher: ohne intervalle

-- try Property
isAclosedAt = (point) -> (
     wA := omegaAat(point);
     if wA===null then null else 0==differentialD wA
     )

bbBetti = bbBetti.rpp isAclosedAt

eBetti.tryProperty "isAclosedAt"

-- consider only closed 2-forms
Iclosed = ideal sub(contract(super basis({0,0,d-2,3},ABRD),differentialD(omegaB)),B)
basisClosed = sub(transpose syz transpose jacobian Iclosed,K)
closedPointBat = (point) -> point*basisClosed
closedOmegaBat = (point) -> omegaBat(closedPointBat(point))
closedBettiAt = (point) -> betti res ideal coeffBat(closedPointBat(point))
closedMat = (point) -> Mat(closedPointBat(point))
closedRankMat = (point) -> rank closedMat(point)
bbClosed = bbClosed.rpp closedRankMat;


bbClosed = blackBoxParameterSpace(rank target basisClosed,K)
bbClosed.numVariables
bbClosed.registerPointProperty closedBettiAt
eClosedBetti = new Experiment from bbClosed
eClosedBetti.setIsInteresting ( (point) -> (closedRankMat(point) < #(gens A)))
eClosedBetti.watchProperties {"closedBettiAt"}
time eClosedBetti.run(10000)
-- used 42.6086 seconds
estimateStratification(eClosedBetti)
--             {total: 1 2 1} => 1  }
--                  0: 1 . .
--                  1: . 2 .
--                  2: . . 1
--                     0 1 2 3
--             {total: 1 3 3 1} => 1
--                  0: 1 . . .
--                  1: . 3 . .
--                  2: . . 3 .
--                  3: . . . 1
--                     0 1 2 3
--             {total: 1 3 3 1} => 1
--                  0: 1 . . .
--                  1: . 3 1 .
--                  2: . . 2 1
--closedOmegaAat = (point) -> omegaAat(closedPointBat(point))
testPoint = (eClosedBetti.pointsByKey((keys eClosedBetti.pointLists())#0))#0

closedIsAclosedAt = (point) -> isAclosedAt(closedPointBat(point))
bbClosed = bbClosed.rpp closedIsAclosedAt

watchedProperties eClosedBetti
time eClosedBetti.run 10000 
-- used 424.527 seconds
estimateStratification eClosedBetti

eClosedBetti.tryProperty closedIsAclosedAt
eClosedBetti.tryProperty closedRankMat

-------------------------
-- reworked up to here --
-------------------------

time betti (J = jacobian I)

closedRankJacobiAt = (point) -> (
     pointB := closedPointBat(point);
     pointA := randomPointAat(pointB);
     tangentB := syz (transpose syz transpose sub(J,pointA|pointB))_{rank source pointA..rank source pointA+rank source pointB-1};
     rank tangentB + rank(basisClosed) - rank (tangentB | transpose basisClosed) 
     )

time eClosedBetti.run(100000) 
-- used 556.586 seconds
estimateStratification2(eClosedBetti)
tryProperty(eClosedBetti,i->(closedIsAclosedAt i,closedRankJacobiAt i))

-- jets in AB-space

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
eBetti.watchProperty "bettiAt"
eBetti.setIsInteresting ( (i)->true)  -- default

time eBetti.run 1000
-- used 9.11658 seconds on my old ThinkPad
time eBetti.run 10000

time eBetti.run 100000
--eBetti.estimateStratification()
estimateStratification2 eBetti
