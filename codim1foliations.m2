-- analyse the modulispace of
-- degree d foliations in AA^3
-- using finite field experiments

needsPackage"BlackBoxIdeals"
load"FiniteFieldExperiments.m2"

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

-- eliminate B
betti (M = matrix entries sub(contract(sub(transpose vars A,AB),mingens I),B))
-- 0: 18 46
assert isHomogeneous M

Mat = (point) -> sub(M,point)
rankMat = (point) -> rank Mat(point)

end
---

restart
load"codim1foliations.m2"

bbRankM = blackBoxIdealFromProperties(#(gens B),K,rankMat)

e = new Experiment from bbRankM
e.run(1000)
e.estimateStratification()
apply(keys e.getPointData(),i->#((e.getPointData())#i))
-- lieber e.points().

sortByFrequency = (t) -> sort apply(keys t,key->(t#key,key))
estimateStratification2 = (e) -> (
     count := e.getCountData();
     trials := e.getTrials();
     charK := char K; -- this must be read from the experimentdata
     print "--";
     apply(sortByFrequency(count),i->(
	       --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
	       print (net(round(1,(log(trials)-log(i#0))/log(char K)))|" <= "|net(i#1)))
	       )
	  ;
     print "--";
     )

estimateStratification2(e)     

-- n/t = 1/p^c = p^-c
-- log n - log t = -c * log p
-- log t - log n =  c*log p
-- c = (log t - log n)/log p




omegaBat = (point) -> sub(sub(omegaB,sub(vars A,ABRD)|point|sub(vars RD,ABRD)),RD)
omegaAat = (point) -> ( 
     allPointsA = transpose syz transpose Mat(point);
     if rank target allPointsA == 0 then return null;
     if rank target allPointsA == 1 then (
     	  randomPointA = allPointsA
     	  ) else (
	  randomPointA = (random(K^1,K^(rank target allPointsA))*allPointsA);
	  );
     sub(sub(omegaA,randomPointA|sub(vars B,ABRD)|sub(vars RD,ABRD)),RD)
     )

-- wA = omegaAat(point);
-- wB = omegaBat(point);
-- assert (wA*wB==0);
-- assert (1==rank source mingens ideal(differentialD(wA),wB));
-- assert (0==differentialD(wB));
-- assert (0==differentialA(wB))

-- sometimes A is closed

-- matrix drops rank generically only by one

-- stratification by betti tableau
betti (coeffB = matrix entries contract(
	  super basis({0,0,0,2},ABRD)
	  ,omegaB))
isHomogeneous coeffB
coeffBat = (point) -> sub(sub(coeffB,sub(vars A,ABRD)|point|sub(vars RD,ABRD)),R)
bettiAt = (point) -> betti res ideal coeffBat(point)

--bbBetti = blackBoxIdealFromProperties(#(gens B),K,i->(rankMat i, bettiAt i))
bbBetti = blackBoxIdealFromProperties(#(gens B),K,i->(bettiAt i))
bbBetti.isZeroAt = (point) -> (rankMat(point) < #(gens A))

eBetti = new Experiment from bbBetti
time eBetti.run(10000)
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

tryProperty = (experiment,property) ->(
     points := experiment.getPointData();
     apply(
     	  apply((keys points),key -> (key,tally apply(points#key,property))),
     	  i->print (net i#0|" => "|net i#1)
     	  )
     )

tryProperty(eBetti,isAclosedAt)
time eBetti.run(10000)
 -- used 396.944 seconds
tryProperty(eBetti,isAclosedAt)

estimateStratification2(eBetti)

-- consider only closed 2-forms
Iclosed = ideal sub(contract(super basis({0,0,d-2,3},ABRD),differentialD(omegaB)),B)
basisClosed = sub(transpose syz transpose jacobian Iclosed,K)
closedPointBat = (point) -> point*basisClosed
closedOmegaBat = (point) -> omegaBat(closedPointBat(point))
closedBettiAt = (point) -> betti res ideal coeffBat(closedPointBat(point))
closedMat = (point) -> Mat(closedPointBat(point))
closedRankMat = (point) -> rank closedMat(point)


bbClosed = blackBoxIdealFromProperties(rank target basisClosed,K,i->(closedBettiAt i))
bbClosed.isZeroAt = (point) -> (closedRankMat(point) < #(gens A))
bbClosed.numVariables()

eClosedBetti = new Experiment from bbClosed
time eClosedBetti.run(10000)
-- used 42.6086 seconds
estimateStratification2(eClosedBetti)
--closedOmegaAat = (point) -> omegaAat(closedPointBat(point))
closedIsAclosedAt = (point) -> isAclosedAt(closedPointBat(point))

time eClosedBetti.run(10000) 
-- used 424.527 seconds
 
tryProperty(eClosedBetti,closedIsAclosedAt)
tryProperty(eClosedBetti,closedRankMat)
tryProperty(eClosedBetti,point->(closedRankMat(point),closedIsAclosedAt(point)))



-- make an approximate black box for the ideal
-- where the matrix drops rank
betti M
randomMinors = apply(12,i->random(K^(rank source M),K^(rank target M)));
-- should only fail in one of
(char K)^#randomMinors
-- cases

evaluateRandomMinorsAt = (point) -> (
     MatPoint = closedMat(point);
     valuesAtPoint = matrix{apply(randomMinors,i->det(MatPoint*i))};
     if valuesAtPoint == 0 then assert (rank(MatPoint) < max(rank source MatPoint,rank target MatPoint));
     --(valuesAtPoint == 0,rank MatPoint)
     valuesAtPoint
     )

evaluateRandomMinorsAt(random(K^1,K^(bbClosed.numVariables())))
tryProperty(eClosedBetti,evaluateRandomMinorsAt)

bbMinors = blackBoxIdealFromEvaluation(bbClosed.numVariables(),K,evaluateRandomMinorsAt)


keys bbClosed
points = eClosedBetti.getPointData()
point = points#((keys points)#0)#0
bbMinors.valuesAt(point)

M0K = closedMat(point)*randomMinors#0
Ke = K[eps]/eps^2
t = matrix {apply(bbClosed.numVariables(),i->if i==0 then 1 else 0_K)}
M0Keps = closedMat(
     point+eps*sub(t,Ke)
     )*randomMinors#0
M0eps = M0Keps-M0K
M0K


time bbMinors.jacobianAt(point)

time bbMinors.run(100)

eMinors = new Experiment from bbMinors

time eMinors.run(100)

tally apply(1000,i->(
	  evaluateRandomMinorsAt(
     	       random(K^1,K^(rank source vars B))
     	       )
	  ))

-- use blackBox and experiment here to find codimension of tangent spaces
-- using epsilon generated jacobi matrix
-- look at syzygies of found points