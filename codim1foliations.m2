-- analyse the modulispace of
-- degree d foliations in AA^3
-- using finite field experiments

restart
K = ZZ/3
--R = K[x,y,z,w]
R = K[x,y,z]
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

-- stratification by rank of M
time tally apply(100,i->(
	  point = random(K^1,K^(rank source vars B));
	  rank sub(M,point)
	  ))
-- 54 => 1   
-- 59 => 2
-- 60 => 9997

-- stratification by rank of M
rank29points = flatten apply(100,i->(
	  point = random(K^1,K^(rank source vars B));
	  if 29==rank sub(M,point) then {point} else {}
	  ))

Mat = (point) -> sub(M,point)
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

point = rank29points#1
wA = omegaAat(point)
wB = omegaBat(point)
assert (wA*wB==0)
assert (1==rank source mingens ideal(differentialD(wA),wB))
assert (0==differentialD(wB))
differentialD(wB)
differentialD(wA)
-- sometimes A is closed

-- matrix drops rank generically only by one

-- consider only closed 2-forms
Iclosed = ideal sub(contract(super basis({0,0,d-2,3},ABRD),differentialD(omegaB)),B)
basisClosed = sub(transpose syz transpose jacobian Iclosed,K)
closedPointBat = (point) -> point*basisClosed
closedOmegaBat = (point) -> omegaBat(colestPointBat(point))

omegaAat = (point) -> 
time tally apply(10000,i->(
	  pointClosed = random(K^1,K^numVars);
	  --assert (0==differentialD omegaBat(pointClosed));
	  rank sub(M,pointBat(pointClosed))
	  ))


numVars = rank source vars B - rank source mingens Iclosed

time tally apply(10000,i->(
	  pointClosed = random(K^1,K^numVars);
	  --assert (0==differentialD omegaBat(pointClosed));
	  rank sub(M,pointBat(pointClosed))
	  ))

	  
-- stratification by betti tableau
betti (coeffB = matrix entries contract(
	  super basis({0,0,0,2},ABRD)
	  ,omegaB))
isHomogeneous coeffB

time tally apply(1000,i->(
	  point = random(K^1,K^(rank source vars B));
	  coeffBpoint = matrix entries sub(sub(coeffB,sub(vars A,ABRD)|point|sub(vars RD,ABRD)),RD);
	  assert isHomogeneous coeffBpoint;
	  (rank sub(M,point),betti syz coeffBpoint)
	  ))

-- make an experiment that can do this
-- and estimate codimension and # of components in each stratum
-- (without the use of rank jacobi)

-- make an approximate black box for the ideal
-- where the matrix drops rank
betti M
randomMinors = apply(10,i->random(K^(rank source M),K^(rank target M)));
-- should only fail in one of
(char K)^10
-- cases

evaluateRandomMinorsAt = (point) -> (
     MatPoint = sub(M,point)
     valuesAtPoint = matrix{apply(randomMinors,i->det(MatPoint*i))};
     if valuesAtPoint == 0 then assert (rank(MatPoint) < max(rank source MatPoint,rank target MatPoint));
     (valuesAtPoint == 0,rank MatPoint)
     )

tally apply(1000,i->(
	  evaluateRandomMinorsAt(
     	       random(K^1,K^(rank source vars B))
     	       )
	  ))

-- use blackBox and experiment here to find codimension of tangent spaces
-- using epsilon generated jacobi matrix
-- look at syzygies of found points