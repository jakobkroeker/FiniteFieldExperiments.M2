-- construct (hopefully new) codim 1-foliations
-- from curves in IP^3

restart
K = ZZ/101
R = K[x,y,z,w]
D = K[dx,dy,dz,SkewCommutative=>true]
RD = R**D

d = 4
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

betti (I = ideal mingens(Iint+Idiff))

-- curve part
indC = flatten apply(flatten entries super basis({0,2},RD),i->(
	       apply(rank source super basis({d-1,0},RD),j->{i,j})));

C = K[apply(indC,i->c_i)]
CRD = C**RD
use R
curve = ideal(x,y)
curve = ideal random(R^{0},R^{2:-2})
curve = ideal random(R^{0},R^{-1,-2})
betti res (curve = minors(2,random(R^{2:0},R^{ -1,-1,-2})))

link33 = (curve) -> (
     (ideal (random(3,curve),random(3,curve))):curve
     )
betti res (curve = link33(curve))

betti res (coeffB = ideal apply(3,i->random(3,curve)))

rank source super basis(7,ker gens coeffB)
rank source super basis(3,R)

betti (gensCurve = super basis(d-1,curve))
domegaC = sub(gensCurve,CRD) * transpose matrix {apply(rank source gensCurve,i->sum apply(flatten entries super basis({0,2},RD),j->sub(c_{j,i},CRD)*sub(j,CRD)))};


-- conditions for d(domegaC) = 0 closed
betti (IclosedC = ideal mingens ideal contract(super basis({0,d-2,3},CRD),differentialD(domegaC)))
pointB = sub(random(K^1,K^(#indC-rank source gens IclosedC))*transpose syz transpose jacobian(sub(IclosedC,C)),K)

domega = sub(domegaC,sub(pointB,RD)|vars RD)
assert (0==differentialD(domega))

betti ideal mingens sub(Idiff,sub(vars A,AB)|pointB)
betti ideal mingens sub(Iint,sub(vars A,AB)|pointB)
betti res (Icoeff = sub(ideal contract(super basis ({0,2},RD),domega),R))
betti super basis(7,ker mingens Icoeff)
betti res sub(ideal contract(super basis ({0,2},RD),domega),R)

