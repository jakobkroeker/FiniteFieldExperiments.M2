-- test decomposition of cyclic-n

-- references:

-- 1) Numerical-symbolic exact irreducible decomposition of cyclic-12
-- Rostam Sabeti
-- LMS J. Comput. Math. 14 (2011) 155–172 Ce2011 Author doi:10.1112/S146115701000001X

-- 2) http://www-sop.inria.fr/saga/POL/BASE/2.multipol/cyclic.html

restart
loadPackage"BlackBoxIdeals"

n = 8
K = ZZ/5
R = K[x_1..x_n]
cycMonomial = (start,len) -> (
     product apply(len,i->(
	       j = start+i;
	       if j>n then j=j-n;
	       x_j
	       )
	  ))

cyc = (len) -> sum apply(n,i->cycMonomial(i+1,len))

cycIdeal = ideal (apply(n-1,i->cyc(i+1))) + ideal((product gens R)-1)

--codim cycIdeal
--degree cycIdeal

--time decompose cycIdeal;
-- used 0.071498 seconds (n=4, with 1)
-- used 0.128491 seconds (n=5, only homogeneous part)

cycBB = blackBoxIdeal cycIdeal

loadPackage"FiniteFieldExperiments"

--apropos "Experiment"
--viewHelp"FiniteFieldExperiments"

cycE = new Experiment from cycBB
time cycE.run(100000)
-- used 4.69959 seconds
apply(unique cycE.points(),P->rank cycBB.jacobianAt(P))

--dCyc = decompose cycIdeal
--tally apply(cycE.points(),P->apply(dCyc,component->sub(component,P)))
--apply(dCyc,i->codim i)

maxDeg = 2
mons = matrix {flatten apply(maxDeg+1,i->flatten entries super(basis(i,R)))}
-- alle monome bis grad maxDeg
-- monomialsToMaxDegree(R,maxDeg)


-- mons: matrix von m=Monomen
-- jetP: jet in backBoxIdeal
-- Ausgabe
-- I: ideal der Polynome mit monomen aus mons die auf jetP verschwinden

interpolate(mons,jetP)
     
L = flatten apply(unique cycE.points(),P->(
	  rank cycBB.jacobianAt(P);
	  time jetP = JetAt(cycBB,P,100,1);
	  -- used 116.95 seconds (200, n=8, p=5)
	  if jetP#"succeeded" then (
	       s = syz sub(last coefficients sub(mons,jetP#"jet"),K);
	       I = ideal mingens ideal(mons*s);
	       {I}
	       ) else {}
     ))


