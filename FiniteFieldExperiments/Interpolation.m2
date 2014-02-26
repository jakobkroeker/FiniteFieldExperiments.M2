

-- find polynomials containing a component
-- via interpolation


-- find linear combinations of monomials containing a given jet
interpolate = (mons,jetList) -> (
     R := ring mons;
     K := coefficientRing R;
     jetsSucceeded := select(jetList,jetP->jetP#"succeeded");
     if #jetsSucceeded != #jetList then error "the point is not smooth";
     -- substitute jets into the monomials and take coefficients
     coeffList := apply(jetList,jetP ->  sub(last coefficients sub(mons,jetP#"jet"),K));
     -- find interpolation solution
     s := syz fold((a,b)->(a||b),coeffList);
     -- make polynomials from the solution
     I := ideal mingens ideal(mons*s);
     --I = ideal(mons*s);
     I
     )

-- find all polynomials of degree smallerEqual than maxDegree containing the component of BB containing the point
-- 
-- this is the most basic simple minded implementation where only one very long jet is considered.
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

-- checks if a point lies on the component of BB defined by an interpolation ideal which might contain only some 
-- of the polynomials defining the component.
isOnComponent = (BB,interpolationIdeal,point,prec) -> (
     jetP := jetAt(BB,point,prec,1);
     if jetP#"succeeded" then (
     	  0==sub(interpolationIdeal,jetP#"jet")
	  ) else error "point not smooth"
     )
-- !! Possibly there is a usefull answer even if the point is not smooth !!


createInterpolatedIdeal (ZZ,BlackBoxIdeal,Matrix) := InterpolatedIdeal => (maxDegree, BB, point)->
(
   createInterpolatedIdeal ( interpolateBB(maxDegree,BB,point), maxDegree, "" )
)




