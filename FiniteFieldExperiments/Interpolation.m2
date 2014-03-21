



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

doc ///
   Key
        interpolate
   Headline
        find polynomials containing a list of jets
   Usage   
        I = interpolate(mons,jetList)
   Inputs  
        mons:Matrix 
            of monomials
	jetList:List
	    of jets    
   Description
        Text
	   Finds those linear combinations of monomials that vanish
	   on the given list of jets.
	   
	   Lets consider a black box that describes
	   a line and a plane intersecting at the origin:
        Example      
           K = ZZ/5
           R = K[x,y,z]
           I = ideal (x*z,y*z)
           bb = blackBoxIdeal I;       
        Text
           \break 
	   Consider a point on the line:
        Example
	   point = matrix{{0,0,1_K}}
	Text
	   \break
	   and a jet of lenght 3 starting at this point and
	   lying on the variety described by the black box ideal
	Example
	   j = jetAt(bb,point,4,1)     
	Text
	   \break
	   Now find linear polynomials containing this jet:
	Example
	   interpolate(matrix{{x,y,z}},{j})   
	Text
	   Notice that polynomials containig the line are found.
	   The surface is invisible to the interpolation.   
   Caveat
      This function does not estimate the lenght of the jet needed to
      get a useful answer. (The jet should be at least as long as the
      number of monomials). This is done by @TO interpolateBB @. 
   SeeAlso
      interpolateBB
      isOnComponent
      createInterpolatedIdeal
      createAllInterpolatedIdeals
      interpolatedIdealKeys      
///

TEST ///
  K = ZZ/5
  R = K[x,y,z]
  I = ideal (x*z,y*z)
  bb = blackBoxIdeal I;     
  -- a point on the line  
  point = matrix{{0,0,1_K}}
  j = jetAt(bb,point,4,1)
  assert (ideal(x,y) == interpolate(matrix{{x,y,z}},{j}))
///

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

doc ///
   Key
        interpolateBB
   Headline
        find polynomials containing a list of jets
   Usage   
        I = interpolateBB(maxDegree,BlackBox,point)
   Inputs  
        maxDegree:ZZ 
            the maximal degree of polynomials considered
	BlackBox:BlackBoxIdeal
	point: Matrix
	    a point where the Blackbox vanishes    
   Description
        Text
	   Finds all polynomials of degree at most maxDegree
	   that contain the component on which the point lies.
	   If the point is not smooth, an error will be produced.
	   
	   Lets consider a black box that describes
	   a line and a plane intersecting at the origin:
        Example      
           K = ZZ/5
           R = K[x,y,z]
           I = ideal (x*z,y*z)
           bb = blackBoxIdeal I;       
        Text
           \break 
	   Consider two points on the variety described 
	   by the blackbox:
        Example
	   pointOnLine = matrix{{0,0,1_K}}
	   pointOnPlane = matrix{{0,1,0_K}}
	Text
	   \break
	   Now find linear equations containing the respective
	   components on which the points lie:
	Example
	   interpolateBB(1,bb,pointOnLine)
	   interpolateBB(1,bb,pointOnPlane)
	Text
	   \break
	   Finding points on the different components can be done
	   by running an @TO Experiment @. Interpolating a component
	   for all points found can be done by @TO createAllInterpolatedIdeals @.   
  Caveat
      This function does not work with multigraded rings.
      At the moment this has to be done by hand with @TO interpolate @. 
      
      At the moment the interpolation is done by producing one
      jet of the appropriate length. Often one could interpolate
      much faster if several shorter jets were used. (Most of the
      time is used when producing the jets)
  SeeAlso
      interpolate
      isOnComponent
      createAllInterpolatedIdeals
      createInterpolatedIdeal
      interpolatedIdealKeys      
///

TEST ///
  K = ZZ/5
  R = K[x,y,z]
  I = ideal (x*(x^2+y^2+z^2+1),y*(x^2+y^2+z^2+1))
  bb = blackBoxIdeal I;     
  -- a point on the line  
  pointOnLine = matrix{{0,0,1_K}}
  pointOnSurface = matrix{{0,2,0_K}}
  intersectionPoint = matrix{{0,0,2_K}}
  notApoint = matrix{{1,1,1_K}}
  pointQQ = matrix{{0,0,1_QQ}}
  assert (ideal(x,y) == interpolateBB(1,bb,pointOnLine))
  -- no polynomials of degree 1 containing the surface
  assert (ideal(0_R) == interpolateBB(1,bb,pointOnSurface))
  assert (ideal(x^2+y^2+z^2+1) == interpolateBB(2,bb,pointOnSurface))
  -- a point where the blackBox does not vanish
  assert (try interpolateBB(1,bb,notApoint) then false else true)
  -- a point where the blackBox is not smooth
  assert (try interpolateBB(1,bb,intersectionPoint) then false else true)
  -- a point over the wrong field
  assert (try interpolateBB(1,bb,pointQQ) then false else true) 
///


isOnComponent = method();
-- checks if a point lies on the component of BB defined by an interpolation ideal which might contain only some 
-- of the polynomials defining the component.
isOnComponent(HashTable,HashTable,Matrix,ZZ) := Boolean => (BB,interpolationIdeal,point,prec) -> (
     jetP := jetAt(BB,point,prec,1);
     if jetP#"succeeded" then (
     	  0==sub(interpolationIdeal,jetP#"jet")
	  ) else error "point not smooth"
     )
-- !! Possibly there is a usefull answer even if the point is not smooth !!



FFEInterpolationData = new Type of MutableHashTable;

createFFEInterpolationData = ( mapfkt) -> (
     iData := new FFEInterpolationData;
     iData.interpolatedIdeals = new MutableHashTable; 
     return iData;
)


createInterpolatedIdeal (ZZ,BlackBoxIdeal,Matrix) := InterpolatedIdeal => (maxDegree, BB, point)->
(
   createInterpolatedIdeal ( interpolateBB(maxDegree,BB,point), maxDegree, "" )
)


FFEInterpolation = new Type of HashTable;

createFFEInterpolation = method();


-- observer observable längst fällig!

createFFEInterpolation(Experiment,Ring, HashTable) := HashTable => (experiment,imageRing, mapdata)->
(
    interpolation := new MutableHashTable;

    interpolation.experiment = () -> experiment;

    interpolation.blackBoxIdeal = () -> (interpolation.experiment()).blackBoxIdeal();
    
    interpolationData := createFFEInterpolationData(mapdata#"matrix");

     -- duplicate code..(otherwise too many parameters...)
     interpolation.isOnComponent = method();
     interpolation.isOnComponent (HashTable,Matrix,ZZ) := Boolean =>  (interpolationIdeal,point,prec) -> (
         jetP := jetAt(interpolation.blackBoxIdeal() ,point,prec,1);
         if jetP#"succeeded" then (
         	  0 == sub( interpolationIdeal, mapdata#"mapfkt"(jetP#"jet") )
	      ) else error "point not smooth"
     );

    interpolation.createAllInterpolatedIdeals  = (maxDegree, prec) -> 
    ( 
        interpolatedIdeals := {};
        for point in interpolation.experiment.points() do
        ( 
            if (interpolation.blackBoxIdeal()).isSingular(point) then continue;
             -- check if point is already on one of the known components
             bIsOnComponent := false;
             for interpolData in interpolatedIdeals do
             (
                 if  interpolation.isOnComponent (  interpolData.ideal, point, 0 ) then
                 (
                     if  interpolation.isOnComponent (  interpolData.ideal, point, prec ) then
                     (
                          bIsOnComponent = true; 
                          break;
                     );
                 );
             );
             if bIsOnComponent then continue;

             interpolatedIdeals = interpolatedIdeals | { createInterpolatedIdeal (maxDegree, interpolation.blackBoxIdeal(), point)  };             
        );
        interpolationData.interpolatedIdeals = new MutableHashTable from 
           apply ( #interpolatedIdeals, idx-> (("ideal_" |toString idx ) => interpolatedIdeals#idx ) ) ;
       
    );

    localMembershipPrecision := 10; -- does this belong to experimentData?


    interpolation.setMembershipPrecision = (prec)->
    (
         localMembershipPrecision = prec;
    );

    interpolation.membershipPrecision = (prec)->
    (
         localMembershipPrecision 
    );

    interpolation.interpolatedIdealKeys = method();
    interpolation.interpolatedIdealKeys (Matrix,ZZ) := Thing => (point,prec)->
    (
       numbers := {};
 
        if (interpolation.blackBoxIdeal()).isSingular(point) then return "is not smooth";    

        for key in keys interpolationData.interpolatedIdeals do
        (
          if  interpolation.isOnComponent (  (interpolationData.interpolatedIdeals#key).ideal, point, 0,  ) then
          (
              if  interpolation.isOnComponent (  (interpolationData.interpolatedIdeals#key).ideal, point, prec,  ) then
              (
                 numbers = numbers | {key};
              );
           );
        );
        return numbers;
 
    );

   interpolation.interpolatedIdealKeys (Matrix) := Thing => (point)->
   (
        return  interpolation.interpolatedIdealKeys(point, localMembershipPrecision );
   );

   iik := (point)-> (return  interpolation.interpolatedIdealKeys(point) ;);

   --(interpolation.blackBoxIdeal()).rpp("interpolatedIdealKeys", iik);



    interpolation.printInterpolatedIdeals = ()->
   (
         apply (keys interpolationData.interpolatedIdeals, key-> ( print (key=>new HashTable from interpolationData.interpolatedIdeals#key ) ) );
   );

    interpolation = newClass( FFEInterpolationData, interpolation );
    return interpolation;
     
);

createFFEInterpolation(Experiment,Ring, Matrix) := HashTable => (experiment,imageRing, pmap)->
(
    mapData := new MutableHashTable;
    mapData#"mapfkt" =  method();
    mapData#"matrix" = pmap;
    rng := (experiment.blackBoxIdeal()).ring ; 
    mapData#"map" = map( rng, imageRing, pmap );
    
    mapData#"mapfkt" (Matrix) := Matrix => (point)->
    (
        return sub(pmap,point);
    );
    
    return createFFEInterpolation(experiment,imageRing,mapData);
    
);

createFFEInterpolation(Experiment) := HashTable => (experiment)->
(
    mapData := new MutableHashTable;
    mapData#"mapfkt" =  method();
    rng := (experiment.blackBoxIdeal()).ring ;
    mapData#"matrix" = matrix{ gens rng };
    mapData#"map" = map( rng, rng, mapData#"matrix" );
    
    mapData#"mapfkt" (Matrix) := Matrix => (point)->
    (
        return point;
    );
    
    return createFFEInterpolation(experiment, rng, mapData);
);

-- todo : new type for mapData. Also support dot access to mapData members.


end
----

restart

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

viewHelp interpolateBB

