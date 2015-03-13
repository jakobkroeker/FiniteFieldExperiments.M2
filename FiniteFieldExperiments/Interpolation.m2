
--needsPackage"BlackBoxIdeals"


-- Type Map contains a matrix and a function
-- to evaluate this matrix on a point
-- This is for example used when projecting onto a
-- subspace (i.e elimination of variables)
Map = new Type of HashTable;

createMap = (mapMatrix, imageRing) -> (
    mapData := new MutableHashTable;
    
    mapData#"imageRing" = imageRing;
    mapData#"matrix" = mapMatrix;
    
    mapData#"valueAt" =  method();    
    mapData#"valueAt" (Matrix) := Matrix => (point)->
    (
        return sub(mapMatrix,point);
    );
   
    mapData#"valueAtJet" = method();
    mapData#"valueAtJet" (HashTable) := HashTable => (jet) -> (
        return new HashTable from {
           "failedJetLength" => jet#"failedJetLength",
           "jet" => (mapData#"valueAt")(jet#"jet"),
           "succeeded" => jet#"succeeded"
           };
     );
   
    return new Map from mapData
    )

--new Map from Matrix := (XXX, mapMatrix) -> (
--     sourceRing := ring mapMatrix;
--     K := coefficientRing sourceRing;
--     m := rank source mapMatrix;
--     imageRing := K[xxx_1..xxx_m]


TEST ///
assert (
     R = QQ[x,y,z];
     mm = createMap(matrix{{x^2,y^3}},QQ[a,b]);
     (mm#"valueAt")(matrix{{1,2,3}}) == matrix{{1,8}}
     )
///


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
interpolateBB = method();
interpolateBB(ZZ,BlackBoxParameterSpace,Matrix,Map) := Ideal => 
             (maxDegree,BB,point,mmap) -> 
     (
     R := mmap#"imageRing";
     mons := matrix {flatten apply(maxDegree+1,i->flatten entries basis(i,R))};
     -- find one jet with precision 10 more then number of monomials
     jetP := jetAt(BB,point,rank source mons+10,2);
     -- !!!this heuristic must be tested!!!
     -- Test: see if interpolated polynomials are in at least one
     -- irreducible component of the BlackBoxIdeal.
     jetPimage :=  (mmap#"valueAtJet")(jetP);
     --new HashTable from {
     --      "failedJetLength" => jetP#"failedJetLength",
     --      "jet" => (mmap#"valueAt")(jetP#"jet"),
     --      "succeeded" => jetP#"succeeded"
     --      };
     interpolate(mons,{jetPimage})
     )

interpolateBB(ZZ,BlackBoxParameterSpace,Matrix,Matrix) := Ideal => 
             (maxDegree,BB,point,mmap) -> 
(
     interpolateBB(maxDegree,BB,point,new Map from mmap)
)

interpolateBB(ZZ,BlackBoxParameterSpace,Matrix) := Ideal => 
             (maxDegree,BB,point) -> 
(
     interpolateBB(maxDegree,BB,point,createMap(vars BB.ring,BB.ring))
)

TEST ///
  K = ZZ/5
  R = K[x,y,z]
  I = ideal (x*z,y*z)
  bb = blackBoxIdeal I;     
  -- a point on the line  
  point = matrix{{0,0,1_K}}
  assert (ideal(x,y) == interpolateBB(1,bb,point))
///


doc ///
    Key
        interpolateBB
        (interpolateBB, ZZ, BlackBoxParameterSpace, Matrix)
        (interpolateBB, ZZ, BlackBoxParameterSpace, Matrix, Map)
    Headline
        find polynomials containing a list of jets
    Usage   
        I = interpolateBB(maxDegree,BlackBox,point)
        I = interpolateBB(maxDegree,BlackBox,point,map)
    Inputs  
        maxDegree:ZZ 
            the maximal degree of polynomials considered
        BlackBox:BlackBoxIdeal
        point: Matrix
            a point where the Blackbox vanishes    
        map: Map
            
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




InterpolatedIdeal = new Type of MutableHashTable;

new InterpolatedIdeal from MutableHashTable :=  (InterpolatedIdealAncestor,l)->
(  
     --type check?
     return l;
);

new InterpolatedIdeal from List :=  (InterpolatedIdealAncestor,l)->
(  
     return new InterpolatedIdeal from l;
);



createInterpolatedIdeal = method();
createInterpolatedIdeal( Ideal,ZZ,String ) := InterpolatedIdeal => 
                       (I,maxDegree,name)->
(
     return new InterpolatedIdeal from {
      "ideal" => I,
      "maxDegree" => maxDegree,
      "name" => name
     };
)

--createInterpolatedIdeal (ZZ,BlackBoxIdeal,Matrix) := InterpolatedIdeal => (maxDegree, BB, point)->
--(
--   createInterpolatedIdeal ( interpolateBB(maxDegree,BB,point), maxDegree, "" )
--)



InterpolatedImage = new Type of HashTable;

createInterpolatedImage = method();




-- observer observable längst fällig!

createInterpolatedImage(Experiment,Map) := HashTable => (experiment,imageRing, mapdata)->
(
    interpolation := new MutableHashTable;

    interpolation.experiment = () -> experiment;

    interpolation.blackBoxIdeal = () -> (interpolation.experiment()).blackBoxIdeal();
    

     -- duplicate code..(otherwise too many parameters...)
     interpolation.isOnComponent = method();
     interpolation.isOnComponent (HashTable,Matrix,ZZ) := Boolean =>  
                         (interpolationIdeal,point,precisionOfSmoothnessTest) -> (
         jetP := jetAt(interpolation.blackBoxIdeal() ,point,precisionOfSmoothnessTest,1);
         if jetP#"succeeded" then (
               0 == sub( interpolationIdeal, ((mapdata#"valueAtJet")(jetP))#"jet" )
          ) else error "point not smooth"
     );
    interpolatedIdeals := {};
    interpolation.createAllInterpolatedIdeals  = (maxDegree, precisionOfSmoothnessTest) -> 
    ( 
        localInterpolatedIdeals := {};
        for point in (interpolation.experiment()).points() do
        ( 
            if (interpolation.blackBoxIdeal()).isCertainlySingularAt(point) then continue;
             -- check if point is already on one of the known components
             bIsOnComponent := false;
             for interpolData in interpolatedIdeals do
             (
                 if  interpolation.isOnComponent (  interpolData.ideal, point, 0 ) then
                 (
                     if  interpolation.isOnComponent (  interpolData.ideal, point, precisionOfSmoothnessTest ) then
                     (
                          bIsOnComponent = true; 
                          break;
                     );
                 );
             );
             if bIsOnComponent then continue;

             localInterpolatedIdeals = localInterpolatedIdeals | { createInterpolatedIdeal (maxDegree, interpolation.blackBoxIdeal(), point)  };             
        );
        interpolatedIdeals = new MutableHashTable from 
           apply ( #localInterpolatedIdeals, idx-> (("ideal_" |toString idx ) => interpolatedIdeals#idx ) ) ;
       
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
 
        if (interpolation.blackBoxIdeal()).isCertainlySingularAt(point) then return "is not smooth";    

        for key in keys interpolatedIdeals do
        (
          if  interpolation.isOnComponent (  (interpolatedIdeals#key).ideal, point, 0,  ) then
          (
              if  interpolation.isOnComponent (  (interpolatedIdeals#key).ideal, point, prec,  ) then
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
         apply (keys interpolatedIdeals, key-> ( print (key=>new HashTable from interpolatedIdeals#key ) ) );
   );

    interpolation = newClass( InterpolatedImage, interpolation );
    return interpolation;
     
);

createInterpolatedImage(Experiment,Ring, Matrix) := HashTable => (experiment,imageRing, pmap)->
(
    mmap := new Map from pmap;
  
    
    return createInterpolatedImage(experiment,imageRing,mmap);
    
);

createInterpolatedImage(Experiment) := HashTable => (experiment)->
(
    mmap := new Map from vars (experiment.blackBoxIdeal()).ring  ;
    rng := (experiment.blackBoxIdeal()).ring ;
    
    return createInterpolatedImage(experiment, rng, mmap);
);

-- todo : new type for mapData. Also support dot access to mapData members.


end
----

restart

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

viewHelp interpolateBB

