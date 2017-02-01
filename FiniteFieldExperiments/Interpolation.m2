
--needsPackage"BlackBoxIdeals"



TEST ///
assert (
     R = QQ[x,y,z];
     mm = createMapHelper(matrix{{x^2,y^3}},QQ[a,b]);
     (mm#"valueAt")(matrix{{1,2,3}}) == matrix{{1,8}}
     )
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


TEST ///
  K = ZZ/5
  R = K[x,y,z]
  I = ideal (x*z,y*z)
  bb = blackBoxIdeal I;     
  -- a point on the line  
  point = matrix{{0,0,1_K}}
  assert (ideal(x,y) == interpolateBB(1,bb,point))
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
  assert (SingularPointException === class catch interpolateBB(1,bb,intersectionPoint))
  -- a point over the wrong field
  assert (try interpolateBB(1,bb,pointQQ) then false else true) 
///





InterpolatedImage = new Type of HashTable;

createInterpolatedImage = method();




-- observer observable längst fällig!

createInterpolatedImage(Experiment,Ring, MapHelper) := HashTable => (experiment,imageRing, mapdata)->
(
    interpolation := new MutableHashTable;

    interpolation.experiment = () -> experiment;

    bb := (interpolation.experiment()).blackBoxIdeal();

    interpolation.blackBoxIdeal = () -> (interpolation.experiment()).blackBoxIdeal();
   
    jets := new MutableHashTable;

     -- duplicate code..(otherwise too many parameters...)
     interpolation.isOnComponent = method();
     interpolation.isOnComponent (HashTable,Matrix,ZZ) := Boolean =>  
                         (interpolationIdeal,point,onComponentPrecision) -> (

         if not (jets#?point) then
         (        
             --jetP := jetAt(interpolation.blackBoxIdeal() ,point,onComponentPrecision,1);
             jets#point = jetAt( bb ,point,onComponentPrecision, 1);
         );
         jetP := jets#point;  

         if jetP#"jetLength" < onComponentPrecision then 
         (
             jetP = jets#point;  
             jets#point = jetAt( bb ,point,onComponentPrecision, 1);
         );

         if jetP#"succeeded" then (
               0 == sub( interpolationIdeal, ((mapdata#"valueAtJet")(jetP))#"jet" )
          ) else error "point not smooth"
     );
    interpolatedIdeals := {}; 
    -- better: interpolatedIdeals:= new MutableHashTable; ? but then we do not know if it was initialized??
    interpolation.createAllInterpolatedIdeals  = (maxDegree, onComponentPrecision) -> 
    ( 
        idealCount := 0;
        localInterpolatedIdeals := {};
        bb := interpolation.blackBoxIdeal();
        -- T := timing 
        for point in (interpolation.experiment()).points() do
        ( 
             --time if bb.isCertainlySingularAt(point) then 
             --(
             --   continue;
             --);
             pointIsSingular := false;
             --
             -- check if point is already on one of the known components
             bIsOnComponent := false;
             for interpolData in localInterpolatedIdeals do
             (
                 try (  
                        if  interpolation.isOnComponent (  interpolData#"ideal", point, 0 ) then
                         (
                             if  interpolation.isOnComponent (  interpolData#"ideal", point, onComponentPrecision ) then
                             (
                                  bIsOnComponent = true; 
                                  break;
                             );
                         ); 
                     )  
                    else (
                    print ("pointIsSingular");
                    pointIsSingular = true;
                    break;
                );
             --       
             );
             if (bIsOnComponent or pointIsSingular) then 
             (
                 continue;
             );
             --debug
             --print(localInterpolatedIdeals);
             --print(point);
             
             interpolatedIdealOrException := catch createInterpolatedIdeal (maxDegree, interpolation.blackBoxIdeal(), point, ("ideal_" |toString idealCount )) ;          
             if ( class interpolatedIdealOrException ===  SingularPointException) then
             (
                print( interpolatedIdealOrException); --debug
                FFELogger.debug("createInterpolatedIdeal: point"| toString point| " was singular");
             )
             else
             (
               interpolatedIdeal := interpolatedIdealOrException;
               localInterpolatedIdeals = localInterpolatedIdeals | { interpolatedIdeal };             
               idealCount = idealCount +1 ;
             );
        );
        -- print "timing for loop", T#0;        
        interpolatedIdeals = new MutableHashTable from 
           apply ( #localInterpolatedIdeals, idx-> (("ideal_" |toString idx ) => localInterpolatedIdeals#idx ) ) ;           
       --print ("created interpolatedIdeals");
       FFELogger.debug("created interpolatedIdeals");
       return interpolation.interpolatedIdeals();
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
          if  interpolation.isOnComponent (  (interpolatedIdeals#key)#"ideal", point, 0,  ) then
          (
              if  interpolation.isOnComponent (  (interpolatedIdeals#key)#"ideal", point, prec,  ) then
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


   interpolation.interpolatedIdeals = ()->
   ( 
       --print (interpolatedIdeals); --debug
       --print (keys interpolatedIdeals); --debug
       ll := apply (keys interpolatedIdeals, key-> ( key=>new InterpolatedIdeal from interpolatedIdeals#key ));
       return new HashTable from ll;
   );

   interpolation.bareIdeals = ()->
   (
         apply (keys interpolatedIdeals, key-> ( (interpolatedIdeals#key)#"ideal" ) )
   );


    interpolation = newClass( InterpolatedImage, interpolation );
    return interpolation;
     
);

createInterpolatedImage(Experiment,Ring, Matrix) := HashTable => (experiment,imageRing, pmap)->
(
    error "this interface is currently broken";
    mMapHelper := new MapHelper from pmap;
     
    return createInterpolatedImage(experiment, imageRing, mMapHelper);
    
);

createInterpolatedImage(Experiment) := HashTable => (experiment)->
(
    rng := (experiment.blackBoxIdeal()).ring ;
 
    mMapHelper  := new MapHelper from vars rng;

    return createInterpolatedImage(experiment, rng, mMapHelper);
);

interpolatedImage = method();

interpolatedImage(Experiment) := InterpolatedImage => (e)->
(
   return createInterpolatedImage(e);
)

new InterpolatedImage from Experiment :=  (InterpolatedImage,e)->
(
    return createInterpolatedImage(e);
)

-- todo : new type for mapData. Also support dot access to mapData members.

doc ///
    Key
        createAllInterpolatedIdeals
    Headline
        find components of a black box ideal
    Usage   
        i.createAllInterpolatedIdeals(maxDeg,prec)
    Inputs  
        i:InterpolatedImage 
        maxDeg:ZZ
            interpolate up to this degree
        prec: ZZ
            pecision of ideal membership test                
    Description
        Text
           For each point found in an experiment find
           all polynomials of degree at most maxDeg
           that contain the component on which the point lies.

           Lets consider a black box that describes
           a conic and a plane intersecting at the origin:
        Example      
           K = ZZ/5
           R = K[x,y,z]
           I = ideal (x*z,(y^2-z)*z)
           bb = blackBoxIdeal I;       
        Text
           Now we create an experiment to find points on 
           all components
        Example
           e = new Experiment from bb;
           e.run(500)
           e.estimateDecomposition()
        Text
           The experiment reveals approximately 1 component of
           codimension 1 and 2 respectively.           
           \break
           Now we want to find linear equations containing the respective
           components on which the points lie. For this we need
           to create an interpolatedImage object which will among
           other things contain the experiment and (later) the equations
           of each component
        Example
           i = createInterpolatedImage(e);
        Text
           Now do the interpolation looking only for linear polynomials
        Example
           i.createAllInterpolatedIdeals(1,1);
        Text 
           Only linear polynomials were found. This is remembered
           behind the scenes.
        Text
           If we also want to find the quadratic polynomial,
           we have to interpolate up to degree 2
        Example
           i.createAllInterpolatedIdeals(2,1);
           i.bareIdeals()
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


end
----

restart

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

viewHelp interpolateBB

