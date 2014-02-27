-- finite field experiments

newPackage(
     "FiniteFieldExperiments",
     Version => "0.1", 
     Date => "15.02.2013",
     Authors => {{
                 Name => "Hans-Christian Graf v. Bothmer", 
           Email => "bothmer@math.uni-hannover.de", 
           HomePage => "http://www.crcg.de/wiki/Bothmer"},
           { Name => "Jakob Kroeker", 
           Email => "kroeker@uni-math.gwdg.de", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"}
      },
     Configuration => {},
     PackageExports => {"BlackBoxIdeals", "IntervalPkg" },
     Headline => "finite field experiments for explicit and implicitly given ideals and parameter spaces",
     DebuggingMode => false,
     AuxiliaryFiles=>true
)





export {
  "estimateDecomposition", 
  "estimateStratification",
  "estimateCodim", 
  "estimateNumberOfComponents",              
  "createInterpolatedIdeal",
  "interpolateBB",
  "interpolate",
  "isOnComponent",   
  "poissonEstimate",
  "Experiment",
  "RandomExperiment",
  "FFELogger",
  "ringCardinality"
}

FiniteFieldExperimentsProtect = ()->
(
  protect reset;
  protect setPointIterator;
  protect printInterpolatedIdeals;
  protect clearRecordList;
  protect testDebug;
  protect next;
  protect begin;
  protect count;
  protect point;
  protect collectedCount;
  protect pointKeys; 
  protect points;
  protect trials;
  protect createAllInterpolationIdeals;
  protect interpolatedIdealKeys;

  protect coefficientRingCardinality;
  protect pointLists;
  protect pointsByKey;
  protect countData;
  protect setIsInteresting;
  protect isInteresting;
  protect interpolatedIdeals;
  protect maxDegree;

  protect getExperimentData;
  protect setRecordedProperties;
  protect recordProperty;
  protect ignoreProperty;
  protect ignoreProperties;
  protect observedProperties;

  protect update;
  protect updateExperiment;

  protect saveData;
  protect loadData;

  protect propertyList;
  protect clear;
  protect jacobianAtKey;
  protect watchProperty;
  protect watchProperties;
  protect recordProperties;

  protect propertyName;
  protect propertyAt;

  protect tryProperty;
  protect bbi;

  protect recordedProperties;
  protect watchedProperties;
  protect useProperties;
  protect useJacobianAt;
  protect usedJacobianAt;

  protect minPointsPerComponent;
  protect setMinPointsPerComponent;
  protect stratificationIntervalView;
  protect countsByCount;  --protect countsByCount;

 
  protect  estimateStratification2;
  protect experimentData; 
  protect isRandom;
  protect compatible;
  protect membershipPrecision;
  protect setMembershipPrecision;
)

 

FiniteFieldExperimentsExport  = ()->
(
  exportMutable(reset);
  exportMutable (setPointIterator);
  exportMutable (printInterpolatedIdeals);
  exportMutable (membershipPrecision);
  exportMutable (setMembershipPrecision);

  exportMutable( clearRecordList);
  exportMutable( testDebug);
  exportMutable( next);
  exportMutable( begin);
  exportMutable( count);
  exportMutable( point);
  exportMutable( collectedCount);
  exportMutable( pointKeys);
  exportMutable( points);
  exportMutable( trials );
  exportMutable( interpolatedIdeals );
  exportMutable( createAllInterpolationIdeals );
  exportMutable(interpolatedIdealKeys);
 
  exportMutable(coefficientRingCardinality);
  exportMutable( pointLists );
  exportMutable( pointsByKey );
  exportMutable( maxDegree ) ;

  exportMutable( countData );

  exportMutable( setIsInteresting);
  exportMutable( isInteresting);

  exportMutable( getExperimentData);
  exportMutable( setRecordedProperties);
  exportMutable( recordProperty);
  exportMutable( ignoreProperty);
  exportMutable( ignoreProperties);
  exportMutable( observedProperties);


  exportMutable( update);
  exportMutable( updateExperiment);
  exportMutable( saveData);
  exportMutable( loadData);

 exportMutable( propertyList);
 exportMutable( clear);
 exportMutable( jacobianAtKey);
 exportMutable( watchProperties );
 exportMutable( watchProperty );
 exportMutable( recordProperties );
 exportMutable( propertyName);
 exportMutable( propertyAt);

 exportMutable( tryProperty );
 exportMutable( bbi );

 exportMutable( recordedProperties);
 exportMutable( watchedProperties);
 exportMutable( useProperties);
 exportMutable( useJacobianAt);
 exportMutable( usedJacobianAt);
 

  exportMutable( minPointsPerComponent);
  exportMutable( setMinPointsPerComponent);
  exportMutable( stratificationIntervalView);

 
  exportMutable( countsByCount);

 exportMutable (  estimateStratification2);
 exportMutable (  experimentData);
 exportMutable (  isRandom );
 exportMutable (  compatible );
 
 exportMutable (  createExperimentData );
) 

needsPackage "SimpleDoc";
needsPackage "Text";


exportMutable(savedExperimentData412398472398473923847 );

FiniteFieldExperimentsExport();

if FiniteFieldExperiments#Options#DebuggingMode then
  errorDepth=0;


FFELogger = Logger("FiniteFieldExperiments");

ffelog := FFELogger;



poissonEstimate = method( Options => {"confidence" => 1.96} );

poissonEstimate(ZZ) := HashTable => opts -> (numPoints) -> 
( 
     -- we use a poisson approximation since there are only 
     -- very few solution points compared to the number of trials
     --
     -- we then use the normal approximatiom of the poissondistribution
     -- to calculate the confidence interval
     --
     -- for a Poisson distribution the standard deviation is
     -- the square root of the expected value

     err := opts#"confidence"*sqrt(numPoints); 
     estimate := new Interval from   (  max(round(1,numPoints-err), 0.0001),    round(1,numPoints+err)  );
     return estimate;
     )

TEST ///
  loadPackage ("FiniteFieldExperiments",Reload=>true)
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert ( (poissonEstimate(16,"confidence" => 2)).min == 8  );
  assert ( (poissonEstimate(16,"confidence" => 2)).max == 24 );
///



 

TEST ///
  loadPackage ("FiniteFieldExperiments",Reload=>true)
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert (poissonEstimate(16,"confidence"=>2) == new Interval from (8,24))
///


estimateNumberOfComponents = method( Options => (options poissonEstimate) );

estimateNumberOfComponents( ZZ, RR, ZZ, ZZ ) := HashTable => opts->
    ( trials, estimatedCodim, numPoints, fieldCardinality ) -> 
(
    return poissonEstimate( numPoints,opts )*( (fieldCardinality^estimatedCodim)/trials*1.0  );
)


estimateNumberOfComponents( ZZ, ZZ, ZZ, ZZ ) := HashTable => opts->    (trials, estimatedCodim, numPoints, fieldCardinality) -> 
(
    estimateNumberOfComponents(trials,estimatedCodim*1.0,numPoints,fieldCardinality,opts)
)


TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  estimate =  estimateNumberOfComponents( 11^2*16, 2, 16, 11, "confidence"=>2 );
  estimate = estimate.round(1);
  assert ( estimate == new Interval from (0.5,1.5))
///

estimateCodim = method( Options => (options poissonEstimate) );

estimateCodim( ZZ, ZZ, ZZ ) := Interval => opts->
    ( trials, numPoints, fieldCardinality ) -> 
(
         est := (1/trials)*poissonEstimate(numPoints,opts);
         logEst := new Interval from (-log est.max,-log est.min);
         codimEst := (1/log fieldCardinality)*logEst;
     new Interval from (round(1,codimEst.min),round(1,codimEst.max))
);

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert (estimateCodim(11^2*10,10,11) == new Interval from (1.8,2.4))
///


 

InterpolatedIdeal = new Type of MutableHashTable;

new InterpolatedIdeal from MutableHashTable :=  (InterpolatedIdealAncestor,l)->
(  
     --type check?
     return l;
);


createInterpolatedIdeal = method();
createInterpolatedIdeal( Ideal,ZZ,String ) := InterpolatedIdeal => (I,maxDegree,name)->
(
     interpolatedIdeal := new MutableHashTable;
     interpolatedIdeal.ideal = I;
     interpolatedIdeal.maxDegree = maxDegree;
     interpolatedIdeal.name = name;
     
     return new InterpolatedIdeal from interpolatedIdeal;
)


load "./FiniteFieldExperiments/Interpolation.m2";

Experiment = new Type of HashTable;

ExperimentData = new Type of MutableHashTable;

PointData = new Type of MutableHashTable;

new PointData from HashTable := (ancestorType, pointData)->( 
    return pointData;
);


createPointData =(pBlackBox, point)->
(
    blackBox := pBlackBox;
    p := new MutableHashTable;
    p.point = point;
    
    localIsSingular := null;
    p.isSingular = ()->
    (
        if (localIsSingular=!=null) then return localIsSingular;
        localIsSingular = blackBox.isSingular(point);
        return localIsSingular;
    );
);

--new ExperimentData from HashTable := (E,coeffRing) -> (
new ExperimentData from Ring := (E,coeffRing) -> (
     print (toString E);
     e := new MutableHashTable;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable;
    -- format: key=>{ideal, maxDegree, name} --later: data type for interpolated ideal
     e.interpolatedIdeals = new MutableHashTable; 
     e.count = new Tally;
     e.trials = 0;
     e.propertyList = {};
     e.isRandom = null;
     return e;
);






estimateNumberOfComponents(Experiment,List) := HashTable => opts->    (experiment,key) -> 
(
      count := experiment.countData();
      posRankJacobianAt := experiment.position( "rankJacobianAt" );
     if posRankJacobianAt === null then error("To estimate number of components, \"rankJacobianAt\" must be watched");
    
     cardinality := experiment.coefficientRingCardinality();
     estimateNumberOfComponents(
      experiment.trials(),
      key#posRankJacobianAt,
      count#key,
      cardinality,opts)
)

createExperimentData = (coeffRing,points,count,trials,propertyList,isRandom) -> (
     e := new ExperimentData;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable from points;
     e.interpolatedIdeals = new MutableHashTable; 
     e.count = count;
     e.trials = trials;
     e.propertyList = propertyList;
     e.isRandom = isRandom;
     return e;
)

new ExperimentData from ExperimentData := (E,ed) -> (
     print (toString E);
     e := new MutableHashTable;
     e.coefficientRing = ed.coefficientRing;
     e.points = copy ed.points;
     e.interpolatedIdeals =  copy ed.interpolatedIdeals;
     e.count = copy ed.count;
     e.trials = ed.trials;
     e.propertyList = copy ed.propertyList;
     e.isRandom = ed.isRandom;
     return e;
     );

     ExperimentData == ExperimentData := Boolean=>(ed1,ed2)->
(

   if ( ed1.coefficientRing === ed2.coefficientRing and
            ed1.propertyList    == ed2.propertyList and
            keys ed1.points      == keys ed2.points and    
            keys ed1.count      == keys ed2.count and  
            ed1.isRandom        == ed2.isRandom  and
            ed1.trials        == ed2.trials  
) then 
  (
   
    
    for key in keys ed2.points do
    ( 
        if not ed1.points#key == ed2.points#key then 
          return false;
    );

     for key in keys ed2.points do
    ( 
        if not ed1.count#key == ed2.count#key then 
          return false;
    );
     return true;
  )
  else return false; 
);

--coeffRing,points,count,trials,propertyList,isRandom

toExternalString(ExperimentData) := String=> (ed)->
( 
   return "createExperimentData " | "(" 
   | toExternalString ed.coefficientRing | ", \n "
   | toString (new HashTable from ed.points) | ", \n"
   | toExternalString ed.count | ",\n"
   | toExternalString ed.trials | ",\n"
   | toExternalString ed.propertyList | ",\n"
   | toExternalString ed.isRandom 
   | ")" ;
)

ExperimentData + ExperimentData := ExperimentData=>(ed1,ed2)->
(

   if ( ed1.coefficientRing === ed2.coefficientRing and
        ed1.propertyList    == ed2.propertyList and
        keys ed1.count      == keys ed2.count and    
        ed1.isRandom        == ed2.isRandom  ) then 
  (
    edNew := new ExperimentData from ed1;
    
    for key in keys ed2.points do
    ( 
      if not edNew.points#?key then 
         edNew.points#key = copy ed2.points#key
      else
         edNew.points#key = edNew.points#key | copy ed2.points#key;
    );
     edNew.count = ed1.count + ed2.count;
     edNew.trials = ed1.trials + ed2.trials;
     return edNew;
  )
  else error ("+: experiment data not compatible");                   
);


createRandomPointIterator = method();

createRandomPointIterator (Function) := HashTable =>( weakRandomPointGenerator )->
(
        rpi := new MutableHashTable;
        randomPoint := null;
  
        currTrial := 0;

        rpi.next=()->
        (
            while(true) do 
            (
                 randomPoint = weakRandomPointGenerator();
                 currTrial = currTrial + 1;
                 if randomPoint =!= null then break;
            );
         
            return true;
        );
        rpi.position = ()->
        (
            return currTrial;
        );
        rpi.reset = () ->
        (
             randomPoint = null;
             currTrial = 0;
        );

        rpi.begin=()->
        (
            ri := createRandomPointIterator(weakRandomPointGenerator);
            ri.next();
            return ri;
        );

        rpi.point=()-> randomPoint;

       return new HashTable from rpi;
)

createRandomPointIterator (ZZ,Ring) := HashTable =>(numVariables, coeffRing )->
    (
        rpi := new MutableHashTable;
        randomPoint := null;
  
        currTrial := 0;

        rpi.next=()->
        (
            randomPoint = random( coeffRing^1, coeffRing^numVariables );         
            currTrial = currTrial + 1;
 
            return true;
        );

        rpi.reset = () ->
        (
             randomPoint = null;
             currTrial = 0;
        );

        rpi.position = ()->
        (
            return currTrial;
        );

        rpi.begin=()->
        (
            ri := createRandomPointIterator(numVariables, coeffRing);
            ri.next();
            return ri;
        );

        rpi.point=()-> randomPoint;

       return new HashTable from rpi;
    );

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  pointIterator = createRandomPointIterator(5,ZZ/7)
  pointIterator.next()
  pointIterator.point()
  apply(99, i-> pointIterator.next() )
///

createPointIterator = ( pPoints )->
    (
        pIterator := new MutableHashTable; 

        points := pPoints  ;
 
        point := null;

        pointCount := #points;
  
        currPosition := 0;

        pIterator.position = ()->
        (
            return currPosition;
        );

        pIterator.reset = () ->
        (
             point = null;
             currPosition = 0;
        );

        pIterator.next=()->
        (
            if (currPosition+1 >= pointCount) then return false;
            point = points#currPosition;
            currPosition = currPosition + 1;
            return (true) ;
        );

        pIterator.begin=()->
        (
            localPointIterator := createPointIterator( pPoints);
            localPointIterator.next();
            return localPointIterator;
        );

        pIterator.point=()-> point;

       return new HashTable from pIterator;
    );




estimateDecompositionOld := (experiment) -> (
       count := experiment.countData();
       posRankJacobianAt := experiment.position( "rankJacobianAt" );
       if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");

       cardinality := experiment.coefficientRingCardinality();

       print "(estimated codim, estimated number of components [confidence interval] <= {watched Properties})";
       print "";
       apply(sort apply(keys count,key->
         (net(
              key#posRankJacobianAt,
              estimateNumberOfComponents(
               experiment.trials(),
               key#posRankJacobianAt,
               count#key,
               cardinality ) )
         ) |" <= " |net key ),
       print);
   );
 

estimateDecomposition =  (experiment) -> (
        posRankJacobianAt := experiment.position( "rankJacobianAt" );
       if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");
       print "(estimated codim, estimated number of components [confidence interval] <= {watched Properties})";
       print "";
       apply(sort apply(keys experiment.countData(),key->
              net( key#posRankJacobianAt, estimateNumberOfComponents(experiment,key) )
         |" <= " | net key),
       print);
)

--needs to be documented
 

 
estimateStratification =  (experiment) -> (
     trials := experiment.trials();
     orderK := experiment.coefficientRingCardinality(); -- this must be read from the experimentdata
     -- (jk): need more advice. Did we want to use a different ring for search that the ideal coefficient ring? If so, 
     print "--";
     print "-- estimated codim <= {wachtched properties}";
     print "--";
     count := experiment.countData();
     -- sort keys by number of occurence
     sortKeysCount := apply(reverse sort apply(keys count,k->(count#k,k)),i->i#1);
     apply(sortKeysCount,k->(
           --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
           print (net(round(1,(log(trials)-log(count#k))/log(orderK)))|" <= "|net(k)))
           )
      ;
     print "--";
)


   stratificationIntervalView := (stratificationData )->
   (
     
      prerearrangedData := new HashTable from apply( keys stratificationData, key -> (stratificationData#key=>(stratificationData#key,key )));
      toSort := apply( keys stratificationData, key -> (stratificationData#key));
      sorted := sort (toSort); 
      rearrangedData := apply (sorted, key-> (prerearrangedData#key) );
      return new HashTable from rearrangedData;
   );



   --# fuer den Fall dass dich die Schlüssel nicht sortieren lassen.
   countsByCount := (experimentData)->
   (
      counts := experimentData.count;
      prerearrangedData := new MutableHashTable;
     
      for key in keys counts do
      (
          if not  prerearrangedData#?(counts#key) then 
             prerearrangedData#(counts#key)= { key }
          else
          (
              prerearrangedData#(counts#key)=  prerearrangedData#(counts#key) | { key };
          );

      );
      --bug fixed (unique missing) todo: test for this bug!
      toSort := unique apply( keys counts, key -> (counts#key));
      sorted := sort (toSort); 
      rearrangedData := {};
      for count in sorted do
      (
          for entry in prerearrangedData#count do
          rearrangedData = rearrangedData | {(count,entry)};
      );
     
      return new List from rearrangedData;
   );


 
ringCardinality = (rng)->
(
   cardinality := null;
   if char rng ==0 then   cardinality  = infinity;
   try { cardinality = (rng).order } then () else();
   return cardinality;
);
 

new Experiment from Thing := (E, thing) -> 
(
   error("not impelemted");
)


RandomExperiment = new Type of Experiment;


new RandomExperiment from Thing := (E, thing) -> 
(
   error("not impelemted");
)


 


RandomExperiment + RandomExperiment := RandomExperiment => (re1,re2)->
(
    bb1 := re1.blackBoxIdeal();
    bb2 := re2.blackBoxIdeal();

    if  re1.compatible(re2)    then
   (
     re := copy re1;
     re.merge(re2);
       return re;
   )
   else error ("experiments not compatible");
);


new Experiment from BlackBoxParameterSpace := (E, pBlackBox) -> 
(


   blackBoxIdeal := pBlackBox;    --black box ideal or black box parameter space so far
   experimentData := new ExperimentData from blackBoxIdeal.coefficientRing;

   coefficientRingCardinality := ringCardinality( blackBoxIdeal.coefficientRing );


   experiment := new MutableHashTable;

   -- todo: issue about naming here. blackBox (what?)
   --experiment.bbi = blackBoxIdeal; -- is this safe?
   --experiment.blackBoxIdeal = blackBoxIdeal;


   -- todo: maype return a (deep) copy 
   experiment.getExperimentData = ()->
   (
      print(" warning : you got a reference to experiment data, do not modify! ");
      return experimentData;
   );

   -- some of the  following values initialized at end ( e.g. propertiesAt initialization depends presence of some functions defined later)
   minPointsPerComponent := 10;
   jacobianAtKey := null;
   jacobianAt := null;
   propertiesAt := (point)->{};
   isInteresting := (point)->true;

   pointIterator := createRandomPointIterator ( blackBoxIdeal.numVariables, experimentData.coefficientRing );

   experiment.experimentData=()->
   (
      return new ExperimentData  from experimentData;
   );

   experiment.blackBoxIdeal=()->
   ( 
      return blackBoxIdeal;
   );


   experiment.estimateDecomposition = () -> (estimateDecomposition(experiment));


   experiment.estimateStratification = () -> (estimateStratification(experiment));

--   experiment.estimateStratification2 = () -> ( estimateStratification2(experiment) );

   experiment.compatible = method();
     experiment.compatible (RandomExperiment) := Boolean =>(re2)->
    (
        bb1 := experiment.blackBoxIdeal();
        bb2 := re2.blackBoxIdeal();

     return (  experiment.minPointsPerComponent() ==  re2.minPointsPerComponent() and
               experiment.coefficientRing()       === re2.coefficientRing() and
               experiment.watchedProperties()     ==  re2.watchedProperties() and
               experiment.jacobianAtKey()         ==  re2.jacobianAtKey() and
               --experiment.isInteresting           ===  re2.isInteresting -- um sicher zu gehen, kann man  alle punkte in re2 durch isInteresting von e jagen und umgekehrt.
               experiment.pointKeys()             ==  re2.pointKeys() and

              bb1.numVariables            ==   bb2.numVariables and
              bb1.numGenerators           ==   bb2.numGenerators and
              bb1.coefficientRing         ===  bb2.coefficientRing 
             );
    );


   experiment.merge = method();

   -- todo: how to prevent from self-merging?
   experiment.merge (RandomExperiment)  := RandomExperiment => (re)->
   (
      if (experimentData==re.experimentData() ) then 
          error ("attempt to merge with itself");

      if   experiment.compatible(re) then 
        experimentData = experimentData + re.experimentData()
      else
       error ("experiments not compatible!");
   );
   
   experiment.stratificationIntervalView = ()->
   (
        stratificationData := experiment.estimateStratification();
        return stratificationIntervalView(stratificationData);
   );

   experiment.coefficientRing=()->
   (
     
       return experimentData.coefficientRing;
   );
  

   experiment.countsByCount = ()->
   (  
       return countsByCount( experimentData ); 
   );



   experiment.coefficientRingCardinality=()->
   (
       -- return    ringCardinality( blackBoxIdeal.coefficientRing );
       return    coefficientRingCardinality;
   );

 
   experiment.testDebug=()->
   (
      a:=5;
      1/0;
      return a;
   );

    
   experiment.isInteresting=(point)->
   ( 
      return isInteresting(point);
   );

   experiment.points=()->
   (
      return flatten  apply (keys experimentData.points, countKey-> experimentData.points#countKey);
   );



   runExperimentOnce := method();




   experiment.setMinPointsPerComponent = (numPointsPerComponent)->
   ( 
     minPointsPerComponent = numPointsPerComponent;
   );


   experiment.minPointsPerComponent = ()->
   ( 
     return minPointsPerComponent ;
   );


  

   -- todo: updateExperiment?

   -- maybe think about writing 'connectProperty' ... (propName, bb.propName); default is 1:1.
   -- then, what should happen if a user requests to watch property xy ?
   
   experiment.useJacobianAt = (jacobianAtName)->
   ( 
       if experiment.trials()=!=0 then error ("cannot change jacobianAt  - experiment was already run! You could clear() the statistics and retry. ");

       if jacobianAtName===null then
       (
         jacobianAtKey = jacobianAtName;
         jacobianAt = jacobianAtName ;
         return;
       );
       if blackBoxIdeal.hasPointProperty(jacobianAtName) then 
       (
          jacobianAtKey = jacobianAtName;
          jacobianAt = blackBoxIdeal.pointProperty(jacobianAtName) ;
       ) 
       else  error ("blackBoxIdeal seems not to have property" | propertyName );
   );


   experiment.position =(property)->
   (
     return position( experiment.watchedProperties(), watchedProperty->watchedProperty == property );
   );




   experiment.clear = ()->
   ( 
        pointIterator.reset();
        --  later  call experimentDataClear() instead.
        experimentData.trials = 0;
        experimentData.points = new MutableHashTable;
        experimentData.count = new Tally;
   );

 
   experiment.jacobianAtKey = ()->
   (
       return jacobianAtKey;
   );

  experiment.usedJacobianAt = ()->
   (
       return jacobianAtKey;
   );
 

    runExperimentOnce(ExperimentData, Matrix,ZZ ) := Tally => (experimentData, point, wantedPoints) -> 
    (  
        K := experimentData.coefficientRing;
        --prime = char K
        numVariables := blackBoxIdeal.numVariables;

        -- if ideal vanishes on the random point do something
        if experiment.isInteresting(point) then   
        (
            countKey := propertiesAt(point);
 
            -- count number of found points for each rank and property
            experimentData.count = experimentData.count + tally {countKey};

          
            if  jacobianAt =!= null then 
            (
                FFELogger.debug( "update wanted points" );
                rankJacobian := rank  jacobianAt(point); 
                upperEstimate := (estimateNumberOfComponents(experiment,countKey)).max;
                ffelog.debug ("upper estimate number of components:  " | toString upperEstimate );
                --upperEstimate := ((estimateDecompositionOld( experimentData ))#countKey).max;
                --upperEstimate := 1; -- test

                wantedPoints = max(1,upperEstimate)*wantedPoints;
            );            

          -- remember some points
            if experimentData.points#?(countKey) then 
            (
                -- collect a fixed number of points per estimated component
                -- use upper limit of estimation for this
                if #(experimentData.points#countKey) < wantedPoints then 
                (
                    FFELogger.debug( "attaching point" );
                    experimentData.points#countKey = experimentData.points#countKey | {point};
                );
            )
            else (
                FFELogger.debug( "attaching first point for some key");
                experimentData.points#countKey = {point};
            );
        );
        -- this trial is counted in runExperiments and not here to allow update() without changing trial num.
    );

   
    experiment.setPointIterator  = (pRpi)->
    (
       if experiment.trials()=!=0 then error ("cannot change random iterator - experiment was already run! You could call clear() and retry.");
        pointIterator = pRpi;
    );

     experiment.trials = ()-> pointIterator.position();

     runExperiment := method();
     runExperiment(ExperimentData, Thing, ZZ) := Tally => (experimentData, pPointIterator, newTrials) -> 
     ( 
       for i in 1..newTrials do
       (
           assert( pPointIterator.next() );
           runExperimentOnce( experimentData, pPointIterator.point(), minPointsPerComponent );
           experimentData.trials =  experiment.trials();
       );
     );

 

     update := method();
     update(ExperimentData) := Tally => opts -> (experimentData) -> 
     ( 
        pointIterator := createPointIterator (  experiment.points() );
        experimentData.points = new MutableHashTable;
        experimentData.count = new Tally;

        while ( pointIterator.next() ) do
        (
            runExperimentOnce( experimentData, pointIterator.point(), minPointsPerComponent );
        );
     );




   setRecordedPropertiesInternal := (propListToObserve)->
   ( 
      for propertyName in propListToObserve do
      (
          if not blackBoxIdeal.hasPointProperty(propertyName) then 
              error ("blackBoxIdeal seems not to have property" | propertyName );
      );
      experimentData.propertyList=propListToObserve;
      propertiesAt = (point)->
      ( 
        apply( experimentData.propertyList, propertyName->( (blackBoxIdeal.pointProperty(propertyName))(point) ) )  
      );   
   );

  experiment.clearRecordList = (   )->
  (
      setRecordedPropertiesInternal({});
  );


   experiment.setRecordedProperties = ( propertyStringList )->
   (  
      if experiment.trials()=!=0 then error ("cannot change watched properties - experiment was already run! You could clear() the statistics and retry.");

      setRecordedPropertiesInternal(propertyStringList);
      update(experimentData);
   );

   UpdateRecordedPropertiesError := "cannot change watched properties - experiment was already run! You could clear() the statistics and retry.";
  

   experiment.watchProperties = experiment.setRecordedProperties;

   experiment.watchProperty=(propertyName)->
   (
       if experiment.trials()=!=0 then error (UpdateRecordedPropertiesError);

          if not blackBoxIdeal.hasPointProperty(propertyName) then 
              error ("blackBoxIdeal seems not to have property" | propertyName );

      experimentData.propertyList = unique (experimentData.propertyList | { propertyName }) ;
      setRecordedPropertiesInternal( experimentData.propertyList );   
      update(experimentData);
   );

   experiment.recordProperty=experiment.watchProperty;

   experiment.ignoreProperty=(propertyName)->
   (
       if experiment.trials()=!=0 then error (UpdateRecordedPropertiesError);

      experimentData.propertyList = delete(propertyName, experimentData.propertyList ) ;
      setRecordedPropertiesInternal( experimentData.propertyList );   
        
      update(experimentData);
   );


   experiment.ignoreProperties = (ignorePropertyStringList)->
   (
       if experiment.trials()=!=0 then error (UpdateRecordedPropertiesError);

      apply( ignorePropertyStringList, propToIgnore-> ( experimentData.propertyList = delete(propToIgnore, experimentData.propertyList ); ));
      setRecordedPropertiesInternal( experimentData.propertyList );    
      update(experimentData);
   );


   -- maybe observed Properties is not a good name - is 'recordedProperties' better ?
   experiment.recordedProperties =  ()->
   (
      return experimentData.propertyList;
   );

   experiment.watchedProperties = experiment.recordedProperties;

   experiment.setIsInteresting = (pIsInteresting)->
   (  
      if ( pIsInteresting=!=isInteresting ) then
      (
         isInteresting = pIsInteresting;
         update(experimentData);
      );
   );

   ---- syntax for the moment too hard (method without parameters)
   --experiment.runExperimentOnce = method(Options => (options runExperimentOnce));
   -- experiment.runExperimentOnce() := Thing => opts->()->
   --(
   --   return runExperimentOnce(experimentData);
   --);

  
  

   experiment.run = method();
   experiment.run(ZZ) := Thing=> (newTrials)->
   (
       runExperiment( experimentData, pointIterator, newTrials );
       return experiment.countData();
   );
  
 
   experiment.pointLists = ()->
   (
      return new HashTable from experimentData.points;
   );

   experiment.pointsByKey = (key)->
   (
      if not (experimentData.points)#?key then 
         error "invalid key";
      return  (experimentData.points)#key;
   );
 
   -- returns a HashTable with recorded BlackBoxIdeal properties as keys and corresponding occured count as values.
   experiment.countData = ()->
   (
      return new Tally from experimentData.count;
   );


   experiment.pointKeys = ()->
   (
      return keys experimentData.points;
   );

   experiment.collectedCount = ()->
   (
      return new HashTable from apply( experiment.pointKeys(), key->( key=> #(experimentData.points)#key ) );
   );
   

   --init part:

   -- todo: test if changing blackBoxIdeal.jacobianAt is transparent (means after blackBoxIdeal.updatePointProperty("jacobianAt") the new one is called)
 
    experiment.watchProperties( {} );

       if  blackBoxIdeal.hasPointProperty("isZeroAt") then
       (
          --print("isInteresting"| toString  (blackBoxIdeal.isZeroAt ));
          isInteresting = blackBoxIdeal.pointProperty("isZeroAt");

          --     experiment.watchProperties( {"isZeroAt"} );
       );

       if blackBoxIdeal.hasPointProperty("jacobianAt") then 
       (
          jacobianAtKey = "jacobianAt";
          jacobianAt =  blackBoxIdeal.pointProperty("jacobianAt");
          if  not blackBoxIdeal.hasPointProperty("isZeroAt") then 
              error(" when for a balck box the property 'jacobianAt' is available, it is also expected that the black box has the point property 'isZeroAt' ");  
       );

 
       if blackBoxIdeal.hasPointProperty("rankJacobianAt") then
       ( 
            --experiment.watchProperties( {"isZeroAt","rankJacobianAt"} );
            experiment.watchProperties( {"rankJacobianAt"} );
       );

   --end init:


  -- to fix the issue that the internal experiment reference
  -- is not of type Experiment, there are two solutions:
  --  either to introduce to different method signatures (ones accept a HashTable and others an Experiment)
  -- or to make the internal experiment variable as an Experiment , too  - done with 'newClass'

    experiment.tryProperty = (tProp) -> (
      print("-- ( " | toString experiment.watchedProperties() |" | " | tProp | " ) => count " );
       pointLists = experiment.pointLists();
       tally flatten apply(keys pointLists, key->apply(pointLists#key, point->( key, (blackBoxIdeal.pointProperty(tProp))(point))) )
     );
   
    experiment.saveData  = (filename) ->
    (
       f := openOut(filename);        
       --maybe not that smart...
       externalString := toExternalString experiment.experimentData(); 
       --also not that smart...
       f << "savedExperimentData412398472398473923847 = " << externalString;
       f << flush;
       f << close;
    );

    experiment.loadData = (filename) ->
    (
      load filename;
      experimentData = savedExperimentData412398472398473923847;
      apply( keys experimentData.points, key-> 
                                    ( (experimentData.points)#key = apply( (experimentData.points)#key , point->
                                                                             sub(point, experimentData.coefficientRing )
                                                                         )
                                    )
               );
    );
    
    experiment.createAllInterpolationIdeals  = (maxDegree, prec) -> 
    ( 
        interpolatedIdeals := {};
        for point in experiment.points() do
        ( 
            if blackBoxIdeal.isSingular(point) then continue;
             -- check if point is already on one of the known components
             bIsOnComponent := false;
             for interpolData in interpolatedIdeals do
             (
                 if  isOnComponent ( blackBoxIdeal, interpolData.ideal, point, 0 ) then
                 (
                     if  isOnComponent ( blackBoxIdeal, interpolData.ideal, point, prec ) then
                     (
                          bIsOnComponent = true; 
                          break;
                     );
                 );
             );
             if bIsOnComponent then continue;

             interpolatedIdeals = interpolatedIdeals | { createInterpolatedIdeal (maxDegree, blackBoxIdeal, point)  };             
        );
        experimentData.interpolatedIdeals = new MutableHashTable from 
           apply ( #interpolatedIdeals, idx-> (("ideal_" |toString idx ) => interpolatedIdeals#idx ) ) ;
       
    );
    
    experiment.interpolatedIdealKeys = method();
    experiment.interpolatedIdealKeys (Matrix,ZZ) := Thing => (point,prec)->
    (
       numbers := {};
 
        if blackBoxIdeal.isSingular(point) then return "is not smooth";    

        for key in keys experimentData.interpolatedIdeals do
        (
          if  isOnComponent ( blackBoxIdeal, (experimentData.interpolatedIdeals#key).ideal, point, 0 ) then
          (
              if  isOnComponent ( blackBoxIdeal, (experimentData.interpolatedIdeals#key).ideal, point, prec ) then
              (
                 numbers = numbers | {key};
              );
           );
        );
        return numbers;
 
    );

    localMembershipPrecision := 10; -- does this belong to experimentData?


    experiment.setMembershipPrecision = (prec)->
    (
         localMembershipPrecision = prec;
    );

    experiment.membershipPrecision = (prec)->
    (
         localMembershipPrecision 
    );


    experiment.interpolatedIdealKeys (Matrix) := Thing => (point)->
   (
        return  experiment.interpolatedIdealKeys(point, localMembershipPrecision );
   );

   iik := (point)-> (return  experiment.interpolatedIdealKeys(point) ;);

   blackBoxIdeal.rpp("interpolatedIdealKeys", iik);



    experiment.printInterpolatedIdeals = ()->
   (
         apply (keys experimentData.interpolatedIdeals, key-> ( print (key=>new HashTable from experimentData.interpolatedIdeals#key ) ) );
   );


   experiment = newClass( Experiment, experiment );
   ffelog.debug ("type of internal experiment variable is : " | toString class experiment );
   return experiment; 
   -- return new HashTable from experiment; 
);

 


-- tryProperty = (experiment, property) ->(
--     pointLists := experiment.pointLists();
--     apply(
--           apply((keys pointLists),key -> (key,tally apply(pointLists#key,property))),
--           i->print (net i#0|" => "|net i#1)
--           )
--     )

doc ///
   Key
        Experiment
   Headline
        type to collect data from a finite field experiment
   Usage   
        e = new Experiment from bb
        e.run(trials)
        e.countData()
        e.pointLists()
        e.estimateStratification()
        e.estimateDecomposition()
   Inputs  
        bb:HashTable 
            a blackBoxIdeal
   Description
        Text
            Creates an experiment from a black box, see @TO BlackBoxIdeals@ . The experiment object will keep
            track of all information created when the ideal is evaluated at random points.
            Here is a typical extremely simple minded application:
                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5
           R = K[x,y,z]
           I = ideal (x*z,y*z)
           bb = blackBoxIdeal I;       
        Text
           \break The ideal describes a line and a plane intersecting at the origin.
       
           The experiment will evaluate the ideal at random points. If
           the point is in the vanishing set of the ideal (i.e. either on the point
           or on the line) it will calculate the rank of the jacobi matrix at that point.
           (2 on the line, 1 on the plane, 0 in the origin). 
        Example
           bb.isZeroAt(matrix{{0_K,0,1}})
           bb.rankJacobianAt(matrix{{0_K,0,1}})
           bb.rankJacobianAt(matrix{{1_K,2,0}})
           bb.rankJacobianAt(matrix{{0_K,0,0}}) 
        Text
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
        Text
           \break If a black box has a property "rankJacobianAt" it is
           automatically watched:
        Example
           e.watchedProperties() 
        Text
           \break Now we run the experiment by evaluating at 1250 random points:  
        Example
           time e.run(1250)      
        Text
           There are 125 points in (F_5)^3 of which 25 are on the
           plane, 5 are on the line and 1 (the origin) is on the line and the plane.
           We therefore expect about 10 points with rankJacobian = 0, 240 with rankJacobian = 1
           and 40 with rankJacobian = 2.
           
           Often it is useful to have explicit points on all components of a variety. Therefore
           the experiment has stored some of the points:
        Example
           e.pointLists()
        Text
           Since one always finds many points found on components of low
           codimension it is not useful to remember all of them. The experiment
           remembers by default only about 10 points per component:
        Example
           e.collectedCount() 
        Text
           The amount of stored points can be adjusted:
        Example
           e.setMinPointsPerComponent(20)
           -- collect about 20 points per component now:
           time e.run(1250);
           e.collectedCount() 
        Text
           Here we have not collected exactly 20 points per component since
           the experiment uses the upper end of the confidence interval for the
           number of components as guide for the number of points to keep. 
           
           Lets now estimate the number and codimension of reduced components
           of the vanishing set:
        Example
           e.estimateDecomposition()       
        Text
           Indeed we see that there is probably no component of codimension 0 and about
           1 component of codimension 1 and 2. The count for the component of 
           codimension 2 is slightly lower since we see only 4 of the 5 points
           on the line - one is the origin that has a tangent space smaller 
           codimension. So we expect the heuristic count to be 4/5 = 0.8. The
           same problem occurs for the plane but there it is not so relevant. We
           expect 24/25 = 0.96 in this case. For higher characteristics this
           systematic error gets smaller.
           
           If one finds the confidence intervals to large one can continue
           to run the experiment:
        Example
           time e.run(1250)
           e.estimateDecomposition()     
        Text
           A doubeling of the number of experiments is expected to devide the
           width of a confidence interval by the square root of 2.         
   Caveat
        does not check if the ideal ring is a quotient ring (not supported)
       
///


TEST ///
    debug FiniteFieldExperiments
    FiniteFieldExperimentsProtect()
    coeffRing := ZZ/3;
    bbRankM = blackBoxParameterSpace( 5 ,coeffRing )
    rankMat := (point)->5
    bbRankM.registerPointProperty("rankJacobianAt",rankMat)


    point := matrix {{1,2,3,4,5}};
    point = sub( point, coeffRing);

 
    bbRankM = bbRankM.getUpdatedBlackBox(bbRankM)
 

    e = new Experiment from bbRankM
    assert (e.coefficientRing()===coeffRing);

    e.setMinPointsPerComponent(20);
    assert( e.minPointsPerComponent()==20);
    FFELogger.setLogLevel(4);
    e.watchProperties {"rankJacobianAt"};
    e.watchedProperties()
    assert( 1== # select( e.watchedProperties() , 
                       (prop)->(prop=="rankJacobianAt") ) 
     )
    e.useJacobianAt("rankJacobianAt");
    e.useJacobianAt(null);
   
    e.countsByCount()
    points := e.points();
    #points
    apply(points,point->rankMat(point))
    FFELogger.setLogLevel(2);
    time e.run(1000)
    assert (e.trials()==1000);

    e.estimateStratification()
    e.estimateDecomposition()
  
    e.collectedCount()
    e.watchedProperties()
    e.jacobianAtKey()

    bbRankM.knownPointProperties()

    -- e.stratificationIntervalView() -- test fails here; d

///


end
---

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

viewHelp FiniteFieldExperiments

restart

loadPackage"BlackBoxIdeals"
load "FiniteFieldExperiments.m2"

R = (ZZ/7)[x_0..x_3]
M = matrix{
     {x_0,x_1,0},
     {x_1,x_2,x_3}
     }
I = minors(2,M)
B = blackBoxIdeal I
 
e = new Experiment from B

e.run(1, "numPointsPerComponentToCollect"=>20 ) 
e.run(1) 

e.run( 3000,  "numPointsPerComponentToCollect"=>20 )

pointData = e.getPointLists()

apply(keys pointData,i->#(pointData#i))

assert (#(pointData#{2}) >= 40)

