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
     Headline => "black boxes for explicit and implicitly given ideals",
     DebuggingMode => true
)


needsPackage "IntervalPkg";
needsPackage "BlackBoxIdeals";


export {
  "estimateDecomposition", 
  "estimateStratification",
  "estimateCodim", 
  "estimateNumberOfComponents",    	    	 
  "poissonEstimate",
  "Experiment",
  "RandomExperiment",
  "FFELogger",
  "ringCardinality"
}

FiniteFieldExperimentsProtect = ()->
(
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

  protect coefficientRingCardinality;
  protect pointLists;
  protect pointsByKey;
  protect countData;
  protect setIsInteresting;
  protect isInteresting;

  protect getExperimentData;
  protect setRecordedProperties;
  protect recordProperty;
  protect ignoreProperty;
  protect ignoreProperties;
  protect observedProperties;

  protect update;
  protect updateExperiment;

  protect propertyList;
  protect clear;
  protect jacobianAtKey;
  protect watchProperties;
  protect recordProperties;

  protect propertyName;
  protect propertyAt;
  protect recordedProperties;
  protect watchedProperties;
  protect useProperties;
  protect useJacobianAt;
  protect usedJacobianAt;

  protect minPointsPerComponent;
  protect setMinPointsPerComponent;
  protect stratificationIntervalView;
  protect countsByCount;  protect countsByCount;

 
  protect  estimateStratification2;
  protect experimentData; 
  protect isRandom;
  protect compatible;
  protect createExperimentData;
)

 

FiniteFieldExperimentsExport  = ()->
(
  
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
 
  exportMutable(coefficientRingCardinality);
  exportMutable( pointLists );
  exportMutable( pointsByKey );

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

 exportMutable( propertyList);
 exportMutable( clear);
 exportMutable( jacobianAtKey);
 exportMutable( watchProperties );
 exportMutable( recordProperties );
 exportMutable( propertyName);
 exportMutable( propertyAt);
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


FiniteFieldExperimentsExport();




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
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert ( (poissonEstimate(16,"confidence" => 2)).min == 8  );
  assert ( (poissonEstimate(16,"confidence" => 2)).max == 24 );
///



 

TEST ///
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

Experiment = new Type of HashTable;

ExperimentData = new Type of MutableHashTable;

--new ExperimentData from HashTable := (E,coeffRing) -> (
new ExperimentData from Ring := (E,coeffRing) -> (
     print (toString E);
     e := new MutableHashTable;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable;
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
   | toExternalString ed.coefficientRing | ","
   | toString (new HashTable from ed.points) | ","
   | toExternalString ed.count | ","
   | toExternalString ed.trials | ","
   | toExternalString ed.propertyList | ","
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





createRandomPointIterator = (numVariables, coeffRing, trials )->
    (
        rpi := new MutableHashTable;
        randomPoint := null;
  
        currTrial := 0;

        assert( trials >= 0 );

        rpi.next=()->
        (
            if (currTrial+1 > trials) then return false;
            randomPoint = random( coeffRing^1, coeffRing^numVariables );         
            currTrial = currTrial + 1;
 
            return true;
        );

        rpi.begin=()->
        (
            ri := createRandomPointIterator(numVariables, trials);
            ri.next();
            return ri;
        );

        rpi.point=()-> randomPoint;

       return new HashTable from rpi;
    );

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  pointIterator = createRandomPointIterator(5,ZZ/7,100)
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
     print "estimated codim <= {wachtched properties}";
     print "--";
     apply(experiment.countsByCount(),i->(
	       --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
	       print (net(round(1,(log(trials)-log(i#0))/log(orderK)))|" <= "|net(i#1)))
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



--new RandomExperiment from BlackBoxParameterSpace := (E, pBlackBoxIdeal) -> 

new Experiment from BlackBoxParameterSpace := (E, pBlackBoxIdeal) -> 
 
(


   blackBoxIdeal := pBlackBoxIdeal;    
   experimentData := new ExperimentData from blackBoxIdeal.coefficientRing;

   coefficientRingCardinality := ringCardinality( blackBoxIdeal.coefficientRing );


   experiment := new MutableHashTable;


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
              bb1.imageRank()             ==   bb2.imageRank() and
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
       if experimentData.trials=!=0 then error ("cannot change jacobianAt  - experiment was already run! You could clear() the statistics and retry. ");

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
    -- count this trial even if no solution is found
    );

   

 
     runRandomExperiment := method();
     runRandomExperiment(ExperimentData, ZZ) := Tally => (experimentData, trials) -> 
     ( 
      
       K := experimentData.coefficientRing;
       --prime = char K
       numVariables := blackBoxIdeal.numVariables;

       apply(trials, trial-> (  randomPoint := random( K^1,  K^numVariables );     
           experimentData.trials = experimentData.trials + 1;    --at first increase trial num, otherwise statistic computations fail ( because trials==0 or reduced by 1) !
           runExperimentOnce( experimentData, randomPoint, minPointsPerComponent );
           
           )
       );
     );

     runExperiment := method();
     --rpi := createRandomPointIterator ( blackBoxIdeal.numVariables, experimentData.coefficientRing, trials );
     -- too slow:
     runExperiment(ExperimentData, Thing) := Tally => (experimentData, pointIterator) -> 
     ( 
       while pointIterator.next() do
       (
           experimentData.trials = experimentData.trials + 1; --at first increase trial num, otherwise statistic computations fail ( because trials==0 or reduced by 1)!
           runExperimentOnce( experimentData, pointIterator.point(), minPointsPerComponent );

       );
     );


     update := method();
     update(ExperimentData) := Tally => opts -> (experimentData) -> 
     ( 
        pointIterator := createPointIterator (  experiment.points() );
        --experimentData.trials = 0;
        experimentData.points = new MutableHashTable;
        experimentData.count = new Tally;

        while ( pointIterator.next() ) do
        (
            runExperimentOnce( experimentData, pointIterator.point(), minPointsPerComponent );
        );
     );


   -- add 'setRecordedProperties' (reinitializes observed properties)  and 'recordProperty' (adds a property to observe)


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
      if experimentData.trials=!=0 then error ("cannot change watched properties - experiment was already run! You could clear() the statistics and retry.");

      setRecordedPropertiesInternal(propertyStringList);
      update(experimentData);
   );

   UpdateRecordedPropertiesError := "cannot change watched properties - experiment was already run! You could clear() the statistics and retry.";
  

   experiment.watchProperties = experiment.setRecordedProperties;

   experiment.recordProperty=(propertyName)->
   (
       if experimentData.trials=!=0 then error (UpdateRecordedPropertiesError);

          if not blackBoxIdeal.hasPointProperty(propertyName) then 
              error ("blackBoxIdeal seems not to have property" | propertyName );

      experimentData.propertyList = unique (experimentData.propertyList | { propertyName }) ;
      setRecordedPropertiesInternal( experimentData.propertyList );   
      update(experimentData);
   );

   experiment.ignoreProperty=(propertyName)->
   (
       if experimentData.trials=!=0 then error (UpdateRecordedPropertiesError);

      experimentData.propertyList = delete(propertyName, experimentData.propertyList ) ;
      setRecordedPropertiesInternal( experimentData.propertyList );   
        
      update(experimentData);
   );


   experiment.ignoreProperties = (ignorePropertyStringList)->
   (
       if experimentData.trials=!=0 then error (UpdateRecordedPropertiesError);

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
   experiment.run(ZZ) := Thing=> (trials)->
   (
      rpi := createRandomPointIterator ( blackBoxIdeal.numVariables, experimentData.coefficientRing, trials );
       runExperiment( experimentData, rpi );
       --runRandomExperiment( experimentData, trials );
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

   experiment.trials = ()->
   (
      return experimentData.trials;
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


   experiment = newClass( Experiment, experiment );
   ffelog.debug ("type of internal experiment variable is : " | toString class experiment );
   return experiment; 
   -- return new HashTable from experiment; 
);

 


-- tryProperty = (experiment, property) ->(
--     pointLists := experiment.pointLists();
--     apply(
--     	  apply((keys pointLists),key -> (key,tally apply(pointLists#key,property))),
--     	  i->print (net i#0|" => "|net i#1)
--     	  )
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
            Creates an experiment from a black box. The experiment object will keep
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
	   \break If a black box hat a property "rankJacobianAt" it is
	   automatically watched:
	Example
	   e.watchedProperties() 
	Text
	   \break Now we run the experiment by evaluating at 1250 random points:  
	Example
	   time e.run(1250)      
	Text
	   \break Let's see what kind of points were found:
	Example
	   e.countData() 
	Text
	   There are 125 points in (F_5)^3 of which 25 are on the
	   plane, 5 are on the line and 1 (the origin) is on the line and the plane.
	   We therefore expect about 10 points with rankJacobian 0, 240 with rankJacobian 1
	   and 40 with rankJacobian 2.
	   
	   Often it is useful to have explicit points on all components of a variety. Therefore
	   the experiment has stored some of the points:
	Example
	   e.pointLists()
	Text
	   Since one always finds many points found on components of low
	   codimension it is not useful to remember all of them. The experiment
	   remembers only about 10 points per component:
	Example
	   apply(keys e.pointLists(),key->print (key => #((e.pointLists())#key))); 
	Text
	   Here we have not collected exactly 10 points per component since
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

 
    bbRankM = bbRankM.getUpdatedBlackBox()
 

    e = new RandomExperiment from bbRankM
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

    e.stratificationIntervalView() -- test fails here 

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

