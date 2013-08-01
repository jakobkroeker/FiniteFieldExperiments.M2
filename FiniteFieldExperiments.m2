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
     PackageExports => {"BlackBoxIdeals"},
     Headline => "black boxes for explicit and implicitly given ideals",
     DebuggingMode => true
)

needsPackage "BlackBoxIdeals";


export {
  "Interval",
  "estimateDecomposition", 
  "estimateStratification",
  "estimateCodim", 
  "estimateNumberOfComponents",    	    	 
  "poissonEstimate",
  "Experiment",
  "FFELogger"
}

FiniteFieldExperimentsProtect = ()->
(
  protect center;
  protect next;
  protect begin;
  protect count;
  protect point;
  protect collectedCount;
  protect pointKeys; 
  protect points;
  protect trials;

  protect pointLists;
  protect countData;
  protect setIsInteresting;
  protect isInteresting;

  protect getExperimentData;
  protect setObservedProperties;
  protect observeProperty;
  protect ignoreProperty;
  protect ignoreProperties;
  protect observedProperties;

  protect update;
  protect updateExperiment;

  protect propertyList;
  protect clear;
  protect jacobianAtKey;
  protect watchProperties;
  protect propertyName;
  protect propertyAt;
  protect watchedProperties;
  protect useProperties;
  protect useJacobianAt;
  protect usedJacobianAt;

  protect minPointsPerComponent;
  protect setMinPointsPerComponent;
  protect stratificationIntervalView;
  protect countsByCount;  protect countsByCount;

  protect setIdealPropertiesAt;
 protect  estimateStratification2;
  
)

 

FiniteFieldExperimentsExport  = ()->
(
  exportMutable( center);
  exportMutable( next);
  exportMutable( begin);
  exportMutable( count);
  exportMutable( point);
  exportMutable( collectedCount);
  exportMutable( pointKeys);
  exportMutable( points);
  exportMutable( trials );
 
  exportMutable( pointLists );

  exportMutable( countData );

  exportMutable( setIsInteresting);
  exportMutable( isInteresting);

  exportMutable( getExperimentData);
  exportMutable( setObservedProperties);
  exportMutable( observeProperty);
  exportMutable( ignoreProperty);
  exportMutable( ignoreProperties);
  exportMutable( observedProperties);


  exportMutable( update);
  exportMutable( updateExperiment);

 exportMutable( propertyList);
 exportMutable( clear);
 exportMutable( jacobianAtKey);
 exportMutable( watchProperties);
 exportMutable( propertyName);
 exportMutable( propertyAt);
 exportMutable( watchedProperties);
 exportMutable( useProperties);
 exportMutable( useJacobianAt);
 exportMutable( usedJacobianAt);
 

  exportMutable( minPointsPerComponent);
  exportMutable( setMinPointsPerComponent);
  exportMutable( stratificationIntervalView);

   exportMutable( setIdealPropertiesAt);
  exportMutable( countsByCount);

 exportMutable (  estimateStratification2);

)


FiniteFieldExperimentsExport();

FFELogger = createLogger();




Interval = new Type of HashTable;

roundInterval = method();


roundInterval ( HashTable, ZZ) := Interval=>(I,n)->
(
     return new Interval from ( round(n,I.min), round(n,I.max) );
)

intervalCenter = method();

intervalCenter (HashTable) := Interval=>(I)->
(
    return (I.min+I.max)/2.0;
)


new Interval from Thing := (E,ll) -> 
(
   error ("this interface is not implemented");
);

--new Interval from HashTable := (E,ll) -> 
new Interval from List := (E,ll) -> 
( 
  ht := new MutableHashTable from ll;
  --print "new Interval from HashTable";
  assert(ht#?(symbol min) );
  assert(ht#?(symbol max) );
  assert(ht.min<=ht.max);
 
  ht.round = method();
  ht.round(ZZ) := Interval => (n)->
  (
     return roundInterval( ht, n);
  );

  ht.center = intervalCenter(ht);

  --print( keys ht);
  return new HashTable from ht;
);

-- maybe Interval coould also contain statistical information?

new Interval from Sequence := (E,seq) -> 
(
    --print "new Interval from Sequence";
    assert(#seq==2);
    ht := { (symbol min)=>seq#0*1.0, (symbol max)=>seq#1*1.0};
    return new Interval from ht;
);



formatIntervalBounds = (num)->
(
  f :=(pnum) ->format(10,2,2,2,pnum);
  str:="";
  if (value f num)==0 then return "0.";
  if (abs(value f num)<1) then 
  str = str|"0";
  origstr := format(10,2,2,2,num);
  if origstr#0=="-" then
  (
    origstr = origstr_(1..#origstr-1);
     str = "-"|str;
  );
  str = str|origstr;
  return str;
)


net (Interval) := Net =>(interval)->
(
   str :=  " " | formatIntervalBounds  interval.center  | " " | "[" | formatIntervalBounds  interval.min | ", " | formatIntervalBounds interval.max | "]";
   return net str;
)

poissonEstimate = method(Options => {"confidence" => 1.96});
poissonEstimate(ZZ) := HashTable => opts -> (numPoints) -> ( 
     -- we use a poisson approximation since there are only 
     -- very few solution points compared to the number of trials
     --
     -- we then use the normal approximatiom of the poissondistribution
     -- to calculate the confidence interval
     --
     -- for a Poisson distribution the standard deviation is
     -- the square root of the expected value
     err := opts#"confidence"*sqrt(numPoints); 
     estimate := new Interval from (max(round(1,numPoints-err),0.0001),round(1,numPoints+err));
     return estimate;
     )

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert ((poissonEstimate(16,"confidence" => 2)).min == 8);
  assert ((poissonEstimate(16,"confidence" => 2)).max == 24);
///

RR * Interval := (scalingFactor,I) -> (
     new Interval from (scalingFactor*I.min,scalingFactor*I.max)
     )

Interval * RR := (I,scalingFactor) -> (
     scalingFactor*I
     )

ZZ * Interval := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

Interval * ZZ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )

QQ * Interval := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

Interval * QQ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )

-- problem: sortierung nicht ohne weiteres austauschbar -> wrappertyp muss eingeführt werden, der eine andere sortierung implementiert
-- und alle objekte müssen in den Wrappertyp konvertiert werden.

Interval ? Interval := (I1,I2) ->
(
   if  (I1.center==I2.center ) then 
   (
        I2.min ? I1.min 
   )
   else 
   (  
         I1.center ? I2.center
   )
)



TEST ///
   debug FiniteFieldExperiments
   FiniteFieldExperimentsProtect()
   assert (2*(new Interval from (1,2)) == new Interval from (2,4))
   assert (2.0*(new Interval from (1,2)) == new Interval from (2,4))
///

Interval == Interval := (I1,I2) -> (
     (I1.min == I2.min) and (I1.max == I2.max)
     )

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  new Interval from { (symbol min)=>0, (symbol max)=>1};
  i2 := new Interval from (1, 1.04);
  i3 := new Interval from (1, 1.045);
  assert( i3.round(2)== i2);
  i4 :=  new Interval from (2.5,3.5);
  assert(i4.center == 3.0);
  i2<i3;
  i2?i3;
  assert ( (i2?i3) === (1?2));
  assert ( (i2?i2) === (1?1));
  i2?i2;
///


TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert (poissonEstimate(16,"confidence"=>2) == new Interval from (8,24))
///

estimateNumberOfComponents = method(Options => (options poissonEstimate));
estimateNumberOfComponents(ZZ,RR,ZZ,ZZ) := HashTable => opts->
    (trials,estimatedCodim,numPoints,fieldCardinality) -> (
    return poissonEstimate(numPoints,opts)*(fieldCardinality^estimatedCodim/trials*1.0);
    )
estimateNumberOfComponents(ZZ,ZZ,ZZ,ZZ) := HashTable => opts->    (trials,estimatedCodim,numPoints,fieldCardinality) -> 
(
    estimateNumberOfComponents(trials,estimatedCodim*1.0,numPoints,fieldCardinality)
)

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  estimate =  estimateNumberOfComponents(11^2*16,2,16,11,"confidence"=>2);
  estimate = estimate.round(1);
  assert ( estimate == new Interval from (0.5,1.5))
///

estimateCodim = method(Options => (options poissonEstimate));
estimateCodim(ZZ,ZZ,ZZ) := Interval => opts->
    ( trials, numPoints, fieldCardinality ) -> (
    	 est := (1/trials)*poissonEstimate(numPoints,opts);
    	 logEst := new Interval from (-log est.max,-log est.min);
    	 codimEst := (1/log fieldCardinality)*logEst;
	 new Interval from (round(1,codimEst.min),round(1,codimEst.max))
    )

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  assert (estimateCodim(11^2*10,10,11) == new Interval from (1.8,2.4))
///

ExperimentData = new Type of MutableHashTable;


new ExperimentData from HashTable := (E,coeffRing) -> (
     print (toString E);
     e := new MutableHashTable;
     --e.blackBoxIdeal = blackBoxIdeal;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable;
     e.count = new Tally;
     e.trials = 0;
     e.propertyList = {};
     return e;
     )

     






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


estimateDecomposition := (experimentData) -> (
       count := experimentData.count;
       cardinality := null;
       if char experimentData.coefficientRing ==0 then   cardinality  = infinity;
       try { cardinality = experimentData.coefficientRing.order } then () else();
       new Tally from apply( keys count, countKey->
	    countKey => estimateNumberOfComponents(
	         experimentData.trials,
	         (estimateCodim( experimentData.trials,  count#countKey, cardinality )).center,
	         count#countKey, 
	         cardinality
	       )
	  )
   );


estimateStratification2 := (experiment) -> (
     trials := experiment.trials();
     orderK := (experiment.coefficientRing()).order; -- this must be read from the experimentdata
     print "--";
     apply(experiment.countsByCount(),i->(
	       --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
	       print (net(round(1,(log(trials)-log(i#0))/log(orderK)))|" <= "|net(i#1)))
	       )
	  ;
     print "--";
     )



  estimateStratification := experimentData -> (
       count := experimentData.count;
       cardinality := null;
       if char experimentData.coefficientRing ==0 then   cardinality  = infinity;
       try { cardinality = experimentData.coefficientRing.order } then () else();

       new Tally from apply( keys count, countKey->
	    countKey => estimateCodim(
	         experimentData.trials,
	         count#countKey,
	         cardinality
	       )
	  )
   );

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

      toSort := apply( keys counts, key -> (counts#key));
      sorted := sort (toSort); 
      rearrangedData := {};
      for count in sorted do
      (
          for entry in prerearrangedData#count do
          rearrangedData = rearrangedData | {(count,entry)};
      );
     
      return new List from rearrangedData;
   );

Experiment = new Type of HashTable;



new Experiment from HashTable := (E, pBlackBoxIdeal) -> 
(


   blackBoxIdeal := pBlackBoxIdeal;    
   experimentData := new ExperimentData from blackBoxIdeal.coefficientRing;
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

   experiment.estimateDecomposition = () -> (estimateDecomposition(experimentData));


   experiment.estimateStratification = () -> (estimateStratification(experimentData));

   experiment.estimateStratification2 = () -> ( estimateStratification2(experiment) );

   
   
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
                FFELogger.log(4, "update wanted points");
                rankJacobian := rank  jacobianAt(point); 
                upperEstimate := ((estimateDecomposition( experimentData ))#countKey).max;
                wantedPoints = max(1,upperEstimate)*wantedPoints;
            ); 	       

          -- remember some points
            if experimentData.points#?(countKey) then 
            (
                -- collect a fixed number of points per estimated component
                -- use upper limit of estimation for this
                if #(experimentData.points#countKey) < wantedPoints then 
                (
                    FFELogger.log(4, "attaching point");
                    experimentData.points#countKey = experimentData.points#countKey | {point}
                );
            )
            else (
                FFELogger.log(4, "attaching first point for some key");
                experimentData.points#countKey = {point}
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
           experimentData.trials = experimentData.trials + 1;    --at first increase trial num!
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
           experimentData.trials = experimentData.trials + 1; --at first increase trial num!
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

   -- add 'setObservedProperties' (reinitializes observed properties)  and 'observeProperty' (adds a property to observe)

   setObservedPropertiesInternal := (propListToObserve)->
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

   experiment.setObservedProperties = ( propertyStringList )->
   (  
      if experimentData.trials=!=0 then error ("cannot change watched properties - experiment was already run! You could clear() the statistics and retry.");

      --if jacobianAtKey =!=null then 
      --propertyStringList = { jacobianAtKey } | propertyStringList;

      for propertyName in propertyStringList do
      (
          if not blackBoxIdeal.hasPointProperty(propertyName) then 
              error ("blackBoxIdeal seems not to have property" | propertyName );
      );

      setObservedPropertiesInternal(propertyStringList);
      update(experimentData);
   );


  

   experiment.watchProperties=experiment.setObservedProperties;

   experiment.observeProperty=(propertyName)->
   (
       if experimentData.trials=!=0 then error ("cannot change watched properties - experiment was already run! You could clear() the statistics and retry.");

          if not blackBoxIdeal.hasPointProperty(propertyName) then 
              error ("blackBoxIdeal seems not to have property" | propertyName );

      experimentData.propertyList = unique (experimentData.propertyList | { propertyName }) ;
      setObservedPropertiesInternal( experimentData.propertyList );   
      update(experimentData);
   );

   experiment.ignoreProperty=(propertyName)->
   (
       if experimentData.trials=!=0 then error ("cannot change observed properties - experiment was already run! You could clear() the statistics and retry.");

      experimentData.propertyList = delete(propertyName, experimentData.propertyList ) ;
      setObservedPropertiesInternal( experimentData.propertyList );   
        
      update(experimentData);
   );

   experiment.ignoreProperties = (ignorePropertyStringList)->
   (
       if experimentData.trials=!=0 then error ("cannot change observed properties - experiment was already run! You could clear() the statistics and retry.");

      apply( ignorePropertyStringList, propToIgnore-> ( experimentData.propertyList = delete(propToIgnore, experimentData.propertyList ); ));
      setObservedPropertiesInternal( experimentData.propertyList );    
      update(experimentData);
   );

   -- maybe observed Properties is not a good name - is 'recordedProperties' better ?
   experiment.observedProperties =  ()->
   (
      return experimentData.propertyList;
   );

   experiment.watchedProperties = experiment.observedProperties;

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
      return  runExperiment( experimentData, rpi );
      --return runRandomExperiment( experimentData, trials );
   );
  
   experiment.pointLists = ()->
   (
      return new HashTable from experimentData.points;
   );
 
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
   

   --init:
   -- todo: test if changing blackBoxIdeal.jacobianAt is transparent (means after blackBoxIdeal.updatePointProperty("jacobianAt") the new one is called)
 
    experiment.watchProperties( {} );

       if  blackBoxIdeal.hasPointProperty("isZeroAt") then
       (
          --print("isInteresting"| toString  (blackBoxIdeal.isZeroAt ));
          isInteresting = blackBoxIdeal.pointProperty("isZeroAt");

               experiment.watchProperties( {"isZeroAt"} );
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
            experiment.watchProperties( {"isZeroAt","rankJacobianAt"} );
       );

   --end init:


   return new HashTable from experiment; 
);


-- tryProperty = (experiment, property) ->(
--     pointLists := experiment.pointLists();
--     apply(
--     	  apply((keys pointLists),key -> (key,tally apply(pointLists#key,property))),
--     	  i->print (net i#0|" => "|net i#1)
--     	  )
--     )


TEST ///
    debug FiniteFieldExperiments
    FiniteFieldExperimentsProtect()
    coeffRing := ZZ/3;
    bbRankM = blackBoxIdeal( 5 ,coeffRing )
    rankMat := (point)->5
    bbRankM.registerPointProperty("rankMat",rankMat)


    point := matrix {{1,2,3,4,5}};
    point = sub( point, coeffRing);

    bbRankM = bbRankM.rebuild()

    e = new Experiment from bbRankM
    assert (e.coefficientRing()===coeffRing);

    e.setMinPointsPerComponent(20);
    assert( e.minPointsPerComponent()==20);
    FFELogger.setLevel(4);
    e.watchProperties {"rankMat"};
    e.watchedProperties()
    assert( 1== # select( e.watchedProperties() , 
                       (prop)->(prop=="rankMat") ) 
     )
    e.useJacobianAt("rankMat");
    e.useJacobianAt(null);
   
    e.countsByCount()
    points := e.points();
    #points
    apply(points,point->rankMat(point))
    FFELogger.setLevel(2);
    time e.run(1000)
    assert (e.trials()==1000);

    e.estimateStratification()
    e.estimateDecomposition()
    e.stratificationIntervalView()
    e.collectedCount()
    e.watchedProperties()
    e.jacobianAtKey()

    bbRankM.knownPointProperties()
///


end
---

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

