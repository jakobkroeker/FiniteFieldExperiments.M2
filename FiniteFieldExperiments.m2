-- finite field experiments





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





TEST ///
///

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





net (Interval) := Net =>(interval)->
(
   str := "[" | toString interval.min | ", " | toString interval.max | "]";
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
   assert (2*(new Interval from (1,2)) == new Interval from (2,4))
   assert (2.0*(new Interval from (1,2)) == new Interval from (2,4))
///

Interval == Interval := (I1,I2) -> (
     (I1.min == I2.min) and (I1.max == I2.max)
     )

TEST ///
restart

loadPackage"BlackBoxIdeals"
load "FiniteFieldExperiments.m2"
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
  assert (poissonEstimate(16,"confidence"=>2) == new Interval from (8,24))
///

estimateNumberOfComponents = method(Options => (options poissonEstimate));
estimateNumberOfComponents(ZZ,ZZ,ZZ,ZZ) := HashTable => opts->
    (trials,estimatedCodim,numPoints,fieldCardinality) -> (
    return poissonEstimate(numPoints,opts)*(fieldCardinality^estimatedCodim/trials*1.0);
    )

TEST ///
  assert (estimateNumberOfComponents(11^2*16,2,16,11,"confidence"=>2) == new Interval from (0.5,1.5))
///

estimateCodim = method(Options => (options poissonEstimate));
estimateCodim(ZZ,ZZ,ZZ) := HashTable => opts->
    (trials,numPoints,fieldCardinality) -> (
    	 est := (1/trials)*poissonEstimate(numPoints,opts);
    	 logEst := new Interval from (-log est.max,-log est.min);
    	 codimEst := (1/log fieldCardinality)*logEst;
	 new Interval from (round(1,codimEst.min),round(1,codimEst.max))
    )

TEST ///
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
            if (currTrial+1 >= trials) then return false;
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
            if (currTrial+1 >= pointCount) then return false;
            point = points#currPosition;
            currPosition = currPosition + 1;
            return (true) ;
        );

        rpi.begin=()->
        (
            localPointIterator := createPointIterator( pPoints);
            localPointIterator.next();
            return localPointIterator;
        );

        rpi.point=()-> point;

       return new HashTable from pIterator;
    );



Experiment = new Type of HashTable;

new Experiment from HashTable := (E, blackBoxIdeal) -> 
(
   experiment := new MutableHashTable;
   experimentData := new ExperimentData from blackBoxIdeal.coefficientRing();
   -- blackBoxIdeal := pblackBoxIdeal;
   
   estimateDecomposition := experimentData -> (
       count = experimentData.count;
       new Tally from apply( keys count, rankJacobian->
	    rankJacobian => estimateNumberOfComponents(
	         experimentData.trials,
	         rankJacobian,
	         count#rankJacobian,
	         char experimentData.coefficientRing
	       )
	  )
   );

   

   experiment.estimateDecomposition = () -> (estimateDecomposition(experimentData));

   estimateStratification := experimentData -> (
       count = experimentData.count;
       new Tally from apply( keys count, rankJacobian->
	    rankJacobian => estimateCodim(
	         experimentData.trials,
	         count#rankJacobian,
	         char experimentData.coefficientRing
	       )
	  )
   );

 

   experiment.estimateStratification = () -> (estimateStratification(experimentData));
   
   stratificationIntervalView := (stratificationData )->
   (
     
      prerearrangedData := new HashTable from apply( keys stratificationData, key -> (stratificationData#key=>(stratificationData#key,key )));
      toSort := apply( keys stratificationData, key -> (stratificationData#key));
      sorted := sort (toSort); 
      rearrangedData := apply (sorted, key-> (prerearrangedData#key) );
      return new HashTable from rearrangedData;
   );

   experiment.stratificationIntervalView = ()->
   (
        stratificationData := experiment.estimateStratification();
        return stratificationIntervalView(stratificationData);
   );


   isInteresting := ()->true;
   
   if  blackBoxIdeal#?(getSymbol "isZeroAt") then
      isInteresting =blackBoxIdeal.IsZeroAt;
    
   experiment.isInteresting=(point)->
   ( 
      return isInteresting(point);
   );

   experiment.getPoints=()->
   (
      return flatten  apply (keys experimentData.points, countKey-> experimentData.points#countKey);
   );



   runExperimentOnce := method();

   numPointsPerComponentToCollect := 10;

   experiment.setNumPointsPerComponentToCollect = (numPointsPerComponent)->
   ( 
     numPointsPerComponentToCollect = numPointsPerComponent;
   );

   experiment.numPointsPerComponentToCollect = ()->
   ( 
     return numPointsPerComponentToCollect ;
   );
 
   runExperimentOnce(ExperimentData, Matrix,ZZ ) := Tally => (experimentData, point, wantedPoints) -> (  
       K := experimentData.coefficientRing;
       --prime = char K
       numVariables := blackBoxIdeal.numVariables();
       -- todo 'numVariables' keine Funkton - in Konflikt mit "Vermeidung von doppeltem code"

       -- choose a random point
    
       -- if ideal vanishes on the random point do something
       if blackBoxIdeal.isZeroAt(point) then   
       (
	  countKey := {};  
	  if blackBoxIdeal#?(symbol jacobianAt) then (
	       countKey = countKey|{rank blackBoxIdeal.jacobianAt(point)}
	       );
	  if blackBoxIdeal#?(symbol propertiesAt) then (
	       countKey = countKey|{blackBoxIdeal.propertiesAt(point)}
	       );	  
	  
          -- count number of found points for each rank and property
     	  experimentData.count = experimentData.count + tally {countKey};
	 
      	  -- remember some points
	  if blackBoxIdeal#?(symbol jaobianAt) then (
	        upperEstimate := ((estimateDecomposition(experimentData))#countKey).max;
		wantedPoints = max(1,upperEstimate)*wantedPoints;
		); 	       
     	  if experimentData.points#?countKey then (
	        -- collect a fixed number of points per estimated component
	        -- use upper limit of estimation for this
	        if #(experimentData.points#countKey) < wantedPoints then 
            	(
	       	    experimentData.points#countKey = experimentData.points#countKey | {randomPoint}
	       	);
	   ) else (
	        experimentData.points#countKey = {randomPoint}
     	  	);
      );
      -- count this trial even if no solution is found
    );

   

 
     runRandomExperiment := method();
     runRandomExperiment(ExperimentData, ZZ) := Tally => (experimentData, trials) -> 
     ( 
      
       K := experimentData.coefficientRing;
       --prime = char K
       numVariables := blackBoxIdeal.numVariables();

       apply(trials, trial-> (  randomPoint := random( K^1,  K^numVariables );         
           runExperimentOnce( experimentData, randomPoint, numPointsPerComponentToCollect );
           experimentData.trials = experimentData.trials + 1;)
       );
     );

     runExperiment := method();
     --rpi := createRandomPointIterator ( blackBoxIdeal.numVariables(), experimentData.coefficientRing, trials );
     -- too slow:
     runExperiment(ExperimentData, Thing) := Tally => (experimentData, pointIterator) -> 
     ( 
       while pointIterator.next() do
       (
           runExperimentOnce( experimentData, pointIterator.point(), numPointsPerComponentToCollect );
           experimentData.trials = experimentData.trials + 1;
       );
     );

     updateExperiment = method();
     updateExperiment(ExperimentData) := Tally => opts -> (experimentData) -> 
     ( 
        pointIterator := createPointIterator (  experiment.getPoints() );
        experimentData.trials = 0;
        experimentData.points = new MutableHashTable;
        experimentData.count = new Tally;

        while ( pointIterator.next() ) do
        (
            runExperimentOnce( experimentData, pointIterator.point(), numPointsPerComponentToCollect );
        );
     );

   experiment.setIsInteresting = (pIsInteresting)->
   (  
      if ( pIsInteresting=!=isInteresting ) then
      (
         isInteresting = pIsInteresting;
         updateExperiment(experimentData);
      );
   );

   experiment.setIdealPropertiesAt =(idealPropertiesAt)->
   (
      blackBoxIdeal.setPropertiesAt( idealPropertiesAt );
      updateExperiment(experimentData);
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
      rpi := createRandomPointIterator ( blackBoxIdeal.numVariables(), experimentData.coefficientRing, trials );
      return  runExperiment( experimentData, rpi );
      --return runRandomExperiment( experimentData, trials );
   );
  
   experiment.getPointData = ()->
   (
      return new HashTable from experimentData.points;
   );
 
   experiment.getCountData = ()->
   (
      return new Tally from experimentData.count;
   );

   experiment.getTrials = ()->
   (
      return experimentData.trials;
   );


   


   return new HashTable from experiment; 
);



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

pointData = e.getPointData()

apply(keys pointData,i->#(pointData#i))

assert (#(pointData#{2}) >= 40)

