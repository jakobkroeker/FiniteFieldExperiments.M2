-- finite field experiments




Interval = new Type of HashTable;

new Interval from HashTable := (E,ht) -> 
(  
  assert(ht?#(symbol min) );
  assert(ht?#(symbol max) );
  assert(ht.min<=ht.max);
  return ht;
);

-- maybe Interval coould also contain statistical information?

new Interval from Sequence := (E,seq) -> 
(
    assert(#seq==2);
    return new Interval from { (symbol min)=>seq#0*1.0, (symbol max)=>seq#1*1.0};
);

TEST ///
new Interval from { (symbol min)=>0, (symbol max)=>1};
new Interval from (1,1);
///

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
     estimate := new Interval from (round(1,numPoints-err),round(1,numPoints+err));
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

QQ * Interval := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

TEST ///
   assert (2*(new Interval from (1,2)) == new Interval from (2,4))
   assert (2.0*(new Interval from (1,2)) == new Interval from (2,4))
///

Interval == Interval := (I1,I2) -> (
     (I1.min == I2.min) and (I1.max == I2.max)
     )

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
    	 logEst := new Interval from (log est.max,log est.min);
    	 codimEst := (-1/log fieldCardinality)*logEst;
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

        rpi.next=()->
        (
            if (currTrial==trials) then error " cannot get next point";
            randomPoint = random( coeffRing^1, coeffRing^numVariables );         
            currTrial = currTrial + 1;
            return (currTrial<trials) ;
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

Experiment = new Type of HashTable;

new Experiment from HashTable := (E,blackBoxIdeal) -> 
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



   runExperimentOnce := method(Options => {"numPointsPerComponentToCollect" => 10});
 
   runExperimentOnce(ExperimentData) := Tally => opts -> experimentData -> ( 
       numberOfPointsToCollect := opts#"numPointsPerComponentToCollect";	   
       K := experimentData.coefficientRing;
       --prime = char K
       numVariables := blackBoxIdeal.numVariables();
       -- todo 'numVariables' keine Funkton - in Konflikt mit "Vermeidung von doppeltem code"

     -- choose a random point
     randomPoint := random(K^1,K^numVariables);
     -- if ideal vanishes on the random point do something
     if bb.isZeroAt(randomPoint) then   -- only if jacobianAt is defined
     (
     	  rankJacobian := rank bb.jacobianAt(randomPoint);
	  
      -- count number of found points for each rank (=codim Tangentspace)
     	  experimentData.count = experimentData.count + tally {rankJacobian};
	 
      -- remember all points
     	  if experimentData.points#?rankJacobian then 
          (
	        estimate = (experimentStatistics(experimentData))#rankJacobian;
	        -- collect a fixed number of points per estimated component
	        -- use upper limit of estimation for this
	        if #(experimentData.points#rankJacobian) < numberOfPointsToCollect*max(1,estimate.max) then 
            (
	       	    experimentData.points#rankJacobian = experimentData.points#rankJacobian | {randomPoint}
	       	);
	      ) else (
	           experimentData.points#rankJacobian = {randomPoint}
     	  );
	  );
      -- count this trial even if no solution is found
      experimentData.trials = experimentData.trials + 1;
      return experimentStatistics (experimentData);
    );

    runExperiment := method(Options => options runExperimentOnce);

     runExperiment(ExperimentData,ZZ) := Tally => opts -> (experimentData,trials) -> ( 
       apply( trials, trialNum->runExperimentOnce( experimentData, opts) );
       return experimentStatistics (experimentData);
     );

   ---- syntax for the moment too hard (method without parameters)
   --experiment.runExperimentOnce = method(Options => (options runExperimentOnce));
   -- experiment.runExperimentOnce() := Thing => opts->()->
   --(
   --   return runExperimentOnce(experimentData);
   --);

  

   experiment.run = method(Options => options runExperimentOnce);
   experiment.run(ZZ) := Thing=> opts->(trials)->
   (
      return runExperiment(experimentData,trials);
   );
  
   experiment.getPointData = ()->
   (
      return new HashTable from experimentData.points;
   );

   return new HashTable from experiment; 
);



end
---

restart
loadPackage"idealBlackBoxes"
-- besser BlackBoxIdeals
load "finiteFieldExperiments.m2"

R = (ZZ/7)[x_0..x_3]
M = matrix{
     {x_0,x_1,0},
     {x_1,x_2,x_3}
     }
I = minors(2,M)
B = blackBoxIdeal I
 
e = new Experiment from B

e.runExperiment(1, "numPointsPerComponentToCollect"=>20 ) 
e.runExperiment(1) 

e.runExperiment( 2000,  "numPointsPerComponentToCollect"=>20 )

pointData = e.getPointData()

apply(keys pointData,i->#(pointData#i))




