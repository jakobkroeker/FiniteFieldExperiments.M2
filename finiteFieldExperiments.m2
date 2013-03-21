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
    return new Interval from { (symbol min)=>seq#0, (symbol max)=>seq#1};
);

TEST ///
new Interval from { (symbol min)=>0, (symbol max)=>1};
new Interval from (1,1);
///

net (Interval) := Net =>(interval)->
(
   str := "(" | toString interval.min | ", " | toString interval.max | ")";
   return net str;
)

ExperimentData = new Type of MutableHashTable;

-- new ExperimentData from BlackBoxIdeal := 
new ExperimentData from HashTable := (E,blackBoxIdeal) -> (
     print (toString E);
     e := new MutableHashTable;
     e.blackBoxIdeal = blackBoxIdeal;
     e.coefficientRing = blackBoxIdeal.coefficientRing();
     e.points = new MutableHashTable;
     e.count = new Tally;
     e.trials = 0;
     return e;
     )
     

 

-- prime is fieldCardinality?
estimateNumberOfComponents = method(Options => {"confidence" => 1.96});
estimateNumberOfComponents(ZZ,ZZ,ZZ,ZZ) := HashTable => opts->( trials,rankJacobian,pointCount, prime) -> (
     confidence := opts#"confidence";
     expected := pointCount*prime^rankJacobian/trials*1.0;
     err := confidence*sqrt(pointCount)*prime^rankJacobian/trials*1.0;
     estimate := new Interval from (round(1,expected-err),round(1,expected+err));
     return estimate;
     )

TEST ///
  estimateNumberOfComponents(1000,1,100,11)
  estimateNumberOfComponents(1000,1,100,11,"confidence"=>1.0)
///


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
   experimentData := new ExperimentData from B;

   experimentStatistics := experimentData -> (
       count = experimentData.count;
       new Tally from apply( keys count, rankJacobian->rankJacobian => estimateNumberOfComponents(
	         experimentData.trials,
	         rankJacobian,
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
     bb := experimentData.blackBoxIdeal;
     numVariables := bb.numVariables();
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

  

   experiment.runExperiment= method(Options => options runExperimentOnce);
   experiment.runExperiment(ZZ) := Thing=> opts->(trials)->
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

e.runExperiment( 2000,  "numPointsPerComponentToCollect"=>20 )

pointData = e.getPointData()

apply(keys pointData,i->#(pointData#i))




