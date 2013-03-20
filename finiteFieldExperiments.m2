-- finite field experiments

R = (ZZ/7)[x_0..x_3]
M = matrix{
     {x_0,x_1,0},
     {x_1,x_2,x_3}
     }
I = minors(2,M)
B = blackBoxIdeal I
 

Experiment = new Type of MutableHashTable 

-- new Experiment from BlackBoxIdeal := 
new Experiment from HashTable := (E,blackBoxIdeal) -> (
     e := new MutableHashTable;
     e.blackBoxIdeal = blackBoxIdeal;
     e.coefficientRing = blackBoxIdeal.coefficientRing();
     e.points = new MutableHashTable;
     e.count = new Tally;
     e.trials = 0;
     return e;
     )
     
    
e = new Experiment from B

estimateNumberOfComponents = (trials,rankJacobian,count,prime) -> (
     confidence := 1.96; -- could be an option
     expected := count*prime^rankJacobian/trials*1.0;
     err := confidence*sqrt(count)*prime^rankJacobian/trials*1.0;
     (round(1,expected-err),round(1,expected+err))
     )

experimentStatistics = experiment -> (
     count = experiment.count;
     new Tally from apply(keys count,i->i => estimateNumberOfComponents(
	       e.trials,
	       i,
	       count#i,
	       char experiment.coefficientRing
	       )
	  )
     )

runExperimentOnce = method(Options => {numberOfPointsToCollect => 10});
 
runExperimentOnce(Experiment) := Tally => opts -> experiment -> ( 
     numberOfPointsToCollect = opts.numberOfPointsToCollect;
     -- this should be an option;	   
     K := experiment.coefficientRing;
     --prime = char K
     bb := experiment.blackBoxIdeal;
     numVariables := bb.sourceRank();
     -- besser keine Funkton
     -- besser: numVariables
     --randomPoints = select(apply(100,i->random(K^1,K^numVariables)),bb.isZeroAt)
     -- choose a random point
     randomPoint := random(K^1,K^numVariables);
     -- if ideal vanishes on the random point do something
     if bb.isZeroAt(randomPoint) then (
     	  -- only if jacobianAt is defined
     	  rankJacobian = rank bb.jacobianAt(randomPoint);
	  -- count number of found points for each rank (=codim Tangentspace)
     	  experiment.count = experiment.count + tally {rankJacobian};
	  -- remember all points
     	  if experiment.points#?rankJacobian then (
	       estimate = (experimentStatistics(experiment))#rankJacobian;
	       -- collect a fixed number of points per estimated component
	       -- use upper limit of estimation for this
	       if #(experiment.points#rankJacobian) < numberOfPointsToCollect*max(1,estimate#1) then (
	       	    experiment.points#rankJacobian = experiment.points#rankJacobian | {randomPoint}
	       	    );
	       ) else (
	       experiment.points#rankJacobian = {randomPoint}
     	       );
	  );
      -- count this trial even if no solution is found
      experiment.trials = experiment.trials + 1;
      experimentStatistics experiment
      )

runExperiment = method(Options => options runExperimentOnce);

runExperiment(Experiment,ZZ) := Tally => opts -> (experiment,trials) -> ( 
--runExperiment = (experiment,trials) -> (
     apply(trials,i->runExperimentOnce(experiment,opts));
     experimentStatistics experiment
     )


end
---

restart
loadPackage"idealBlackBoxes"
-- besser BlackBoxIdeals
load"./finiteFieldExperiments.m2"


runExperimentOnce(e,numberOfPointsToCollect=>20)
keys runExperimentOnce     

runExperiment(e,2000,numberOfPointsToCollect=>20)
apply(keys e.points,i->#(e.points#i))
