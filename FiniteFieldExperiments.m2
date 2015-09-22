-- finite field experiments

newPackage(
     "FiniteFieldExperiments",
     Version => "0.1", 
     Date => "25.09.2014",
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




-- Map eigentlich Morphism
export {
  "estimateDecomposition", 
  "estimateStratification",
  "estimateCodim", 
  "estimateNumberOfComponents",              
  "createInterpolatedIdeal",
  "createIterator",
  "createRandomPointIterator",
  "Map",
  "interpolateBB",
  "interpolate",
  "isOnComponent",
  "poissonEstimate",
  "Experiment",
  "RandomExperiment",
  "InterpolationImage",
  "InterpolatedImage",
  "FFELogger",
  "ringCardinality",
  "createInterpolatedImage",
  "netEstimatedDecomposition"
}

FiniteFieldExperimentsProtect = ()->
(
  protect bareIdeals;
  protect experiment;
  protect reset;
  protect setPointIterator;
  protect setPointGenerator;
  protect printInterpolatedIdeals;
  protect clearWatchedProperties;
  protect testDebug;
  protect next;
  protect begin;
  protect count;
  protect countinfo;
  protect point;
  protect collectedCount;
  protect pointKeys; 
  protect points;
  protect trials;
  protect createAllInterpolatedIdeals;
  protect interpolatedIdealKeys;

  protect coefficientRingCardinality;
  protect pointLists;
  protect pointsByKey;
  protect countData;
  protect setIsInteresting;
  protect isInteresting;
  protect interpolatedIdeals;

  protect getExperimentData;
  protect setRecordedProperties;
  protect recordProperty;
  protect ignoreProperty;
  protect ignoreProperties;

  protect update;
  protect updateExperiment;

  protect saveData;
  protect loadData;

  protect propertyList;
  protect clear;
  protect rankJacobianAtKey;
  protect watchProperty;
  protect watchProperties;
  protect recordProperties;

  protect propertyName;
  protect propertyAt;

  protect tryProperty;

  protect recordedProperties;
  protect watchedProperties;

  protect useRankJacobianAt;
  protect usedRankJacobianAt;


  protect minPointsPerComponent;
  protect setMinPointsPerComponent;
  protect stratificationIntervalView;
  protect countsByCount;  --protect countsByCount;

  protect experimentData; 
  protect isRandom;
  protect compatible;
  protect membershipPrecision;
  protect setMembershipPrecision;
  --protect createMapHelper;
)

 

FiniteFieldExperimentsExport  = ()->
(
    exportMutable(bareIdeals);
  exportMutable(experiment);
  exportMutable(reset);
  exportMutable(setPointIterator);
  exportMutable(setPointGenerator);
  exportMutable(printInterpolatedIdeals);
  exportMutable(membershipPrecision);
  exportMutable(setMembershipPrecision);

  exportMutable(clearWatchedProperties);
  exportMutable(testDebug);
  exportMutable(next);
  exportMutable(begin);
  exportMutable(count);
  exportMutable(countinfo);
  exportMutable(point);
  exportMutable(collectedCount);
  exportMutable(pointKeys);
  exportMutable(points);
  exportMutable(trials);
  exportMutable(interpolatedIdeals);
  exportMutable(createAllInterpolatedIdeals);
  exportMutable(interpolatedIdealKeys);
 
  exportMutable(coefficientRingCardinality);
  exportMutable(pointLists);
  exportMutable(pointsByKey);

  exportMutable(countData);

  exportMutable(setIsInteresting);
  exportMutable(isInteresting);

  exportMutable(getExperimentData);
  exportMutable(setRecordedProperties);
  exportMutable(recordProperty);
  exportMutable(ignoreProperty);
  exportMutable(ignoreProperties);


  exportMutable(update);
  exportMutable(updateExperiment);
  exportMutable(saveData);
  exportMutable(loadData);

 exportMutable(propertyList);
 exportMutable(clear);
 exportMutable(rankJacobianAtKey);
 exportMutable(watchProperties);
 exportMutable(watchProperty);
 exportMutable(recordProperties);
 exportMutable(propertyName);
 exportMutable(propertyAt);

 exportMutable(tryProperty);

 exportMutable(recordedProperties);
 exportMutable(watchedProperties);

 exportMutable(useRankJacobianAt);
 exportMutable(usedRankJacobianAt);
 

  exportMutable(minPointsPerComponent);
  exportMutable(setMinPointsPerComponent);
  exportMutable(stratificationIntervalView);

 
  exportMutable(countsByCount);

 exportMutable(estimateStratification2);
 exportMutable(experimentData);
 exportMutable(isRandom);
 exportMutable(compatible);
 
 exportMutable(createExperimentData);
 exportMutable(createMapHelper);
) 

undocumented {
estimateStratification2, -- deprecated
propertyList,              --internal variable
propertyName,              --internal variable
countData,                 --internal variable
createExperimentData,      --internal, only used for IO
createIterator,            --document in random point iterator, later.
createRandomPointIterator, --document in random point iterator, later.
begin,                         --document in random point iterator, later.
next,                          --document in random point iterator, later.
point,                         --document in random point iterator, later.
reset,                         --iterator
compatible,           --internal method
createMapHelper,
estimateNumberOfComponents,    --internal
estimateCodim,      --deprecated
ringCardinality,    --internal
experiment,         --internal
experimentData,     --internal
getExperimentData,
isInteresting,      --internal
interpolatedIdeals, --internal
isRandom,           --internal
rankJacobianAtKey,      --document later, redesign
loadData,           --IO; in development
propertyAt,         --unnecessary/deprecated
saveData,           --IO; in development
savedExperimentData412398472398473923847, --IO; in development
testDebug,
update,              --intern
updateExperiment,    -- newFeature not ready.
FFELogger,            -- internal for debug.
recordProperties,     -- replace with watchProperties
recordProperty,       -- replace with watchProperty
recordedProperties,   -- replace with watchedProperties
setRecordedProperties  -- replace with watchProperties
}


needsPackage "SimpleDoc";
needsPackage "Text";


exportMutable(savedExperimentData412398472398473923847);

FiniteFieldExperimentsExport();

if FiniteFieldExperiments#Options#DebuggingMode then
  errorDepth=0;


FFELogger = Logger("FiniteFieldExperiments");
FFELogger.setLogLevel(2);

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


 




Experiment = new Type of HashTable;

load "./FiniteFieldExperiments/Interpolation.m2";

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
    
    localIsCertainlySingularAt := null;
    p.isCertainlySingularAt = ()->
    (
        if (localIsCertainlySingularAt=!=null) then return localIsCertainlySingularAt;
        localIsCertainlySingularAt = blackBox.isCertainlySingularAt(point);
        return localIsCertainlySingularAt;
    );
);

--new ExperimentData from HashTable := (E,coeffRing) -> (
new ExperimentData from Ring := (E,coeffRing) -> (
     ffelog.debug (toString E);
     e := new MutableHashTable;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable;
    -- format: key=>{ideal, maxDegree, name} --later: data type for interpolated ideal
     e.countData = new Tally;
     e.trials = 0;
     e.propertyList = {};
     e.isRandom = null;
     return e;
);






estimateNumberOfComponents(Experiment,List) := HashTable => opts->    (experiment,key) -> 
(
      countData := experiment.count();
      posRankJacobianAt := experiment.position( "rankJacobianAt" );
     if posRankJacobianAt === null then error("To estimate number of components, \"rankJacobianAt\" must be watched");
    
     cardinality := experiment.coefficientRingCardinality();
     estimateNumberOfComponents(
      experiment.trials(),
      key#posRankJacobianAt,
      countData#key,
      cardinality,opts)
)

createExperimentData = (coeffRing,points,countData, trials, propertyList,isRandom) -> (
     e := new ExperimentData;
     e.coefficientRing = coeffRing;
     e.points = new MutableHashTable from points;
     e.countData = countData;
     e.trials = trials;
     e.propertyList = propertyList;
     e.isRandom = isRandom;
     return e;
)

new ExperimentData from ExperimentData := (E,ed) -> (
     -- print (toString E);
     e := new MutableHashTable;
     e.coefficientRing = ed.coefficientRing;
     e.points = copy ed.points;
     e.countData = copy ed.countData;
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
            keys ed1.countData      == keys ed2.countData and  
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
        if not ed1.countData#key == ed2.countData#key then 
          return false;
    );
     return true;
  )
  else return false; 
);

--coeffRing,points,countData,trials,propertyList,isRandom

toExternalString(ExperimentData) := String=> (ed)->
( 
   return "createExperimentData " | "(" 
   | toExternalString ed.coefficientRing | ", \n "
   | toString (new HashTable from ed.points) | ", \n"
   | toExternalString ed.countData | ",\n"
   | toExternalString ed.trials | ",\n"
   | toExternalString ed.propertyList | ",\n"
   | toExternalString ed.isRandom 
   | ")" ;
)

ExperimentData + ExperimentData := ExperimentData=>(ed1,ed2)->
(

   if ( ed1.coefficientRing === ed2.coefficientRing and
        ed1.propertyList    == ed2.propertyList and
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
     edNew.countData = ed1.countData + ed2.countData;
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

TEST ///

 rng = QQ[x];
 bbI = new BlackBoxIdeal from ideal(x)
 e = new Experiment from bbI
 
  weakPoint =()->
  (
     num := random(ZZ);
     if odd num then 
     return random(QQ^1,QQ^1)
     else
     return null;
  );
  wrpi = createRandomPointIterator(weakPoint);
  e.setPointIterator(wrpi);
  e.run(100)
  e.trials();
  assert(e.trials()>100);
  wrpi.reset();
  assert(wrpi.position()==0);
  assert(wrpi.point()===null);
  
///

createRandomPointIterator (Ring,ZZ) := HashTable =>( coeffRing,numVariables )->
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
            ri := createRandomPointIterator(coeffRing, numVariables );
            ri.next();
            return ri;
        );

        rpi.point=()-> randomPoint;

       return new HashTable from rpi;
    );

TEST ///
  debug FiniteFieldExperiments
  FiniteFieldExperimentsProtect()
  pointIterator = createRandomPointIterator(ZZ/7, 5)
  pointIterator.next()
  pointIterator.point()
  apply(99, i-> pointIterator.next() )
  assert(pointIterator.position()==100);
  pointIterator.reset();
  assert(pointIterator.position()==0);
  assert(pointIterator.point()===null);
///


doc ///
   Key
        (createRandomPointIterator, Ring, ZZ)
   Headline
        create  iterator over random points given by ground ring and number of coordinates.
   Usage   
        pointIterator = createRandomPointIterator(R, dim)
   Inputs  
         rng: Ring
            the coefficient ring 
         n:ZZ
            dimension of the vector space over rng 
   Description
        Text
           Create an iterator over random points \in R^{n}, given as matrices \break
           See also PointIterator.
        Example
           rng = QQ
           numCoordinates := 4_ZZ
           pointIterator = createRandomPointIterator(rng, numCoordinates);
        Text
           now we are able to generate random points by calling next():
        Example
           pointIterator.next()
           pointIterator.point()
           pointIterator.position()
           pointIterator.next()
           pointIterator.point()
           pointIterator.position()       
///

doc ///
   Key
        (createRandomPointIterator, Function)
   Headline
        create iterator over points given by a point generator.
   Usage   
        pointIterator = createRandomPointIterator(weakPointGenerator)
   Inputs  
        weakPointGenerator: Function
            a function which generates a random point of type Matrix or returns null
   Description
        Text
           Create an iterator over random points using a given random point generator \break
           See also PointIterator.
        Example
           weakPointGenerator := ()-> (if odd random(ZZ) then  random(QQ^1,QQ^3) );
           pointIterator = createRandomPointIterator(weakPointGenerator);
        Text
           now we are able to generate random points by calling next():
        Example
           pointIterator.next()
           pointIterator.point()
           pointIterator.position()
           pointIterator.next()
           pointIterator.point()
           pointIterator.position()        
///

createIterator = method();
createIterator (List) := HashTable =>( pPoints )->
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
            localPointIterator := createIterator( pPoints);
            localPointIterator.next();
            return localPointIterator;
        );

        pIterator.point=()-> point;

       return new HashTable from pIterator;
    );

doc ///
   Key
        (createIterator)
   Headline
        create iterator over elements given by a list
   Usage   
        pointIterator = createIterator(list)
   Inputs  
        list: List
            a list of ite
   Description
        Text
           Create an iterator over random points using a given random point generator \break
           See also PointIterator.
        Example
           pointList = apply (4, i -> random( ZZ^1, ZZ^3 ) )
           pointIterator = createIterator(pointList);
        Text
           now we are able to iterate over all points:
        Example
           while  pointIterator.next() do pointIterator.point()
           pointIterator.next()
///


EstimatedDecomposition = new Type of HashTable;


createEstimatedDecomposition := (bbi, estimate, watchedProperties)->
(
    ht :=   new HashTable from {  
        "blackBox" =>bbi,
        "estimate"=> estimate,
        "watchedProperties" =>watchedProperties,
    
    };
    return new EstimatedDecomposition from ht;
)

netEstimatedDecomposition = (ed)->
(   
   formatinfo :=  {
                   -- net "-- format: ",
                      net ("--(estimated codim, estimated number of components [confidence interval]) <= "
                        | toString ed#"watchedProperties" | " )") 
                  };
   estimate := ed#"estimate";

   sss := sort apply (estimate, entry ->  net( entry#0)  | " <= " | net (entry#1 ) );

   sss = formatinfo |sss;
   return stack sss;
)

net (EstimatedDecomposition) := Net =>(estimatedDecomposition)->
(
    return (net netEstimatedDecomposition(estimatedDecomposition));
)



estimateDecompositionOld := (experiment) -> (
       countData := experiment.count();
       posRankJacobianAt := experiment.position( "rankJacobianAt" );
       if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");

       cardinality := experiment.coefficientRingCardinality();

       print( "(estimated codim, estimated number of components [confidence interval] <= {watched Properties})");
       print "";
       apply(sort apply(keys countData,key->
         (net(
              key#posRankJacobianAt,
              estimateNumberOfComponents(
               experiment.trials(),
               key#posRankJacobianAt,
               countData#key,
               cardinality ) )
         ) |" <= " |net key ),
       print);
   );
 

estimateDecomposition =  (experiment) -> (
       posRankJacobianAt := experiment.position( "rankJacobianAt" );
       if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");

       estimate := flatten apply(keys experiment.count(), key-> ( (key#posRankJacobianAt, estimateNumberOfComponents(experiment,key)), key) );
        
       return createEstimatedDecomposition(experiment.blackBoxIdeal(), estimate, experiment.watchedProperties() );
)

--needs to be documented
 
EstimatedStratification = new Type of HashTable;

createEstimatedStratification := (bbi, estimate, watchedProperties)->
(
    ht :=   new HashTable from {  
        "blackBox" =>bbi,
        "estimate"=> estimate,
        "watchedProperties" =>watchedProperties,
    
    };
    return new EstimatedStratification from ht;
)

netEstimatedStratification = (es)->
(
  
   formatinfo :=  {
                   -- net "-- format: ",
                      net ("-- estimated codim <= "
                        | toString es#"watchedProperties" | " )") 
                  };

   estimate := es#"estimate";

   sss := sort apply (estimate, entry ->  net( entry#0)  | " <= " | net (entry#1) );

   sss = formatinfo |sss;
   return stack sss;
)

net (EstimatedStratification) := Net =>(es)->
(
    return (net netEstimatedStratification(es));
)

possibleCodimComponents = (found,trials,orderK) -> (
     estimate := 1/trials*poissonEstimate(found);
     possibleAnswers := flatten apply(10,c->apply(4,d->(c,d+1,c-log(d+1)/log(orderK))));
     apply(
          select(possibleAnswers,a->((-log(estimate.min)/log(orderK)>a#2) and (-log(estimate.max)/log(orderK)<a#2))),
          a->(a#0,a#1)
          )
     )

-- round a real number to one decimal place
-- and add a zero before and after the decimal point
-- if necessary (works only for positive numbers)
round1withZeros = (rr) -> (
     if rr<0 then error;
     if rr==0 then return "0.0";
     if rr<1 then return "0"|toString(round(1,rr));
     if round(1,rr) == round(0,rr) then return toString(round(1,rr))|".0";
     return toString(round(1,rr))
     )

roundCodim = (trials,found,orderK) -> (
     estimate := 1/trials*poissonEstimate(found);
     if (log(estimate.max)/log(orderK) - log(estimate.min)/log(orderK))<0.8 
     then return round1withZeros((log(trials)-log(found))/log(orderK))
     else return "..."
     )
 
estimateStratification =  (experiment) -> (
     trials := experiment.trials();
     orderK := experiment.coefficientRingCardinality(); -- this must be read from the experimentdata
     -- (jk): need more advice. Did we want to use a different ring for search that the ideal coefficient ring? If so, 

     countData := experiment.count();

     -- sort keys by number of occurence
     sortedKeysByOccurence := apply(reverse sort apply(keys countData,k->(countData#k,k)), i->i#1);

     --apply(sortedKeysByOccurence,k->(
     --      --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
     --      print (net(round(1,(log(trials)-log(countData#k))/log(orderK)))|" <= "|net(k)))
     --      );

     --estimate := flatten apply(sortedKeysByOccurence,k->( (round(1,(log(trials)-log(countData#k))/log(orderK))), (k)) );
     estimate := flatten apply(sortedKeysByOccurence,k->( (round1withZeros((log(trials)-log(countData#k))/log(orderK))), (k)) );
     --estimate := flatten apply(sortedKeysByOccurence,k->( (roundCodim(trials,countData#k,orderK))), (k)) );
     return  createEstimatedStratification(experiment.blackBoxIdeal(), estimate, experiment.watchedProperties() );
)

-- deprecated
estimateStratification2 = (e) -> (
     --count := e.countData();
     trials := e.trials();
     orderK := (e.coefficientRing()).order; -- this must be read from the experimentdata
     -- print "--";
     -- apply(e.countsByCount(),i->(
     --       --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
     --       print (net(round(1,(log(trials)-log(i#0))/log(orderK)))|" <= "|net(i#1)))
     --       )
     --  ;
     estimate := flatten apply(e.countsByCount(),k->( round1withZeros((log(trials)-log(k#0))/log(orderK)), (k#1)) );
     return  createEstimatedStratification(e.blackBoxIdeal(), estimate, e.watchedProperties() );
)



   stratificationIntervalView := (stratificationData )->
   (
     
      prerearrangedData := new HashTable from apply( keys stratificationData, key -> (stratificationData#key=>(stratificationData#key,key )));
      toSort := apply( keys stratificationData, key -> (stratificationData#key));
      sorted := sort (toSort); 
      rearrangedData := apply (sorted, key-> (prerearrangedData#key) );
      return new HashTable from rearrangedData;
   );



   --#??? fuer den Fall dass dich die Schlüssel nicht sortieren lassen.

   -- countsByCount() sorts the statistic by Number
   countsByCount := (experimentData)->
   (
      counts := experimentData.countData;
      prerearrangedData := new MutableHashTable;
     
      -- 1. use the count number as key and the values of corresponding watched properties as values.
      for key in keys counts do
      (
          if not  prerearrangedData#?(counts#key) then 
             prerearrangedData#(counts#key)= { key }
          else
          (
              prerearrangedData#(counts#key)=  prerearrangedData#(counts#key) | { key };
          );

      );
      -- 2. now create a list of the point count and sort it
      ---- bug fixed (unique missing) todo: test for this bug!
      toSort := unique apply( keys counts, key -> (counts#key));
      sorted := sort (toSort); 
      rearrangedData := {};

      -- 3.  asseble the result
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
      ffelog.warning("-- warning : you got a reference to experiment data, do not modify! ");
      return experimentData;
   );

   -- some of the  following values initialized at end ( e.g. propertiesAt initialization depends presence of some functions defined later)
   minPointsPerComponent := 10;

   rankJacobianAtKey := null;
   rankJacobianAt := null;

   propertiesAt := (point)->{};
   isInteresting := (point)->true;

   pointIterator := createRandomPointIterator ( experimentData.coefficientRing, blackBoxIdeal.numVariables );

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
               experiment.rankJacobianAtKey()         ==  re2.rankJacobianAtKey() and
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
   

   experiment.useRankJacobianAt = (rankJacobianAtName)->
   ( 
       if experiment.trials()=!=0 then error ("cannot change rankJacobianAt  - experiment was already run! You could clear() the statistics and retry. ");

       if rankJacobianAtName===null then
       (
         rankJacobianAtKey = null;
         rankJacobianAt = null ;
         return;
       );
       if blackBoxIdeal.hasPointProperty(rankJacobianAtName) then 
       (
          rankJacobianAtKey = rankJacobianAtName;
          rankJacobianAt = blackBoxIdeal.pointProperty(rankJacobianAtName) ;
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
        experimentData.countData = new Tally;
   );
 
   experiment.rankJacobianAtKey = ()->
   (
       return rankJacobianAtKey;
   );

   experiment.usedRankJacobianAt = ()->
   (
       return rankJacobianAtKey;
   );


    runExperimentOnce(ExperimentData, Matrix,ZZ ) := Tally => (experimentData, point, wantedPoints) -> 
    (  
        K := experimentData.coefficientRing;
        --prime = char K
        numVariables := blackBoxIdeal.numVariables;

        -- if ideal vanishes on the random point do something
        if experiment.isInteresting(point) then   
        (
            countKey := propertiesAt(point); -- todo: countKey: better naming?
 
            -- countData number of found points for each rank and property
            experimentData.countData = experimentData.countData + tally {countKey};

          
            if  rankJacobianAt =!= null then 
            (
                FFELogger.debug( "update wanted points" );
                rankJacobian := rankJacobianAt(point);  
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

    experiment.setPointGenerator  = (pGen) ->
    (
        wrip := createRandomPointIterator(pGen);
        experiment.setPointIterator(wrip);
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
        pointIterator := createIterator (  experiment.points() );
        experimentData.points = new MutableHashTable;
        experimentData.countData = new Tally;

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

  experiment.clearWatchedProperties = (   )->
  (
      experiment.clear(); 
      setRecordedPropertiesInternal({});
  );


   experiment.setRecordedProperties = ( propertyStringList )->
   (  
      if experiment.trials()=!=0 then error ("cannot change watched properties - experiment was already run! Clear statistics and retry.");

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
       return experiment.count();
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
   experiment.count = ()->
   (
      return new Tally from experimentData.countData;
   );

   experiment.countinfo = ()->
   (
     strcountinfo := "-- count() structure: \n-- values of " |toString experimentData.propertyList  | " =>  number of random points ";
     return strcountinfo;
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

   -- todo: test if changing blackBoxIdeal.rankJacobianAt is transparent (means after blackBoxIdeal.updatePointProperty("rankJacobianAt") the new one is called)
 
    experiment.watchProperties( {} );

       if  blackBoxIdeal.hasPointProperty("isZeroAt") then
       (
          --print("isInteresting"| toString  (blackBoxIdeal.isZeroAt ));
          isInteresting = blackBoxIdeal.pointProperty("isZeroAt");

          --     experiment.watchProperties( {"isZeroAt"} );
       );

       if blackBoxIdeal.hasPointProperty("jacobianAt") then 
       (
          assert ( blackBoxIdeal.hasPointProperty("rankJacobianAt") );
       );

       if blackBoxIdeal.hasPointProperty("rankJacobianAt") then 
       (
          assert ( blackBoxIdeal.hasPointProperty("rankJacobianAt") );
          rankJacobianAtKey = "rankJacobianAt";
          rankJacobianAt =  blackBoxIdeal.pointProperty("rankJacobianAt");


            experiment.watchProperties( {"rankJacobianAt"} );

       );



   --end init:


  -- to fix the issue that the internal experiment reference
  -- is not of type Experiment, there are two solutions:
  --  either to introduce to different method signatures (ones accept a HashTable and others an Experiment)
  -- or to make the internal experiment variable as an Experiment , too  - done with 'newClass'

    experiment.tryProperty = (tProp) -> (
       ffelog.info("-- ( " | toString experiment.watchedProperties() |" | " | tProp | " ) => count " );
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
    
 
 

   experiment = newClass( Experiment, experiment );
   ffelog.debug ("type of internal experiment variable is : " | toString class experiment );
   return experiment; 
   -- return new HashTable from experiment; 
);

TEST ///

 rng = QQ[x];
 bbI = new BlackBoxIdeal from ideal(x)
 e = new Experiment from bbI
 
  weakPoint =()->
  (
     num := random(ZZ);
     if odd num then 
     return random(QQ^1,QQ^1)
     else
     return null;
  );
  wrpi = createRandomPointIterator(weakPoint);
  e.setPointGenerator(weakPoint);

  
  
///

 


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
        an unified interface to an experiment
   Description
         Text
            With an {\tt  Experiment } it is possible to check point properties of an @TO BlackBoxParameterSpace@ or @TO BlackBoxIdeal@ 
            at random points and collect user-defined statistics. \break
            If the black box from supports evaluation, then at each point the jacobian can be computed
            and jets at smooth ones. From the collected statistics a heuristic decomposition can be estimated and finally performed using interpolation methods, see @TO "Experiment example"@
            
            The {\tt  Experiment } objects implements the following interface    \break 

            construction :\break
            \,\, \bullet \, @TO "(NewFromMethod, Experiment, BlackBoxParameterSpace)" @ \break
            \break \break

            properties \break
            \,\, \bullet \,{\tt coefficientRing}:  . \break
            \,\, \bullet \,{\tt coefficientRingCardinality}:  . \break
            \,\, \bullet \,{\tt collectedCount}: . \break
            \,\, \bullet \,{\tt count}:  \break
            \,\, \bullet \,{\tt countsByCount}: \break
            \,\, \bullet \,{\tt membershipPrecision}:  \break
            \,\, \bullet \,{\tt minPointsPerComponent}: \break
            \,\, \bullet \,{\tt points}:  \break
            \,\, \bullet \,{\tt pointKeys}: \break
            \,\, \bullet \,{\tt pointLists}:  \break
            \,\, \bullet \,{\tt pointsByKey}:  \break
            \,\, \bullet \,{\tt trials}:  \break
            \,\, \bullet \,{\tt watchedProperties}: \break
            \,\, \bullet \,{\tt stratificationIntervalView}:  \break
            \break \break

            methods: \break
            \,\, \bullet \,{\tt setIsInteresting}: set a filter for points to consider. \break
            \,\, \bullet \,{\tt setMembershipPrecision}:  \break
            \,\, \bullet \,{\tt setMinPointsPerComponent}:  \break
            \,\, \bullet \,{\tt setPointGenerator}: \break
            \,\, \bullet \,{\tt setPointIterator}:  \break
            \,\, \bullet \,{\tt useJacobianAt}:  \break
            \,\, \bullet \,{\tt watchProperties, watchProperty }:  \break
            \,\, \bullet \,{\tt ignoreProperties,  ignoreProperty }:  \break
            \,\, \bullet \,{\tt clearWatchedProperties}:  \break
            \,\, \bullet \,{\tt run}: find interesting points \break
            \,\, \bullet \,{\tt tryProperty}:  \break
            \,\, \bullet \,{\tt clear}: \break
           
        Text
            \break  For an example see @TO "Experiment example"@
///


doc ///
   Key
        InterpolatedImage
   Headline
        blabla
   Description
        Text
                interpolatedIdealKeys 
                InterpolationImage 
                isOnComponent 
                printInterpolatedIdeals 
///

doc ///
   Key
        "new Experiment"
        (NewFromMethod, Experiment, BlackBoxParameterSpace)
        (NewFromMethod, Experiment, Thing)
   Headline
        type to collect data from a finite field experiment
   Usage   
        e = new Experiment from bb
        e.run(trials)
        e.count()
        e.pointLists()
        e.estimateStratification()
        e.estimateDecomposition()
   Inputs  
        bb:HashTable 
            a blackBoxIdeal
   Description
        Text
            Creates an @TO Experiment@ from a black box, see @TO BlackBoxIdeals@.
                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;       
        Text
           \break Now we are able to create the Experiment object
       
        Example
           e = new Experiment from bb;
   Caveat
        does not check if the ideal ring is a quotient ring (not supported)
   SeeAlso
        "Experiment example"
        "run Experiment"
        estimateDecomposition
        createAllInterpolatedIdeals        
       
///


---    Outputs
---        :HashTable
---            observed point number for each value tuple of the watched properties

doc ///
   Key
        "run Experiment"
   Headline
        run an experiment
   Usage   
        e.run(trials)
   Inputs  
        e:Experiment 
            a finite field Experiment
        trials:ZZ
            number of evaluations at random points
   Description
        Text
           To run an experiment we first create an BlackBoxIdeal or a BlackBoxParameterSpage  we want to analyse
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;       
        Text
           Now we create the experiment object
        Example
           e = new Experiment from bb;
        Text
           If a black box has a property "rankJacobianAt" it is
           automatically watched:
        Example
           e.watchedProperties() 
        Text
           Now we run the experiment by evaluating at 1250 random points 
        Example
           time e.run(1250)    
        Text
           As return value we get the collected point number grouped by the values of the watched properties.
           Later these statistics are acceccible via count() 
        Example
           e.count()            
           e.countinfo()
        Text
            dfdf
   SeeAlso
        watchProperty
        watchProperties
        watchedProperties


///

doc ///
    Key
        FiniteFieldExperiments
    Headline
          heuristic decomposition of black box ideals
    Description
        Text
            With an {\tt  Experiment } it is possible to check point properties of an @TO BlackBoxParameterSpace@ or @TO BlackBoxIdeal@ 
            at random points and collect them. \break
            From the collected statistics a heuristic decomposition can be estimated in case the black box supports the computation of the 
            jacobian rank at a point.
            \break
            Finally a herustic decomposition can be computed using interpolation methods if the black box supports jet calculations.
            \break \break
            See @TO "Experiment example"@ for a tutorial.
                        
      
    Caveat
         The package development is at alpha status and the package is not threadsafe. 
         The interpolation is not time-optimized.
         The documentation is not finished.
///


doc ///
   Key
        "Experiment example"
   Headline
        typical finite field experiment example
   Usage   
        e = new Experiment from bb
        e.run(trials)
        e.count()
        e.pointLists()
        e.estimateStratification()
        e.estimateDecomposition()
   Inputs  
        bb:HashTable 
            a blackBoxIdeal
   Description
        Text
            With an {\tt  Experiment } it is possible to check point properties of an @TO BlackBoxParameterSpace@ or @TO BlackBoxIdeal@ 
            at random points and collect user-defined statistics.
            If the black box supports evaluation, then at each point the jacobian can be computed
            and jets at smooth ones. 
            From the collected statistics a heuristic decomposition can be estimated and finally performed using interpolation methods.

            Here is a typical extremely simple minded application:
                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;       
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
        Text
           Our black box constructed several default point properties, including @TO isZeroAt@.
           The user may add new properties, see 
           running the experiment will trigger evaluation the ideal at random points.
        Example
           bb.hasPointProperty("isZeroAt")
           bb.isZeroAt(matrix{{0_K,0,1}})
        Text 
           If the point is in the vanishing set of the ideal (i.e. either on the point
           or on the line), the point will be considered as insteresting, see @TO setIsInteresting@.
           For all interesting points the experiment will compute the @TO watchedProperties@
           and keep that statistics.
        Text
           \break If a black box has a property "rankJacobianAt" it is
           automatically watched.
        Example
           e.watchedProperties() 
        Text 
           The experiment calculate the rank of the jacobi matrix at interesting points.
           (2 on the line, 1 on the plane, 0 in the origin). 
        Example
           bb.isZeroAt(matrix{{0_K,0,1}})
           bb.rankJacobianAt(matrix{{0_K,0,1}})
           bb.rankJacobianAt(matrix{{1_K,2,0}})
           bb.rankJacobianAt(matrix{{0_K,0,0}}) 
        Text
           \break Now we run the experiment by evaluating at 1250 random points:  
        Example
           time e.run(1250)      
        Text
           As return value we get the collected point number grouped by the values of the watched properties.
           Later these statistics are acceccible via count() :
        Example
           e.count()            
           e.countinfo()
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
           Points with a particular set of properties can be selected
           like this:
        Example
           e.pointsByKey({2})
        Text
           Since one always finds many points found on components of low
           codimension it is not useful to remember all of them. The experiment
           remembers by default only about 10 points per component:
        Example
           e.minPointsPerComponent()
           e.collectedCount() 
        Text
           Here we have not collected exactly 10 points per component since the experiment uses the upper end of the confidence interval for the number of components ( see @TO estimateNumberOfComponents@) as guide for the number of points to keep.
           The amount of stored points can be adjusted:
        Example
           e.setMinPointsPerComponent(20)
           -- collect about 20 points per component now:
           time e.run(1250);
           e.collectedCount() 
        Text          
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
           A doubling of the number of experiments is expected to divide the
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
    e.useRankJacobianAt("rankJacobianAt");
    e.useRankJacobianAt(null);
   
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
    e.rankJacobianAtKey()

    bbRankM.knownPointProperties()

    -- e.stratificationIntervalView() -- test fails here; d

///

doc ///
   Key
        "clear"
   Headline
        deletes data of an experiment to start from a clean slate
   Usage   
        e.clear()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            Sometimes one wants to restart an experiment from a clean slate,
            for example if one has changed some properties or has implemented
            new properties one wants to watch. 
            
            This can be done with clear. 
            
            More precisely clear erases the collected points, the statistics.
            The number of trials is set to zero. If a iterator is used to 
            generate the points, this is also reset to the starting point.
            This is useful if the iterator is 
            used to enumerate a given set of non random points (e.g all points
            in IP^n). 
            
            The black box itself and the list of watched properties
            are not changed.
            
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)          
           e.trials()
           e.count()
           e.pointLists()
           e.watchedProperties()
        Text
           \break Now we clear the statistics and point lists:
        Example
           e.clear()
           e.trials()
           e.count()
           e.pointLists()
           e.watchedProperties()
///

doc ///
   Key
        "clearWatchedProperties"
   Headline
        deletes the list of watched properties
   Usage   
        e.clearWatchedProperties()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            This is the same as  @TO2{clear,"e.clear()"}@ with the additional effect that
            the watched property list is erased.
          
            The black box itself is not changed.
            
            It is not usesful to erase only the watched property list
            without setting the statistics to zero, since the statistics
            count the number of times a particular property has occured.
            
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)          
           e.trials()
           e.count()
           e.pointLists()
           e.watchedProperties()
        Text
           \break Now we clear the statistics, the point lists and the watched properties.
        Example
           e.clearWatchedProperties()
           e.trials()
           e.count()
           e.pointLists()
           e.watchedProperties()
///

doc ///
   Key
        "collectedCount"
   Headline
        counts the number of collected points
   Usage   
        e.collectedCount()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            On strata of low codimension many points are found. To
            avoid memory problems only a small number of points
            are collected (see setMinPointsPerComponent). Therefore
            the number of collected points is usually smaller than
            the number of found points. 
            
            collectedCount() returns an HashTable containing
            the number of collected points
            for each combination of properties.

            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)          
           e.collectedCount()
           e.pointLists()
           e.pointsByKey({2})
           e.minPointsPerComponent()
        Text
           \break Notice that the number of collected points can be larger than
           the number minPointsPerComponent() since the experiment
           tries to estimate the number of components for each combination
           of properties. In the beginning where only a few points have
           been found the statistics might be so errorprone that some extra
           points are collected. 
   SeeAlso
          pointLists
          pointsByKey
          minPointsPerComponent
          setMinPointsPerComponent             
///

doc ///
   Key
        "count"
   Headline
        shows how often points with particular properties were found by an experiment
   Usage   
        e.count()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            Returns a Tally with the number of points that
            an experiment has found so far for each combination
            of wached properties.

            This gives a rough overview over which properties 
            occur in large/small strata. For an estimation of the
            codimension and number of components for each such stratum
            (see estimateDecomposition).
            
            e.count() is called automatically when e.run is finished.
            
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(125) 
           e.count()
        Text
           \break There are 125 point over F_5 in A^3. Of these 25 lie on 
           the plane and 5 on the line. As one can see this is also about
           the number found by a random search.  
        Example
           e.estimateDecomposition()
        Text
           The expected fraction of points on a codim c component is 
           1/p^c. This is used by estimateDecomposition to guess the 
           number of components in each codimension. We see that up to 
           a margin of error this gives the correct answer of 1 component
           in codimension 1 and 2 each.       
///

doc ///
   Key
        "estimateDecomposition"
   Headline
        shows how often points with particular properties were found by an experiment
   Usage   
        e.estimateDecomposition()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
           Returns a HashTable with an estimate of the number of components 
           in each codimension for all combinations of properties found. 
 
           This is done by assuming that for each component of codimension 
           c the fraction of points found is the expected 1/p^c and 
           by further assuming that the components are smooth at most 
           of there rational points P. In this case
           c is equal to the rank of the jacobi matrix at this point.
                        
           Lets see how this works in an example.                
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(500) 
           e.estimateDecomposition()
        Text
           We see that up to 
           a margin of error this gives the correct answer of 1 component
           in codimension 1 and 2 each and no component of codimension 0.
           
           (the point where the line and the plane intersect has
           rankJacobiMatrix equal to 0. Fortunately there are not enough
           such points to fool the algorithm.)
           
           estimateDecomposition works only if the black box used in this experiment
           provides a property rankJacobianAt and this property is
           watched by the experiment. As explained above,
           this function uses the rank of 
           the jacobi matrix at each point as an estimate for the codimension
           of the component the point lies on. In effect assuming that
           the point lies in the smooth part of the component.      
///

doc ///
   Key
        "countsByCount"
   Headline
        sorts the statistics by number of points found
   Usage   
        e.countsByCount()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            Returns a Tally with the number of points that
            an experiment has found so far for each combination
            of wached properties sorted by the number of points found.
            (not by alpha numeric sort of the properties)

            This gives a rough overview over which properties 
            occur in large/small strata. 
            
            This is particularily usefull when working with an
            black box parameter space that provides no propery
            rankJacobianAt and therefore no direct estimation of
            the codimension of each component.
                   
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(500) 
           e.countsByCount()
///


TEST ///
-- test issue #57
prime = 11
K = ZZ/prime
R = K[x,y,z,w]
betti res (I = minors(2,random(R^{2:0},R^{3:-1})))
projI = eliminate(w,I)

bb = blackBoxIdeal I;
e = new Experiment from bb;
f = new Experiment from bb;
///



doc ///
   Key
        "ignoreProperty"
   Headline
        deletes a property from the list of watched properties
   Usage   
        e.ignoreProperty(propertyName)
   Inputs  
        e:Experiment 
            an Experiment
        propertyName: String
            the name of a property
   Description
        Text
           This removes a property from the list of watched properties.
           This works only if the experiment has been reset with e.clear()
           since otherwise the statistics would be inconsistent.
                                 
           Lets see how this works in an example.                
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.watchedProperties()
        Text
           \break Lets not only watch the codimension of the tangentspace
           at a random point, but also whether the point is probably smooth
        Example
           e.watchProperty("isProbablySmoothAt")
           e.watchedProperties()
           e.run(250)
        Text
           \break Lets assume that from the last experiment we conclude,
           that the smoothness of a point does not yield any interesting 
           information for us.
           In this case we would cease to watch this property in 
           follow up experiments. Before we can do that, we have to
           clear the statistics.
        Example
           e.clear()
           e.ignoreProperty("isProbablySmoothAt")
           e.watchedProperties()
           e.run(250)  
///

doc ///
   Key
        "ignoreProperties"
   Headline
        deletes a property from the list of watched properties
   Usage   
        e.ignoreProperty(L)
   Inputs  
        e:Experiment 
            an Experiment
        L:List
            a list of property names. 
   Description
        Text
           This removes several property from the list of watched properties.
           This works only if the experiment has been reset with e.clear()
           since otherwise the statistics would be inconsistent.
                                 
           Lets see how this works in an artificial but instructive example.                
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.watchedProperties()
        Text
           \break Lets not only watch the codimension of the tangentspace
           at a random point, but also whether the point is probably smooth
           and what the values at the random points are.
        Example
           e.watchProperty("isProbablySmoothAt")
           e.watchProperty("valuesAt")
           e.watchedProperties()
           e.run(250)
        Text
           \break Lets assume that from the last experiment we conclude,
           that the smoothness of a point and the value of
           the polynomial at the point does not yield any interesting 
           information for us.
           In this case we would cease to watch these properties in 
           follow up experiments. Before we can do that, we have to
           clear the statistics.
        Example
           e.clear()
           e.ignoreProperties({"isProbablySmoothAt","valuesAt"})
           e.watchedProperties()
           e.run(250)  
///

doc ///
   Key
        "minPointsPerComponent"
   Headline
        the number of points an experiment tries to collect on each component.
   Usage   
        e.minPointsPerComponent()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            On strata of low codimension many points are found. To
            avoid memory problems only a small number of points
            are collected.
            
            This function returns the number of points that an
            experiment tries to collect on each component. Since
            the number of components is estimated heuristically, the 
            number of points collected is often different from the
            goal given by this function. 
 
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)   
        Text
           \break The number of points the experiment tries to collect per component:
        Example
           e.minPointsPerComponent()
        Text 
           \break The number of point the experiment has collected
        Example      
           e.collectedCount()
           e.run(100)
           e.collectedCount()
        Text
           \break Lets now increase the number of points we want to collect
        Example
           e.setMinPointsPerComponent(20)
           e.minPointsPerComponent()    
           e.run(100)
           e.collectedCount()
           e.run(100)
           e.collectedCount()
   SeeAlso
        setMinPointsPerComponent
        collectedCount
        pointLists
        pointsByKey
///

doc ///
   Key
        "pointLists"
   Headline
        the points an experiment has collected for closer inspection
   Usage   
        e.pointLists()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            On strata of low codimension many points are found. To
            avoid memory problems only a small number of points
            are collected.
            
            This function returns a hashTable containing a list
            of collected points for each combination of properties.
  
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)   
        Text
            \break The number of points the experiment has collected:
        Example      
           e.collectedCount()
        Text
           \break The points themselves:
        Example
           e.pointLists()
        Text
           \break The points for a particular set of properites:
        Example   
           (e.pointLists())#{2}
        Text
           The brackets used in this example are NOT redundant. Better
           readable is the following syntax
        Example
           e.pointsByKey({2})
   SeeAlso
        minPointsPerComponent
        setMinPointsPerComponent
        collectedCount        
        pointsByKey
///        

doc ///
   Key
        "pointsByKey"
   Headline
        the points an experiment has collected for a particular set of properties
   Usage   
        e.pointsByKey(key)
   Inputs  
        e:Experiment 
            an Experiment
        key: List
            of property values
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            On strata of low codimension many points are found. To
            avoid memory problems only a small number of points
            are collected.
            
            This function returns a List 
            of collected points for a given combination of properties.
  
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)   
        Text
            \break The number of points the experiment has collected:
        Example      
           e.collectedCount()
        Text
           \break The points themselves:
        Example
           e.pointLists()
        Text
           \break The points for a particular set of properites:
        Example
           e.pointsByKey({2})
   SeeAlso
        minPointsPerComponent
        setMinPointsPerComponent
        collectedCount        
        pointLists
        pointKeys
///        

doc ///
   Key
        "setMinPointsPerComponent"
   Headline
        change the number of points an experiment tries to collect on each component.
   Usage   
        e.setMinPointsPerComponent(num)
   Inputs  
        e:Experiment 
            an Experiment
        num: ZZ
            the number of points an experiment should try to collect on
            each component.
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            On strata of low codimension many points are found. To
            avoid memory problems only a small number of points
            are collected.
            
            This function changes the number of points that an
            experiment tries to collect on each component. Since
            the number of components is estimated heuristically, the 
            number of points collected is in the end often different from 
            (but close to) the
            goal given by this function. 
 
            Lets see how this works in an example.                
            First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(100)   
        Text
           \break The number of points the experiment tries to collect per component:
        Example
           e.minPointsPerComponent()
        Text 
           \break The number of point the experiment has collected
        Example      
           e.collectedCount()
           e.run(100)
           e.collectedCount()
        Text
           \break Lets now increase the number of points we want to collect
        Example
           e.setMinPointsPerComponent(20)
           e.minPointsPerComponent()    
           e.run(100)
           e.collectedCount()
           e.run(100)
           e.collectedCount()
   SeeAlso
        minPointsPerComponent
        collectedCount
        pointLists
        pointsByKey
///

doc ///
   Key
        "pointKeys"
   Headline
        a list of property combinations found by an experiment
   Usage   
        e.pointKeys()
   Inputs  
        e:Experiment 
            an Experiment
        num: ZZ
            the number of points an experiment should try to collect on
            each component.
   Description
        Text
            An experiment collects a limited number of points
            for each combination of properies it encounters. 
            
            This is useful if one wants to inspect
            points with special properties in more detail. If the 
            points are on a moduli space one can create the
            corresponding object and study it in detail.
            
            
            This function returns a list of those property combinations
            that it has encountered. This is sometimes useful if
            each property is very long or difficult to enter.
            
            Lets see how this works in an example. Here we look
            for matrices with special syzygies
                           
            First we create a black box that makes random matrices of quadrics.
        Example      
           K = ZZ/2;
           R = K[x,y,z,w];
        Text
           We want to study the parameter space of 2x3 matrices
           with quadratic entries. Such matrices are defined by
           60 coefficients.
        Example   
           bb = blackBoxParameterSpace(80,K);
        Text
           We start by making quadrics from 10 coefficiets
        Example
           mons2 = matrix entries transpose super basis(2,R)
           quadricAt = (point10) -> (point10*mons2)_0_0;
           quadricAt(matrix{{1,0,0,0,1,0,0,1,0,1}})
        Text
           Now we make a 2x3 matrix of quadrics from 80 coefficients
        Example
           matrixAt = (point80) -> matrix apply(2,i->
                apply(4,j->(
                          quadricAt(point80_{i*10+j..i*10+j+9})
                          )
                     )
                )          
           matrixAt(random(K^1,K^80))
        Text
           For later use we register this function in the 
           black box parameter space
        Example
           bb = bb.registerPointProperty("matrixAt",matrixAt);
        Text
           Now we look at the syzygies of such a matrix
        Example
           bettiAt = (point80) -> betti res coker matrixAt(point80)
           bettiAt(random(K^1,K^80))
           bb = bb.rpp("bettiAt",bettiAt);
        Text
           Now we 
           make an experiment to study this parameter space
        Example   
           e = new Experiment from bb;
        Text
           We are interested in the betti tableau ot the minimal free resolution.
        Example   
           e.watchProperty("bettiAt")
           e.run(100)
           e.countsByCount()
        Text
           We now want to look at the collected point with special
           betti tableaus
        Example
           e.pointKeys()
           (e.pointKeys())#0
           e.pointsByKey((e.pointKeys())#0)
   SeeAlso
        pointLists
        pointsByKey
        collectedCount
        minPointsPerComponent
        setMinPointsPerComponent
///


doc ///
   Key
        "e.run"
   Headline
        runs an experiment for a given number of trials
   Usage   
        e.run(trials)
   Inputs  
        e:Experiment 
            an Experiment
        trials: ZZ
            the number of random point the experiment shall try out
   Description
        Text
           This is the central function of this package that does all 
           the work. It 
           starts an experiment to evaluate the watched properties
           at a given number of random points. The experiment then
           counts the number of times a particular property is encontered.
        
           If the experiment is made from a black box ideal, only points
           on which the generators of the ideal vanish are counted. If
           the experiment is made form a black box parameter space, all
           points are counted.
        
           Also the experiment stores some points for later inspection. 
           For documentationa about how this works see collectedCount
        
           Lets see how this works in an example.                
        
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.watchedProperties()
           e.run(100)
        Text
           \break The experiment has evaluated the black box ideal at 
           100 random points and for each point where the ideal vanished it calculated
           the rank of the jacobi matrix at this point. Then the
           number of points for each rank are counted.
       
           If one wants better statistics one can run the experiment for more
           trials
        Example
           e.trials()   
           e.run(100)
           e.trials()
        Text
           The experiment has autmatically collected some points
           for each combination of properties
        Example
           e.pointLists()
           e.pointsByKey({2})
           e.collectedCount()
        Text
           Notice that the experiment has not collected all points it
           found. This is done to save memory space.
   SeeAlso
       trials
       watchedProperties
       pointLists
       pointsByKey
       collectedCount
///                 
         
doc ///
   Key
        "trials"
   Headline
        the number of trials an experiment has been run so far
   Usage   
        e.trials()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
           Returns the number of trials an experment has been run so far.
           This might be important to interpret the statistics of 
           the experiment. It is used for example in estimateStratification
           and estimateDecomposition.  
             
           Lets see how this works in an example.                
        
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.watchedProperties()
           e.run(100)
           e.trials()   
           e.run(100)
           e.trials()
///                 
 
doc ///
   Key
        "tryProperty"
   Headline
        evaluate a new property on the collected points.
   Usage   
        e.tryProperty(name)
   Inputs  
        e:Experiment 
            an Experiment
        name:String
            the name of a property.
   Description
        Text
           Evaluates a given property on all collected points.
           
           This is useful if in the course of experimenting one thinks
           of a new property that might help in the analysis of 
           the problem at hand. Before running the complete experiment
           again an watching the new property it is much faster to
           begin by evaluating the new property on the collected points.
           Since the experiment collects points on all interesting strata
           this gives a good first overview of what the new property
           will do.
             
           Lets see how this works in an example.                
        
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break The ideal describes a line and a plane intersecting at the origin. \break     
           \break Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.watchedProperties()
           e.run(250)
           e.collectedCount()
        Text 
           \break Perhaps we are wondering is points with a certain
           rank of the Jacobi matrix are allways singular:
        Example
           e.tryProperty("isCertainlySingularAt")
        Text
           It seems that all points with rank 0 are singular. Indeed
           this is true in this example since rank 0 only occurs at
           the intersection point of the line and plane.
           
           For a more realistic example see (Singularities of cubic surfaces)
   SeeAlso
       watchedProperties
       watchProperty
       watchProperties
///                 
 
end
---


quit -- F11 F11 F12

path = append(path,"/Users/bothmer/Desktop/projekte/strudel/Jakob2010/GitHub/padicLiftM2/")

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"

viewHelp FiniteFieldExperiments
viewHelp BlackBoxIdeals

restart

loadPackage"BlackBoxIdeals"
load "FiniteFieldExperiments.m2"
needsPackage"FiniteFieldExperiments"

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

