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
           Email => "jakobkroeker.academic@spaceship-earth.net", 
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
  "Counts",
  "counts",
  "sortedCounts",
  "richPoints",
  "pointsWithProperties",
  "stats",
  "observedValues",
  "realizedValues",
  "realizedProperties",
  "observedValue",
  "realizedValue",
  "realizedProperty",
  "realizedPointProperties",
  "interpolatedIdeals",
  "bareIdeals",
  "estimateDecomposition", 
  "estimateStratification",
  "estimateCodim", 
  "estimateNumberOfComponents",              
  "createIterator",
  "createRandomPointIterator",
  "Map",
  "isOnComponent",
  "poissonEstimate",
  "Experiment",
  "InterpolationImage",
  "InterpolatedImage",
  "FFELogger",
  "ringCardinality",
  "createInterpolatedImage"
}

FiniteFieldExperimentsProtect = ()->
(
    protect setPosition;
    protect setPoint;
    protect pointData;
    protect setDecompositionConfidenceInterval;
    protect allInterpolatedIdeals;
    protect experiment;
    protect reset;
    protect pointKey;
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
    --protect collectedCount;
    protect pointKeys; 
    --protect points;
    protect trials;
    protect smoothPoints;
    protect createAllInterpolatedIdeals;
    protect interpolatedIdealKeys;

    protect coefficientRingCardinality;
    protect pointLists;
    protect pointsByKey;
    protect countData;
    protect setIsInteresting;
    protect isInteresting;
     

    protect getExperimentData;
    protect ignoreProperty;
    protect propertyIsWatched;
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
    protect setWatchedProperties;

    protect propertyName;
    protect propertyAt;

    protect tryProperty;

   
    --protect watchedProperties;

    protect useRankJacobianAt;
    protect usedRankJacobianAt;


    protect pointsPerComponent;
    protect setPointsPerComponent;
    protect stratificationIntervalView;
    protect countsByCount;  --protect countsByCount;

    protect experimentData; 
    protect isRandom;
    protect compatible;
    protect membershipPrecision;
    protect setMembershipPrecision;
    --protect createMapHelper;
);

 

FiniteFieldExperimentsExport  = ()->
(
    exportMutable("setPosition");
    exportMutable("setPoint");
    exportMutable("pointData");
    exportMutable("allInterpolatedIdeals");
    exportMutable ("setDecompositionConfidenceInterval");
    exportMutable("smoothPoints");

   
    exportMutable("experiment");
    exportMutable("reset");
    exportMutable("pointKey");
    exportMutable("setPointIterator");
    exportMutable("setPointGenerator");
    exportMutable("printInterpolatedIdeals");
    exportMutable("membershipPrecision");
    exportMutable("setMembershipPrecision");

    exportMutable("clearWatchedProperties");
    exportMutable("testDebug");
    exportMutable("next");
    exportMutable("begin");
    exportMutable("count");
    exportMutable("countinfo");
    exportMutable("point");
    exportMutable("collectedCount");
    exportMutable("pointKeys");
    exportMutable("points");
    exportMutable("trials");
 
    exportMutable("createAllInterpolatedIdeals");
    exportMutable("interpolatedIdealKeys");

    exportMutable("coefficientRingCardinality");
    exportMutable("pointLists");
    exportMutable("pointsByKey");

    exportMutable("countData");

    exportMutable("setIsInteresting");
    exportMutable("isInteresting");

    exportMutable("getExperimentData");
    exportMutable("ignoreProperty");
    exportMutable("propertyIsWatched");
    exportMutable("ignoreProperties");


    exportMutable("update");
    exportMutable("updateExperiment");
    exportMutable("saveData");
    exportMutable("loadData");

    exportMutable("propertyList");
    exportMutable("clear");
    exportMutable("rankJacobianAtKey");
    exportMutable("watchProperties");
    exportMutable("watchProperty");
    exportMutable("propertyName");
    exportMutable("propertyAt");

    exportMutable("tryProperty");

    
    exportMutable("watchedProperties");
    exportMutable("setWatchedProperties");

    exportMutable("useRankJacobianAt");
    exportMutable("usedRankJacobianAt");


    exportMutable("pointsPerComponent");
    exportMutable("setPointsPerComponent");
    exportMutable("stratificationIntervalView");


    exportMutable("countsByCount");

    exportMutable("estimateStratification2");
    exportMutable("experimentData");
    exportMutable("isRandom");
    exportMutable("compatible");

    exportMutable("createExperimentData");
);

undocumented {
propertyList,              --internal variable
propertyName,              --internal variable
countData,                 --internal variable
createExperimentData,      --internal, only used for IO
--createIterator,            --document in random point iterator, later.
createRandomPointIterator, --document in random point iterator, later.
begin,                     --document in random point iterator, later.
next,                      --document in random point iterator, later.
point,                     --document in random point iterator, later.
points,                        --a list of all points not sorted into list. 
                               --not used anymore since tryProperty has
                               --been implemented
reset,                         --iterator
compatible,           --internal method

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
saveData,           --IO; in development
savedExperimentData412398472398473923847, --IO; in development
propertyAt,         --unnecessary/deprecated
testDebug,
update,              --intern
updateExperiment,    -- newFeature not ready.
FFELogger,            -- internal for debug.
--watchedProperties,   -- replace with watchedProperties
}


needsPackage "SimpleDoc";
needsPackage "Text";


exportMutable("savedExperimentData412398472398473923847");

FiniteFieldExperimentsExport();

if FiniteFieldExperiments#Options#DebuggingMode then
    errorDepth=0;


FFELogger = Logger("FiniteFieldExperiments");
FFELogger.setLogLevel(2);

ffelog := FFELogger;


--95 % aller Messwerte haben eine Abweichung von höchstens 1,96 * sigma (Normalverteilung)

-- rename 'confidence' to 'confidenceInSigma' ?

poissonEstimate = method( Options => {"confidence" => 1.96} );

poissonEstimate(ZZ) := HashTable => opts -> (numPoints) -> 
( 
    -- we use a poisson approximation since there are only 
    -- very few solution points compared to the number of trials
    --
    -- we then use the normal approximatiom of the poissondistribution
    -- to calculate the confidence interval
    -- 
    -- (jk): to use the normal approximation of the poisson distribution numPoints should be >30 ?
    --
    --
    -- for a Poisson distribution the standard deviation is
    -- the square root of the expected value
    --
    err := opts#"confidence"*sqrt(numPoints); 
    estimate := new Interval from   (  max(round(1,numPoints-err), 0.0001),    round(1,numPoints+err)  );
    return estimate;
);

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


--
-- estimateCodim() is not used anywhere. deprecated?
--
estimateCodim = method( Options => (options poissonEstimate) );

estimateCodim( ZZ, ZZ, ZZ ) := Interval => opts->
    ( trials, numPoints, fieldCardinality ) -> 
(
    -- hints:
    -- estimatedInterval = fraction of fp-rational points (interval estimate) 
    --                   ≈ r * (1/cardinality)^codim + higher order terms in 1/cardinality
    -- rules: log_b(x) = - log_(1/b)(x). 
    -- with this we get. log_(1/cardinality) gamma_P = -log_(cardinality) gamma_P 
    -- = -log_natural(gamma_p)/log_natural(cardinality)
    --
    estimatedInterval := (1/trials)*poissonEstimate(numPoints,opts);     
    logEst := new Interval from (-log estimatedInterval.max,-log estimatedInterval.min);
    codimEst := (1/log fieldCardinality)*logEst;
    return new Interval from (round(1,codimEst.min),round(1,codimEst.max));
);

TEST ///
    debug FiniteFieldExperiments
    FiniteFieldExperimentsProtect()
    assert (estimateCodim(11^2*10,10,11) == new Interval from (1.8,2.4))
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


Experiment = new Type of MutableHashTable;

load "./FiniteFieldExperiments/Interpolation.m2";

--
-- ExperimentData is the data part of an Experiment. 
--  Experiment data (findings) and Experiment methods are separated to make storing more simple (or possible at all)
--
ExperimentData = new Type of MutableHashTable;

PointData = new Type of MutableHashTable;

new PointData from HashTable := (ancestorType, pointData)->( 
    return pointData;
);




--------------------------------------------------------------------------------------------------------------
---- observedValues, realizedValues realizedProperties and realizedPointProperties do the same. Choose the preferred naming!
--------------------------------------------------------------------------------------------------------------
observedValues = method();
observedValues (Experiment) := Thing => (experiment)->
(
    return experiment.observedValues();
);


realizedValues = method();
realizedValues (Experiment) := Thing => (experiment)->
(
    return experiment.realizedValues();
);


realizedProperties = method();
realizedProperties (Experiment) := Thing => (experiment)->
(
    return experiment.realizedProperties();
);


realizedPointProperties = method();
realizedPointProperties (Experiment) := Thing => (experiment)->
(
    return experiment.realizedPointProperties();
);

--------------------------------------------------------------------------------------------------------------
---- realizedValue, realizedProperty realizedPointProperty do the same. Choose the preferred naming 
--------------------------------------------------------------------------------------------------------------
realizedValue = method();
realizedValue (Experiment,ZZ) := Thing => (experiment,pos)->
(
    return experiment.realizedValue(pos);
);


realizedProperty = method();
realizedProperty (Experiment,ZZ) := Thing => (experiment,pos)->
(
    return experiment.realizedProperty(pos);
);


realizedPointProperty = method();
realizedPointProperty (Experiment, ZZ) := Thing => (experiment,pos)->
(
    return experiment.realizedPointProperty(pos);
);



--------------------------------------------------------------------------------------------------------------
---- realizedValue, realizedProperty realizedPointProperty do the same. Choose the preferred naming 
--------------------------------------------------------------------------------------------------------------



-- return the statistics for the watched properties : as Counts (derived from Tally)
-- where the key is the realized value set of properties and the value is the number of observed occurances
counts = method(); 
counts (Experiment) := Counts => (experiment)->
(
    return experiment.counts();
);

-- counts the number of collected points. Return type is Counts (derived from Tally and Tally from HashTable)
-- where the key is the realized value set of properties and the value is the number of collected points
collectedCount = method(); 
collectedCount (Experiment) := Counts => (experiment)->
(
    return experiment.collectedCount();
);

--------------------------------------------------------------------------------------------------------------
---- realizedValue, realizedProperty realizedPointProperty do the same. Choose the preferred naming 
--------------------------------------------------------------------------------------------------------------



points = method();

points(Experiment) := Thing => (experiment)->
(
    return experiment.points();
)


-- should return a points with their properties. (a HashTable where points are the keys)
richPoints = method();
 
richPoints (Experiment) := Thing => (experiment)->
(
    return experiment.pointsWithProperties();
);

-- same as richPoints 
pointsWithProperties = richPoints;


-------------------------------------------

createPointData = (pBlackBox, point)->
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


new ExperimentData from Ring := (E, coeffRing) -> 
(
     ffelog.debug (toString E);
     e := new MutableHashTable;
     e.coefficientRing = coeffRing;
     e.pointData = new MutableHashTable; -- (JK) unfortunate naming
    -- format: key=>{ideal, maxDegree, name} --later: data type for interpolated ideal (RichInterpolatedIdeal?)
     e.countData = new Counts;
     e.trials = 0;
     e.propertyList = {};
     e.isRandom = null;
     return e;
);
 
-- estimateNumberOfComponents()
--
-- estimates the interval for number of components for points with given set of properties 
-- uses jacobian rank as codim estimator.
-- the interval is estimated using poissonEstimate
--
-- (jk)  later :m'estimateNumberOfComponents'  could be encapsulated in a ComponentEstimator. 
--       Purpose: replace componentNumberEstimator if required (we could have different estimators)
--
-- interface: ComponentNumberEstimator.estimate(realizedPointPropertyTuple, estimatedCodimension) ;
--
estimateNumberOfComponents(Experiment,List) := HashTable => opts->    (experiment,realizedPointPropertyTuple) -> 
(
    countData := experiment.counts();
    posRankJacobianAt := experiment.position( experiment.usedRankJacobianAt() );
    if posRankJacobianAt === null then error("To estimate number of components, \"rankJacobianAt\" must be watched");
    
    cardinality := experiment.coefficientRingCardinality();
    
    jacobianRank := realizedPointPropertyTuple#posRankJacobianAt;
    estimatedCodimension:= jacobianRank;
    
    numPointsWithGivenPropertyTuple := countData#realizedPointPropertyTuple;
    
    return estimateNumberOfComponents(  experiment.trials(),
                                        estimatedCodimension,
                                        numPointsWithGivenPropertyTuple,
                                        cardinality,
                                        opts
                                     );
)

createExperimentData = (coeffRing,pointData,countData, trials, propertyList,isRandom) -> (
    e := new ExperimentData;
    e.coefficientRing = coeffRing;
    e.pointData = new MutableHashTable from pointData;
    e.countData = countData;
    e.trials = trials;
    e.propertyList = propertyList;
    e.isRandom = isRandom;
    return e;
)

new ExperimentData from ExperimentData := (E,ed) -> 
(
    -- print (toString E);
    e := new MutableHashTable;
    e.coefficientRing = ed.coefficientRing;
    e.pointData = copy ed.pointData;
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
            keys ed1.pointData      == keys ed2.pointData and    
            keys ed1.countData      == keys ed2.countData and  
            ed1.isRandom        == ed2.isRandom  and
            ed1.trials        == ed2.trials  ) then 
    (
        for key in keys ed2.pointData do
        ( 
            if not ed1.pointData#key == ed2.pointData#key then 
            return false;
        );

        for key in keys ed2.pointData do
        ( 
            if not ed1.countData#key == ed2.countData#key then 
            return false;
        );
        return true;
    )
    else return false; 
);

--coeffRing,pointData,countData,trials,propertyList,isRandom

toExternalString(ExperimentData) := String=> (ed)->
( 
    return "createExperimentData " | "(" 
    | toExternalString ed.coefficientRing | ", \n "
    | toString (new HashTable from ed.pointData) | ", \n"
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
        
        for key in keys ed2.pointData do
        ( 
        if not edNew.pointData#?key then 
            edNew.pointData#key = copy ed2.pointData#key
        else
            edNew.pointData#key = edNew.pointData#key | copy ed2.pointData#key;
        );
        edNew.countData = ed1.countData + ed2.countData;
        edNew.trials = ed1.trials + ed2.trials;
        return edNew;
    )
    else error ("+: experiment data not compatible");                   
);


PointIterator = new Type of HashTable;
new PointIterator from Thing := (E, thing) -> 
(
    error("not impelemted");
);

RandomPointIterator = new Type of PointIterator;
new RandomPointIterator from Thing := (E, thing) -> 
(
    error("not impelemted");
);


createRandomPointIterator = method();

createRandomPointIterator (Function) := RandomPointIterator =>( weakRandomPointGenerator )->
(
    -- todo improvement: use own seek and remember initial seek to get reproducible    
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
        
    rpi.setPosition =(trials)->
    (
        currTrial= trials;
    );
    
    -- usually only for loading from file.    
    rpi.setPoint = (point)->
    (
        randomPoint = point;
    );
    
    rpi.reset = () ->
    (
            randomPoint = null;
            currTrial = 0;
    );

    rpi.begin = ()->
    (
        -- jk why do we this???
        ri := createRandomPointIterator(weakRandomPointGenerator);
        ri.next();
        return ri;
    );

    rpi.point = ()-> randomPoint;
    rpi = new HashTable from rpi;
    rpi = newClass(RandomPointIterator,rpi);
    return rpi;
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
    
    rpi = new HashTable from rpi;
    rpi = newClass(RandomPointIterator,rpi);
    return rpi;
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
createIterator (List) := PointIterator =>( pPoints )->
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

    pIterator.next = ()->
    (
        if (currPosition+1 >= pointCount) then return false;
        point = points#currPosition;
        currPosition = currPosition + 1;
        return (true) ;
    );

    pIterator.begin = ()->
    (
        localPointIterator := createIterator( pPoints);
        localPointIterator.next();
        return localPointIterator;
    );

    pIterator.point = ()-> point;
    pIterator = new HashTable from pIterator;
    pIterator = newClass(PointIterator,pIterator);
    return pIterator;
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

-- todo 1. why do we not save the experiment?
-- todo 2: createEstimatedDecomposition should be an internal function.

createEstimatedDecomposition := (bbi, estimate, watchedProperties, usedConfidence)->
(
    ht :=   new HashTable from {  
        "blackBox" =>bbi,
        "estimate"=> estimate,
        "watchedProperties" =>watchedProperties,
        "usedConfidence" => usedConfidence
    
    };
    return new EstimatedDecomposition from ht;
)


net (EstimatedDecomposition) := Net =>(ed)->
(
   formatinfo :=  {
                    -- net "-- format: ",
                        net ("--(estimated codim, estimated number of components [confidence interval " | toString ed#"usedConfidence" |"*σ]) <= "
                            | toString ed#"watchedProperties" | " )")
                    };
    estimate := ed#"estimate";

    -- estimate#0 is (estimated codim, estimated number of components [confidence interval])
    -- and estimate#1 is the list of watched properties.
    sss := sort apply (estimate, entry ->  net( entry#0)  | " <= " | net (entry#1 ) );

    sss = formatinfo |sss;
    return net stack sss;
)


estimateDecomposition = method (Options =>  (options poissonEstimate));

  estimateDecomposition (Experiment) := EstimatedDecomposition  => opts -> (experiment) ->
    (
        posRankJacobianAt := experiment.position(  experiment.usedRankJacobianAt() );
        if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");

        estimate := flatten apply(keys experiment.counts(), 
                                          valuesOfProperties-> (
                                                           (valuesOfProperties#posRankJacobianAt,
                                                            estimateNumberOfComponents(experiment,valuesOfProperties, opts)
                                                           ),
                                                           valuesOfProperties
                                                                )
                                    );
        
        return createEstimatedDecomposition(experiment.blackBoxIdeal(), estimate, experiment.watchedProperties(), opts#"confidence" );
    );


createDecompositionEstimator = method(Options => (options poissonEstimate));


createDecompositionEstimator (Experiment) := Thing => opts-> (experimentParameter) ->
(
    decompEstimator := new MutableHashTable;    
               
    decompEstimator.estimateDecomposition =  (experiment) -> 
    (
        posRankJacobianAt := experiment.position(  experiment.usedRankJacobianAt() );
        if posRankJacobianAt === null then error("To estimate the decomposition, \"rankJacobianAt\" must be watched");

        -- estimate data is  a list of pairs:
        -- first entry is (rank at points, interval for estimated number of components)
        -- and second entry is the value array of watched properties.
        estimate := flatten apply(keys experiment.counts(), valuesOfWatchedProperties-> (
                                                                    (valuesOfWatchedProperties#posRankJacobianAt,
                                                                    estimateNumberOfComponents(experiment,valuesOfWatchedProperties, opts)),
                                                                    valuesOfWatchedProperties
                                                                )
                                    );
        
        return createEstimatedDecomposition(experiment.blackBoxIdeal(), estimate, experiment.watchedProperties(), opts#"confidence" );
    );
    decompEstimator = new HashTable from decompEstimator;
    
    return decompEstimator;     
);


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

possibleCodimComponents = (found, trials, orderK) -> 
(
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
round1withZeros = (rr) -> 
(
    if rr<0 then error;
    if rr==0 then return "0.0";
    if rr<1 then return "0"|toString(round(1,rr));
    if round(1,rr) == round(0,rr) then return toString(round(1,rr))|".0";
    return toString(round(1,rr))
);


--deprecated?
roundCodim = (trials, found, orderK) -> 
(
    estimate := 1/trials*poissonEstimate(found);
    if (log(estimate.max)/log(orderK) - log(estimate.min)/log(orderK))<0.8 
    then return round1withZeros((log(trials)-log(found))/log(orderK))
    else return "..."
);
 
 estimateStratification = method();
 
estimateStratification (Experiment) := EstimatedStratification =>  (experiment) -> 
(
    trials := experiment.trials();
    orderK := experiment.coefficientRingCardinality(); -- this must be read from the experimentdata
    -- (jk): need more advice. Did we want to use a different ring for search that the ideal coefficient ring? If so, 

    countData := experiment.counts();

    -- sort keys by number of occurence (descending)
    sortedKeysByOccurence := apply(reverse sort apply(keys countData,k->(countData#k,k)), i->i#1);

    --apply(sortedKeysByOccurence,k->(
    --      --print (net((log(trials)-log(i#0))/log(charK))|" <= "|net(i#1)));
    --      print (net(round(1,(log(trials)-log(countData#k))/log(orderK)))|" <= "|net(k)))
    --      );

    estimate := flatten apply(sortedKeysByOccurence,k->( (round1withZeros((log(trials)-log(countData#k))/log(orderK))), (k)) );
    --estimate := flatten apply(sortedKeysByOccurence,k->( (roundCodim(trials,countData#k,orderK))), (k)) );
    return  createEstimatedStratification(experiment.blackBoxIdeal(), estimate, experiment.watchedProperties() );
);




stratificationIntervalView := (stratificationData )->
(
    prerearrangedData := new HashTable from apply( keys stratificationData, key -> (stratificationData#key=>(stratificationData#key,key )));
    toSort := apply( keys stratificationData, key -> (stratificationData#key));
    sorted := sort (toSort); 
    rearrangedData := apply (sorted, key-> (prerearrangedData#key) );
    return new HashTable from rearrangedData;
);



   --#??? fuer den Fall dass dich die Schlüssel nicht sortieren lassen.

   
-- todo: naming 'countsByCount' is still unfortunate

   -- countsByCount() sorts the statistic by Number
   
   
sortedCounts = method();
sortedCounts ( Experiment ) := (experiment)->
(
    return experiment.sortedCounts();
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


 
 


Experiment + Experiment := Experiment => (re1, re2)->
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


watchedProperties = method();

watchedProperties (Experiment) := List => (experiment)->
(
    return experiment.watchedProperties();
);




--purpose of new Type Counts: to install sort as a method for Counts, but not for all Tally types.
--
Counts = new Type of Tally; 


--net (Counts) := Net =>(cs)->
--(
--    strcountinfo := "-- count() structure: \n-- values of watched properties" => count ";   
--    sss = stac (strcountinfo, net cs);
--    return net sss;
--)


SortableCounts = new Type of List;

new SortableCounts from Counts := (E, thing) -> 
(
    sortableCounts:= apply (keys thing, key->(key,thing#key));
    return sortableCounts;
)


net (SortableCounts) := Net =>(cs)->
(
    strcountinfo := "-- count structure: (values of watched properties) =>  count ";
    L := apply(cs, entry->  stack ( horizontalJoin(net entry#0 ," => " , net entry#1), " "));    
    L2 := stack L;
    L3 := horizontalJoin(toString class cs,"{",L2, "}");
    L4 := stack (strcountinfo, L3);
    return net L4;
)

-- sortfkt is a either sort or rsort.
-- internal method
--
sortCounts := method();
sortCounts (SortableCounts, MethodFunctionWithOptions) := SortableCounts => (countDataList, sortfkt)->
(    
    counts := new Counts from countDataList;
    
    --print (toString counts);
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
    sorted := sortfkt (toSort); 
    rearrangedData := {};

    -- 3.  asseble the result
    for count in sorted do
    (
        for entry in prerearrangedData#count do
        --rearrangedData = rearrangedData | {(count,entry)};
        rearrangedData = rearrangedData | {(entry,count)};
    );
    --print (toString rearrangedData);
    --return new List from rearrangedData;
    return new SortableCounts from rearrangedData;
);

sort (SortableCounts) := SortableCounts =>opts-> (countData)->
(
    return sortCounts(countData, sort);
);

rsort (SortableCounts) := SortableCounts =>opts-> (countData)->
(
    return sortCounts(countData, rsort);
);

-- sorts experiment statistics by count ascending
sort (Counts) := Counts =>opts-> (countData)->
(
    --  a HashTable is not sortable; => converting counts to a list of pairs (realizedpropertyValueTuple, count) 
    countDataList := new SortableCounts from countData;
    return sortCounts(countDataList, sort);
);

-- sorts experiment statistics by count descending

rsort (Counts) := Counts =>opts-> (countData)->
( 
    --  a HashTable is not sortable; => converting counts to a list of pairs (realizedpropertyValueTuple, count) 
    countDataList := new SortableCounts from countData;
    return sortCounts(countDataList, rsort);
);

new Experiment from BlackBoxParameterSpace := (E, pBlackBox) -> 
(
    blackBoxIdeal := pBlackBox;    --black box ideal or black box parameter space so far
    experimentData := new ExperimentData from blackBoxIdeal.coefficientRing;

    coefficientRingCardinality := ringCardinality( blackBoxIdeal.coefficientRing );

    experiment := new MutableHashTable;
    experiment = newClass( Experiment, experiment );

    -- todo: issue about naming here. blackBox (what?)
    --experiment.bbi = blackBoxIdeal; -- is this safe?
    --experiment.blackBoxIdeal = blackBoxIdeal;


    -- todo: maype return a (deep) copy 
    experiment.getExperimentData = ()->
    (
        ffelog.warning("-- warning : you got a reference to experiment data, do not modify! ");
        return experimentData;
    );

    -- some of the  following values initialized at end 
    -- ( e.g. propertiesAt initialization depends presence of some functions defined later)
    pointsPerComponent := 10;

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
    
    decompositionConfidenceInterval := 1.96; -- value times sigma (measured in sigma)
    
    -- later the decompositionEstimator should be public and modifiable 
    -- (provide a method to modify the confidence parameter)    
    --decompositionEstimator := createDecompositionEstimator(experiment,"confidence"=>1.96);

    experiment.estimateDecomposition = () -> 
    (
        estimateDecomposition(experiment, "confidence"=> decompositionConfidenceInterval)
    );
    
    --experiment.estimateDecomposition = () -> (decompositionEstimator.estimateDecomposition(experiment));

    -- (jk) later setDecompositionConfidenceInterval must be dropped, because we never may know the config parameters for   
    -- a custom decomposition parameter (maybe provided by the user => learn the plugin interface!)
    -- the confidence interval is measured in sigma (not in percent)
    --
    experiment.setDecompositionConfidenceInterval = (confidenceInterval)->
    (
        decompositionConfidenceInterval = confidenceInterval;
    );

    experiment.estimateStratification = () -> (estimateStratification(experiment));

--   experiment.estimateStratification2 = () -> ( estimateStratification2(experiment) );

    experiment.compatible = method();
    experiment.compatible (Experiment) := Boolean =>(re2)->
    (
        bb1 := experiment.blackBoxIdeal();
        bb2 := re2.blackBoxIdeal();

     return (  experiment.pointsPerComponent() ==  re2.pointsPerComponent() and
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
    experiment.merge (Experiment)  := Experiment => (re)->
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

    experiment.coefficientRing = ()->
    (
        
        return experimentData.coefficientRing;
    );
    
    -- todo: rename to sortedCounts or similar
    experiment.countsByCount = ()->
    (   
        --print ("--warning countsByCount deprecated. Use 'sort counts experiment' or 'sort experiment.counts()' or  experiment.sortedCounts()");
        return sort new SortableCounts from experimentData.countData ; 
    );
    
    -- same as countsByCount
    experiment.sortedCounts = ()->
    (           
        return sort  experimentData.countData ; 
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

    experiment.points = ()->
    (
        return flatten  apply (keys experimentData.pointData, propertiesAsKey -> experimentData.pointData#propertiesAsKey);
    );
    
  

    experiment.smoothPoints = (precision, trials)->
    (
        plist :=  experiment.points();
        smothPoints :=  {};
        jet := null;
        for point in plist do
        (
                jet = jetAt(blackBoxIdeal,point,precision,trials);
                if (jet#"succeeded") then
                (
                    smothPoints = smothPoints | {point};
                );
                
        );
        return smothPoints;
    );



    runExperimentOnce := method();


    experiment.setPointsPerComponent = (numPointsPerComponent)->
    ( 
        pointsPerComponent = numPointsPerComponent;
    );


    experiment.pointsPerComponent = ()->
    ( 
        return pointsPerComponent ;
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
        experimentData.pointData = new MutableHashTable;
        experimentData.countData = new Counts;
    );
    
    experiment.rankJacobianAtKey = ()->
    (
        return rankJacobianAtKey;
    );

    experiment.usedRankJacobianAt = ()->
    (
        return rankJacobianAtKey;
    );

    -- question: why the hell we pass the number of wanted points per component and do not acces the 
    -- experiment => to see the dependency
    --
    runExperimentOnce(ExperimentData, Matrix,ZZ ) := Tally => (experimentData, point, wantedPointsPerComponent) -> 
    (  
        K := experimentData.coefficientRing;
        --prime = char K
        numVariables := blackBoxIdeal.numVariables;
        wantedPoints := wantedPointsPerComponent;

        -- if ideal vanishes on the random point do something
        if experiment.isInteresting(point) then   
        (
            valuesTuple := propertiesAt(point); -- todo: valuesTuple: better naming?
 
            -- countData number of found points for each rank and property
            experimentData.countData = experimentData.countData + new Counts from tally {valuesTuple};
          
            if  rankJacobianAt =!= null then 
            (
                FFELogger.debug( "update wanted points" );
                rankJacobian := rankJacobianAt(point);  
                upperEstimate := (estimateNumberOfComponents(experiment,valuesTuple)).max;
                ffelog.debug ("upper estimate number of components:  " | toString upperEstimate );
                --upperEstimate := 1; -- test

                wantedPoints = max(1,upperEstimate)*wantedPointsPerComponent;
            );            

            -- remember some points
            if experimentData.pointData#?(valuesTuple) then 
            (
                -- collect a fixed number of points per estimated component
                -- use upper limit of estimation for this             
                if #(experimentData.pointData#valuesTuple) < wantedPoints then 
                (
                    FFELogger.debug( "attaching point" );
                    experimentData.pointData#valuesTuple = experimentData.pointData#valuesTuple | {point};
                );
            )
            else (
                FFELogger.debug( "attaching first point for some valuesTuple key");
                experimentData.pointData#valuesTuple = {point};
            );
        );
        -- this trial is counted in runExperiments and not here to allow update() without changing trial num.
    );

   
    experiment.setPointIterator  = (pRpi)->
    (
        if experiment.trials()=!=0 then 
            error ("cannot change point iterator - experiment was already run! You could call clear() and retry.");
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
            runExperimentOnce( experimentData, pPointIterator.point(), pointsPerComponent );
            experimentData.trials =  experiment.trials();
        );
    );


    positionsOfProperties := (newWatchedList)->
    ( 
        return  apply( newWatchedList, 
                        targetProperty -> position( experiment.watchedProperties() , (prop)->prop==targetProperty)
            );
    );

    updateWatchedListIsProjection := (newWatchedList)->
    (
        propertyPositions :=  positionsOfProperties(newWatchedList);

        return( 0 == #(select(propertyPositions, (pos)->pos===null)) ); 
    );

    projectionUpdate :=  ( newWatchedList) -> 
    (     
        -- assume: new watched property list is projection of old one.

        -- 1. get the tuple of positions of new watchedProperties in propertyList.
        propertyPositions :=  positionsOfProperties(newWatchedList);


        assert( 0 == #(select(propertyPositions, (pos)->pos===null)) );        

        --    experimentData.propertyList 

        newCountData := new MutableHashTable;

        newKey := null;

        -- sum up counts
        for key in keys experimentData.countData do
        (
            newKey = apply( propertyPositions, pos-> key#pos);
            if (not newCountData#?newKey) then  (    newCountData#newKey =   experimentData.countData#key     )
                                            else (    newCountData#newKey =   newCountData#newKey + experimentData.countData#key     );

        );
        -- reclassify collected points

        newPoints := new MutableHashTable;

        for key in keys experimentData.pointData do
        (
            newKey = apply( propertyPositions, pos-> key#pos);
            if (not newPoints#?newKey)    then  (    newPoints#newKey =   experimentData.pointData#key     )
                                            else (    newPoints#newKey =   newCountData#newKey |  experimentData.pointData#key     );

        );

        experimentData.countData = new Counts from new Tally from newCountData;
        experimentData.pointData    = newPoints;
        
        --experimentData.propertyList = newWatchedList;
    );


    setWatchedPropertiesInternal := (propListToObserve)->
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

    experiment.useRankJacobianAt = (rankJacobianAtName)->
    ( 
        if experiment.trials()=!=0 then 
            error ("cannot change rankJacobianAt  - experiment was already run! You could clear() the statistics and retry. ");

        if rankJacobianAtName===null then
        (
            rankJacobianAtKey = null;
            rankJacobianAt = null ;
            return;
        );

        if (not blackBoxIdeal.hasPointProperty(rankJacobianAtName)) then 
            error ("blackBoxIdeal seems not to have property " | rankJacobianAtName );
    
        if (rankJacobianAtKey=!=null) then
        (
            if (experiment.propertyIsWatched(rankJacobianAtKey)) then
            (
                experimentData.propertyList = delete(rankJacobianAtKey, experimentData.propertyList ) ;
                setWatchedPropertiesInternal( experimentData.propertyList );   
            );
        );

        rankJacobianAtKey = rankJacobianAtName;
        rankJacobianAt = blackBoxIdeal.pointProperty(rankJacobianAtName) ;
        experiment.watchProperty(rankJacobianAtKey);
    
    );



    experiment.clearWatchedProperties = (   )->
    (
        experiment.clear(); 
        setWatchedPropertiesInternal({});
    );

    experiment.reset = (   )->
    (
        --redundant: experiment.clear(); 
        experiment.clear(); 
        clearWatchedProperties();
    );

    --newPropertyListIsProjection = (newPL)->
    --(
    --     
    --);

    experiment.setWatchedProperties = ( propertyStringList )->
    (  

        if ( updateWatchedListIsProjection(propertyStringList) ) then 
        (
            projectionUpdate(propertyStringList);
            setWatchedPropertiesInternal(propertyStringList);
        )
        else
        (

            if experiment.trials()=!=0 then error ("cannot change watched properties - experiment was already run! Clear statistics and retry.");

            -- improvement: allow in case the new properties is a projection
            setWatchedPropertiesInternal(propertyStringList);
        );
    );

    UpdateWatchedPropertiesError := "cannot change watched properties - experiment was already run! You could clear() the statistics and retry.";
  

    experiment.propertyIsWatched = method();
    experiment.propertyIsWatched (String) := Boolean => (propertyName)->
    (
        if ( #(select( experiment.watchedProperties(), (prop)->propertyName==prop)) >0 ) then
        ( 
            return true;
        );    
        return false;
    );

    watchProperty := method();
    watchProperty ( String ) := List => (propertyName)->
    (
        if not blackBoxIdeal.hasPointProperty(propertyName) then 
            error ("blackBoxIdeal seems not to have property" | propertyName );

        if (not experiment.propertyIsWatched(propertyName)) then
        (
            if (experiment.trials()=!=0 ) then error (UpdateWatchedPropertiesError);
        );

        experimentData.propertyList = unique (experimentData.propertyList | { propertyName }) ;
        setWatchedPropertiesInternal( experimentData.propertyList );   
      
    );
    
    watchProperty ( Function ) := List => (propertyMethod)->
    (
         watchProperty(toString propertyMethod);
    );
    
    experiment.watchProperty = method();
    experiment.watchProperty ( String ) := List => (propertyName)->
    (
        watchProperty(propertyName);
        return experiment.watchedProperties();
    );
    
    experiment.watchProperty ( Function ) := List => (propertyFunction)->
    (
        return experiment.watchProperty(toString propertyFunction);
    );

    experiment.watchProperties = (propertyNameList)->
    (
        for propertyName in propertyNameList do
        (
            watchProperty(propertyName);
        );
        return experiment.watchedProperties();
    );

    

    assertPropertyIsWatched := (propertyName) ->
    (
        if (not experiment.propertyIsWatched(propertyName)) then 
            error ("given property '" |propertyName| "' is not watched !");
    );

    assertIgnorePropertyAllowed := (propertyName) ->
    (
        if (propertyName == "rankJacobianAt") then
        (
            if (blackBoxIdeal.type===BlackBoxIdeal) then
            (
                error (" removing 'rankJacobianAt' from watched properties for a " | toString blackBoxIdeal.type | " not allowed !")
            );
        );
    );


    experiment.ignoreProperty = (propertyName)->
    (
        
        assertPropertyIsWatched(propertyName);
        assertIgnorePropertyAllowed(propertyName);

    
        newPropertyList := experimentData.propertyList ;
    
        newPropertyList = delete(propertyName, newPropertyList ) ;
    
        projectionUpdate( newPropertyList );

        setWatchedPropertiesInternal( newPropertyList );   
    );


    experiment.ignoreProperties = (ignorePropertyStringList)->
    (
        for propertyName in ignorePropertyStringList do
        (
            assertPropertyIsWatched(propertyName);
            assertIgnorePropertyAllowed(propertyName);
        );

        newPropertyList := experimentData.propertyList ;
        apply( ignorePropertyStringList, propToIgnore-> ( newPropertyList= delete(propToIgnore, newPropertyList ); ));
        projectionUpdate( newPropertyList );
        setWatchedPropertiesInternal(newPropertyList);    
    );

 
    experiment.interpolateComponents = (maxDeg, onComponentPrecision)->
    (
            blackBoxIdeal.interpolateComponents( experiment.smoothPoints(10,10), maxDeg, onComponentPrecision);
    );


    experiment.watchedProperties =  ()->
    (
        return experimentData.propertyList;
    );

    
    experiment.setIsInteresting = (pIsInteresting)->
    (  
       if experiment.trials()=!=0 then 
            error ("cannot change isInteresting - experiment was already run! You could call clear() and retry.");
        if ( pIsInteresting=!=isInteresting ) then
        (
            isInteresting = pIsInteresting;
        );
    );

    ---- syntax for the moment too hard (method without parameters)
    --experiment.runExperimentOnce = method(Options => (options runExperimentOnce));
    -- experiment.runExperimentOnce() := Thing => opts->()->
    --(
    --   return runExperimentOnce(experimentData);
    --);  
  

    experiment.run = method();
    experiment.run(ZZ) := Thing => (newTrials)->
    (
        runExperiment( experimentData, pointIterator, newTrials );
        return sort experiment.counts();
    );
    
    --
    -- returns a hashtable with 
    -- 
    experiment.pointLists = ()->
    (
        return new HashTable from experimentData.pointData;
    );

    experiment.pointsByKey = (key)->
    (
        if not (experimentData.pointData)#?key then 
            error "invalid key";
        return  (experimentData.pointData)#key;
    );
    
      -- returns a HashTable with watched BlackBoxIdeal properties as keys and corresponding occured count as values.
    -- naming is unfortunate
    experiment.counts = ()->
    (
        return new Counts from experimentData.countData;
    );    
    
    -- count the number of collected points!
    experiment.collectedCount = ()->
    (
        return new Counts from apply( experiment.pointKeys(), key->( key=> #(experimentData.pointData)#key ) );
    );
    
    
    experiment.realizedPointProperties = ()->
    (
        return keys experimentData.pointData;
    );


    -- deprecated
    experiment.pointKeys = ()->
    (
        print("--warning: .pointKeys() is deprecated. use .realizedProperties()");
        return experiment.realizedPointProperties();
    );
    
    experiment.realizedProperties = experiment.realizedPointProperties;    
    experiment.realizedValues = experiment.realizedPointProperties;
    experiment.observedValues = experiment.realizedPointProperties;
    
    experiment.realizedProperty = method();
    experiment.realizedProperty (ZZ) := Thing =>(index)->
    (
        propertyvalues := experiment.realizedProperties();
        if ( (index<0) or (index>= #propertyvalues) ) then error "invalid point position ";
        return propertyvalues#index;
    );
    
    experiment.pointKey = method();
    experiment.pointKey(ZZ) := Thing => (index)->
    (
       print("--warning: .pointKey() is deprecated. use .realizedProperty()");
       return experiment.realizedProperty(index); 
    );
    
   
    
    experiment.realizedValue = realizedProperty;    
    experiment.observedValue = realizedProperty;    
    experiment.realizedPointProperty = realizedProperty;

   

   --init part:

        -- todo: test if changing blackBoxIdeal.rankJacobianAt is transparent
        -- (means after blackBoxIdeal.updatePointProperty("rankJacobianAt") the new one is called)
 
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

    experiment.tryProperty = (tProp) -> 
    (
        ffelog.info("-- ( " | toString experiment.watchedProperties() |" | " | tProp | " ) => count " );
        pointListsCopy := experiment.pointLists();
        tally flatten apply(keys pointListsCopy, --keys are the value tuples of the watched properties,
                       valueTuple->apply(pointListsCopy#valueTuple, -- for a value tuple set we have a list of stored points.
                                     point->( valueTuple, (blackBoxIdeal.pointProperty(tProp))(point))
                                 ) 
                 )
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
        apply( keys experimentData.pointData, key-> 
                                        ( (experimentData.pointData)#key = apply( (experimentData.pointData)#key , point->
                                                                                sub(point, experimentData.coefficientRing )
                                                                            )
                                        )
                );
    );
    
 
 

 
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
            With an @TO{Experiment}@ it is possible to check point properties of an @TO BlackBoxParameterSpace@ or @TO BlackBoxIdeal@ 
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
            \,\, \bullet \,{\tt pointsPerComponent}: \break
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
            \,\, \bullet \,{\tt setPointsPerComponent}:  \break
            \,\, \bullet \,{\tt setPointGenerator}: \break
            \,\, \bullet \,{\tt setPointIterator}:  \break
            \,\, \bullet \,{\tt useRankJacobianAt}:  \break
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
        e.counts()
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
           sort e.counts()            
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

            {\bf QuickStart } \break \break
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
        e.counts()
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
           sort e.counts()                     
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
           e.pointsPerComponent()
           e.collectedCount() 
        Text
           Here we have not collected exactly 10 points per component since the experiment uses the upper end of the confidence interval for the number of components ( see @TO estimateNumberOfComponents@) as guide for the number of points to keep.
           The amount of stored points can be adjusted:
        Example
           e.setPointsPerComponent(20)
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
    -- loadPackage "FiniteFieldExperiments"
    debug FiniteFieldExperiments
    FiniteFieldExperimentsProtect()
    coeffRing := ZZ/3;
    bbRankM = blackBoxParameterSpace( 5 ,coeffRing )
    rankMat := (point)->5
    bbRankM.registerPointProperty("rankJacobianAt",rankMat)


    point := matrix {{1,2,3,4,5}};
    point = sub( point, coeffRing);

 
    --bbRankM = bbRankM.getUpdatedBlackBox(bbRankM)
 

    e = new Experiment from bbRankM
    assert (e.coefficientRing()===coeffRing);

    e.setPointsPerComponent(20);
    assert( e.pointsPerComponent()==20);
    FFELogger.setLogLevel(4);
    e.watchProperties {"rankJacobianAt"};
    e.watchedProperties()
    assert( 1== # select( e.watchedProperties() , 
                       (prop)->(prop=="rankJacobianAt") ) 
     )
    e.useRankJacobianAt("rankJacobianAt");
    --e.useRankJacobianAt(null);
   
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
           e.counts()
           e.pointLists()
           e.watchedProperties()
        Text
           \break Now we clear the statistics and point lists:
        Example
           e.clear()
           e.trials()
           e.counts()
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
           e.counts()
           e.pointLists()
           e.watchedProperties()
        Text
           \break Now we clear the statistics, the point lists and the watched properties.
        Example
           e.clearWatchedProperties()
           e.trials()
           e.counts()
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
            are collected (see setPointsPerComponent). Therefore
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
        Text
           \break We set the initial random seed to produce deterministic results
           for the documentation:
        Example
           setRandomSeed(42); 
           e.run(100)          
           e.collectedCount()
           e.pointLists()
           e.pointsByKey({2})
           e.pointsPerComponent()
        Text
           \break Notice that the number of collected points can be larger than
           the number pointsPerComponent() since the experiment
           tries to estimate the number of components for each combination
           of properties. In the beginning where only a few points have
           been found the statistics might be so errorprone that some extra
           points are collected. 
   SeeAlso
          pointLists
          pointsByKey
          pointsPerComponent
          setPointsPerComponent             
///

doc ///
   Key
        "counts"
   Headline
        shows how often points with particular properties were found by an experiment
   Usage   
        e.counts()
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
            
            e.counts() is called automatically when e.run is finished.
            
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
           e.counts()
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
           sort e.counts()
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
   SeeAlso
      ignoreProperties
      watchProperty
      watchProperties
      watchedProperties  
///

doc ///
   Key
        "ignoreProperties"
   Headline
        deletes a property from the list of watched properties
   Usage   
        e.ignoreProperties(L)
   Inputs  
        e:Experiment 
            an Experiment
        L:List
            a list of property names. 
   Description
        Text
           This removes several properties from the list of watched properties.
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
   SeeAlso
      ignoreProperty
      watchProperty
      watchProperties
      watchedProperties   
///

doc ///
   Key
        "pointsPerComponent"
   Headline
        the number of points an experiment tries to collect on each component.
   Usage   
        e.pointsPerComponent()
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
           e.pointsPerComponent()
        Text 
           \break The number of point the experiment has collected
        Example      
           e.collectedCount()
           e.run(100)
           e.collectedCount()
        Text
           \break Lets now increase the number of points we want to collect
        Example
           e.setPointsPerComponent(20)
           e.pointsPerComponent()    
           e.run(100)
           e.collectedCount()
           e.run(100)
           e.collectedCount()
   SeeAlso
        setPointsPerComponent
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
        pointsPerComponent
        setPointsPerComponent
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
        pointsPerComponent
        setPointsPerComponent
        collectedCount        
        pointLists
        pointKeys
///        

doc ///
   Key
        "setPointsPerComponent"
   Headline
        change the number of points an experiment tries to collect on each component.
   Usage   
        e.setPointsPerComponent(num)
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
           e.pointsPerComponent()
        Text 
           \break The number of point the experiment has collected
        Example      
           e.collectedCount()
           e.run(100)
           e.collectedCount()
        Text
           \break Lets now increase the number of points we want to collect
        Example
           e.setPointsPerComponent(20)
           e.pointsPerComponent()    
           e.run(100)
           e.collectedCount()
           e.run(100)
           e.collectedCount()
   SeeAlso
        pointsPerComponent
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
           sort e.counts()
        Text
           We now want to look at the collected point with special
           betti tableaus
        Example
           e.pointKeys()
           e.pointKey(0)
           e.pointsByKey e.pointKey(0)
   SeeAlso
        pointLists
        pointsByKey
        collectedCount
        pointsPerComponent
        setPointsPerComponent
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
           
           For a more realistic example see @TO "Experiment for singularities of cubic surfaces" @
   SeeAlso
       watchedProperties
       watchProperty
       watchProperties
///                 

doc ///
   Key
        "usedRankJacobianAt"
   Headline
        the name of the black box property used to calculate the codimension of the tangent space at a point.
   Usage   
        e.usedRankJacobianAt()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
           Sometimes it is useful to implement a custom method for
           calculating the codimension of the tangent space at a point.

           This method must be first be registered as a property of the black box 
           used in the experiment. Then one tells the experiment to
           use this new property for caculating tangent spaces.
           
           This function documented here can then be used to see which property
           is currently used for calculating the codimension of a tangent
           space at a point.
            
           Lets see how this works in an example.                
        
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break We now create (in this case stupid) new property:
        Example
           rankAllways5At = (point) -> 5;
           bb = bb.rpp("rankAllways5At",rankAllways5At);
           bb.knownPointProperties()
        Text
           \break Now make an experiment from the blackbox:
        Example        
           e = new Experiment from bb;
           e.usedRankJacobianAt()
           e.run(100)
        Text 
           \break Now we change the method for calculating the
           rank of jacobi matrices. Before doing this we must
           clear the statistics.
        Example
           e.clear()
           e.useRankJacobianAt("rankAllways5At")
           e.usedRankJacobianAt()
           e.run(100)
        Text
           A more realistic application is for example the case
           where the equations of our ideal can be written 
           as A*B = 0 with A and B matrices with polynomial entries.
           We can then use the product rule 
           (A*B)' = A'*B + A*B' to differentiate. This is often faster.
           
           An other case where this might be used is when we have a morphism
           X -> Y and look at random points in X but are interested in the
           tangent space after projecting to Y.
   SeeAlso
      useRankJacobianAt
///                 

doc ///
   Key
        "useRankJacobianAt"
   Headline
       change the black box property used to calculate the codimension of the tangent space at a point.
   Usage   
        e.useRankJacobianAt(name)
   Inputs  
        e:Experiment 
            an Experiment
        name: String
            name of the new property to be used
   Description
        Text
           Sometimes it is useful to implement a custom method for
           calculating the codimension of the tangent space at a point.

           This method must be first be registered as a property of the black box 
           used in the experiment. Then one uses the function
           documented here to tell the experiment to
           use this new property for caculating tangent spaces.
           
           It is important to tell the experiment explicitly which property
           calculates the codimension of the tangenspace at a point since
           this is used in calculating estimates of decompositions. 
            
           Lets see how this works in an example.                
        
           First we create an ideal we want to analyse and put it into a blackbox:
        Example      
           K = ZZ/5;
           R = K[x,y,z];
           I = ideal (x*z,y*z);
           bb = blackBoxIdeal I;
        Text
           \break We now create (in this case stupid) new property:
        Example
           rankAllways5At = (point) -> 5;
           bb = bb.rpp("rankAllways5At",rankAllways5At);
           bb.knownPointProperties()
        Text
           \break Now make an experiment from the blackbox:
        Example        
           e = new Experiment from bb;
           e.usedRankJacobianAt()
           e.run(100)
           e.estimateDecomposition()
        Text 
           \break Now we change the method for calculating the
           rank of jacobi matrices. Before doing this we must
           clear the statistics.
        Example
           e.clear()
           e.useRankJacobianAt("rankAllways5At")
           e.watchProperty("rankAllways5At")
           e.usedRankJacobianAt()
           e.run(100)
           e.estimateDecomposition()
        Text
           A more realistic application is for example the case
           where the equations of our ideal can be written 
           as A*B = 0 with A and B matrices with polynomial entries.
           We can then use the product rule 
           (A*B)' = A'*B + A*B' to differentiate. This is often faster.
           
           An other case where this might be used is when we have a morphism
           X -> Y and look at random points in X but are interested in the
           tangent space after projecting to Y.
   SeeAlso
      usedRankJacobianAt
///                 
 

doc ///
   Key
        "watchProperty"
   Headline
        add a property the list of watched properties
   Usage   
        e.watchProperty(name)
   Inputs  
        e:Experiment 
            an Experiment
        name:String
            the name of the black box property to be added. 
   Description
        Text
           This add a property to the list of watched properties.
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
           at a random point, but also whether the point is probably smooth.
        Example
           e.watchProperty "isProbablySmoothAt"
           e.watchedProperties()
           e.run(500)
   SeeAlso
      ignoreProperty
      ignoreProperties
      watchProperties
      watchedProperties   
///
 
doc ///
   Key
        "watchProperties"
   Headline
        add a list of properties the list of watched properties
   Usage   
        e.watchProperties(L)
   Inputs  
        e:Experiment 
            an Experiment
        L:String
            of name of property to be added. 
   Description
        Text
           This add several properties to the list of watched properties.
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
           at a random point, but also whether the point is probably smooth.
        Example
           e.watchProperties({"isProbablySmoothAt","isCertainlySingularAt"})
           e.watchedProperties()
           e.run(300)
   SeeAlso
      ignoreProperty
      ignoreProperties
      watchProperty
      watchedProperties   
///
 
doc ///
   Key
        "watchedProperties"
   Headline
        returns the list of properties that are watched by an experiment
   Usage   
        e.watchedProperties()
   Inputs  
        e:Experiment 
            an Experiment
   Description
        Text
           During an experiment a black box is evaluated in random points.
           The experiment automatically looks at certain properties 
           associated to the point. This could be the codimension of 
           the tangent space at this point or wether the variety considered
           is singular or smooth.
           
           More interesting applications are possible if the variety
           considered is a parameter space of certain algebraic objects 
           (curves, matrices, etc). In this case the object parametrized
           by a point can have interestering properties wich we want to 
           study (eg. number of singularities in the case of curve or
           betti numbers of their kernels in the case of matrices).
      
           The experiment keeps automatically counts the number of
           times each combination of properties was encountered
           during the experiment so far. This gives some heuristic
           information about the size of the strata in which these occur.
           
           The function documented shows which properties are currently
           watched.
           
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
           \break rankJacobianAt is automatically watched since the rank of
           the Jacobi matrix is the codimension of the tangent space to
           the variety at this point. This gives an upper bound to the
           codimension of the component of the variety on which the 
           point lies. This is used in the heuristic estimates given
           by this package. 
        Example
           e.run(200)
           e.estimateDecomposition()
        Text   
           \break Lets now also watch whether the points considered
           are probably smooth. For this the statitics have to be cleared
           first:
        Example
           e.clear()
           e.watchProperty("isProbablySmoothAt")
           e.watchedProperties()
           e.run(300)
   SeeAlso
      ignoreProperty
      ignoreProperties
      watchProperty
      watchProperties
///

doc ///
   Key
        "setPointIterator"
   Headline
        changes the way random points are chosen
   Usage   
        e.setPointIterator(iterator)
   Inputs  
        e:Experiment 
            an Experiment
        itearator: HashTable
            defining a point iterator            
   Description
        Text
           Sometimes it is usefull to change the way random points of an 
           experiment are chosen. 
           
           Lets look at a typical but somewhat artificial example:           
        Example
           K = ZZ/11;
           R = K[a,b,c,x,y,z];
           I = ideal (a*x^2+b*y^2+c*z^2,a*x+b*y+c*z);
           bb = blackBoxIdeal I;
        Text     
           Notice that the equations of the ideal are linear in a,b,c.
           Therefore to find a point in the vanishing set, 
           we can choose x,y,z randomly first and solve the resulting
           linear equations for a,b,c later:
        Example
           pointX = random(K^1,K^3)
        Text    
           We substitute these values into the equations of the
           ideal to get linear equations for a,b,c:
        Example
           gensA = matrix{{a,b,c}}
           IA = sub(I,gensA|pointX)
        Text
           We now extract the coefficients of these linear equations
        Example
           coefficientsLinearEquation = transpose sub(diff(transpose gensA,gens IA),K)
        Text
           and look at a basis of the linear subspace defined by these equations
        Example
           basisOfLinearSpace = syz coefficientsLinearEquation
        Text
           Notice that this space can sometimes be more than 1 dimensional
           if the coefficient vectors are linearly dependent.
           
           Now we choose a random point inside the linear space by
           taking a random linear combination of the basis vectors
        Example
           pointA = transpose (basisOfLinearSpace * random(source basisOfLinearSpace,K^1))
        Text
           The values for a,b,c and for x,y,z give our random point
        Example
           point = pointA|pointX
        Text
           indeed it lies in the vanishing set of our ideal:
        Example
           sub(I,point)
        Text
           Lets put this into a function:
        Example
           smartRandomPoint = () -> (
                -- choose x,y,z randomly
                pointX = random(K^1,K^3);
                -- substituting these values we get linear equations for a,b,c
                IA = sub(I,gensA|pointX);
                -- we now extract the coefficients of these linear equations
                coefficientsLinearEquation = transpose sub(diff(transpose gensA,gens IA),K);
                -- and look at a basis of the linear subspace defined by these equations
                basisOfLinearSpace = syz coefficientsLinearEquation;
                -- now choose a random point inside this linear space
                pointA = transpose (basisOfLinearSpace * random(source basisOfLinearSpace,K^1));
                point = pointA|pointX
                );
        Text
           Indeed this gives points in the vanishing set of I:
        Example
           sub(I,smartRandomPoint())
        Text
           We now want to use this function to choose random points
           in an experiment. For this we frist have to create an 
           iterator from it:
        Example
           smartIterator = createRandomPointIterator(smartRandomPoint);
        Text
           Now we create the experiment:
        Example
           e = new Experiment from bb;
           e.run(200)
        Text
           running the experiment with the usual random generator
           finds only a few points on the variety. (about 200/11^2 = 8 since
           two equations have to vanish).
      
           We now change the random generator to the function above.
           Before doing this we must clear the statistics since they
           would be useless if points found by different methods are mixed
        Example
           e.clear()
           e.setPointIterator(smartIterator)
           e.run(100)
        Text
           We see that many more points are found in the same time. 
           This reduces the time needed to find interesting points
           in high codimension.          
///



TEST ///
    -- test watching user defined rankJacobianAt
    K = ZZ/5;
    R = K[x,y,z];
    I = ideal (x*z,y*z);
    bb = blackBoxIdeal I;

    rankAllways5At = (point) -> 5;
    bb = bb.rpp("rankAllways5At",rankAllways5At);

    assert(bb.hasPointProperty("rankAllways5At"));
    e = new Experiment from bb;

    e.run(100)

    e.clear()
    e.useRankJacobianAt("rankAllways5At")
    assert(e.propertyIsWatched( "rankAllways5At"));
    assert(1 == #(e.watchedProperties()) );
    e.run(100)
    assert( 1 == #( e.collectedCount() ) );
///
 
doc ///
    Key
        "Experiment for singularities of cubic surfaces"
    Headline
        use an Experiment to study the space of cubic surfaces
    Description
        Text
            A black box parameter space is used to implement parameter spaces
            with their universal families in a pointwise fashion.
            
            Let us build the parameter space of cubic surfaces with a view of 
            studying its stratification with respect to singularity type.
            
            We work in charateristic 7.            
        Example    
            K = ZZ/7
        Text
            The coordinate ring of IP^3
        Example
            R = K[x,y,z,w]
        Text  
            Make an empty blackbox which will later contain our 
            describtion of the parameter space of cubic surfaces.
            It will depend on 20 parameters, since acubic polynomial 
            in 4 variables has 20 coefficients.
        Example
            bbC = blackBoxParameterSpace(20,K);
            bbC.knownPointProperties()
        Text
            We now build the cubics from the coefficents, i.e. we
            construct the member of the universal familiy over 
            a given parameter point:
        Example
            mons3 = matrix entries transpose super basis(3,R)
            cubicAt = (point) -> matrix entries (point*mons3)
        Text 
            register this function in the black box 
        Example
            bbC = bbC.registerPointProperty("cubicAt",cubicAt);
            bbC.knownPointProperties()
        Text
            Lets test this functionality with some special cubics.
            The first example is the cubic cone. It is singular
            at (0:0:0:1):                   
        Example    
            cubicCone = matrix{{x^3+y^3+z^3}}
            coeffCubicCone = contract(transpose mons3,cubicCone)
            bbC.cubicAt(coeffCubicCone)
        Text
            The second example is the Fermat cubic. It is smooth everywhere    
        Example
            cubicFermat = matrix{{x^3+y^3+z^3+w^3}}
            coeffCubicFermat = contract(transpose mons3,cubicFermat)
            bbC.cubicAt(coeffCubicFermat)
        Text
            Now we want to implement the stratification by singularity type.
            For this we first determine the singular locus of a cubic surface:
        Example
            singularLocusAt = (bb,point) -> ideal jacobian bb.cubicAt(point)
            bbC = bbC.rpp("singularLocusAt",singularLocusAt);
            bbC.knownPointProperties()
            bbC.singularLocusAt(coeffCubicCone)   
            bbC.singularLocusAt(coeffCubicFermat)
        Text
            As a first approximation of the singularity type we use
            the degree of the singular locus
        Example
            degreeSingularLocusAt = (bb,point) -> (
                 s := bb.singularLocusAt(point);
                 if dim s == 0 then return 0;                
                 if dim s == 1 then return degree s;                 
                 if dim s >= 2 then return infinity;
              )
            bbC = bbC.rpp("degreeSingularLocusAt",degreeSingularLocusAt);
            bbC.knownPointProperties()
        Text
            Calculate the degree of the singular locus for our examples
        Example
            bbC.degreeSingularLocusAt(coeffCubicCone)
            bbC.degreeSingularLocusAt(coeffCubicFermat)
        Text
            Now the BlackBoxParameterspace hat a number of point properties
        Example
            bbC.knownPointProperties()
        Text
            These properties can now be used in a finite field experiment
            that studies the statification of our parameter space. Here is a
            simple minded version of such an experiment:
        Example
            tally apply(100,i->bbC.degreeSingularLocusAt(random(K^1,K^20))) 
        Text
            We see that there is an open stratum of smooth cubics. The
            largest closed stratum consists of those cubics with a A1 singularity.
            The package finiteFieldExperiments helps to do the bookkeeping 
            for such experiments and also provides more detailed interpretation
            of the results.

 
///
end
---

--update := method();
--update(ExperimentData) := Tally => opts -> (experimentData) -> 
--( 
--   pointIterator := createIterator (  experiment.points() );
--   experimentData.pointData = new MutableHashTable;
--   experimentData.countData = new Tally;

--   while ( pointIterator.next() ) do
--   (
--       runExperimentOnce( experimentData, pointIterator.point(), pointsPerComponent );
--   );
--);


-- deprecated
estimateStratification2 = (e) -> 
(
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



estimateDecompositionOld := (experiment) -> 
(
    countData := experiment.counts();
    posRankJacobianAt := experiment.position(  experiment.usedRankJacobianAt() );
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
-- todo: updateExperiment?

    -- maybe think about writing 'connectProperty' ... (propName, bb.propName); default is 1:1.
    -- then, what should happen if a user requests to watch property xy ?
    

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

