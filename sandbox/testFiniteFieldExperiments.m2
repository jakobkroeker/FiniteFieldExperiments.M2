restart

uninstallPackage "FiniteFieldExperiments"
installPackage "FiniteFieldExperiments"
loadPackage("FiniteFieldExperiments",Reload=>true)

check FiniteFieldExperiments;

 

 debug FiniteFieldExperiments
    FiniteFieldExperimentsProtect()
    coeffRing := ZZ/3;
    bbRankM = blackBoxIdeal( 5 ,coeffRing )
    rankMat := (point)->5
    bbRankM.registerPointProperty("rankMat",rankMat)


    point := matrix {{1,2,3,4,5}};
    point = sub( point, coeffRing);

    bbRankM = bbRankM.rebuild()

    e = new RandomExperiment from bbRankM
    keys e
    assert (e.coefficientRing()===coeffRing);

    e.setPointsPerComponent(20);
    assert( e.pointsPerComponent()==20);
    FFELogger.setLevel(4);
    e.watchProperties {"rankMat"};
     e.watchedProperties()
    assert( 1== # select( e.watchedProperties(), 
                       (prop)->(prop=="rankMat") ) 
     )


    e.useRankJacobianAt("rankMat");
    e.useRankJacobianAt(null);
    e.pointKeys()
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


    e2 = new Experiment from bbRankM
    e2.setPointsPerComponent(20);
    e2.watchProperties {"rankMat"};
    e2.run(1000);
    e.merge(e2);



