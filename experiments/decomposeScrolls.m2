-- decompose a union of scrolls using 
-- finite field experiments

quit -- F11 F11 F12

restart
path = append(path,"/Users/bothmer/Desktop/projekte/strudel/Jakob2010/GitHub/padicLiftM2/")

loadPackage "BlackBoxIdeals"
viewHelp BlackBoxIdeals
loadPackage "FiniteFieldExperiments"
viewHelp FiniteFieldExperiments


prime = 11
K =  ZZ/prime
n = 3
R = K[x_0..x_n]

randomScroll = (c) -> minors(2,random(R^{2:0},R^{c+1:-1}))
betti res randomScroll(3)

randomCI = (degs) -> ideal apply(degs,d->random(d,R))
betti res randomCI{2,5,5}

m = 3
betti (I = intersect apply(m,i->randomScroll(2)))
--time betti (I = intersect (apply(m,i->randomCI{2,2,5})))
--time dI = decompose I;
-- used   1.2584 seconds (m=2)
-- used  13.008 seconds (m=2)
-- used  10.454 seconds (m=2)
-- used 194.264 seconds (m=3)
-- used 420.376 seconds (m=3)
-- 17:44 - 21:54 nicht fertig (m=4)

bbI = blackBoxIdeal I;

-- experiment to find points on the components
e = new Experiment from bbI;

time e.run(1000)
-- used 0.584484 seconds (1000)
-- used 4.87354 seconds (2000, m=2, {2,5,5})
e.estimateDecomposition()

bbI.setSingularityTestOptions(5,1)

interpolProjection = createInterpolatedImage(e);
time interpolProjection.createAllInterpolatedIdeals(2,1)
-- used 6.11904 seconds (m=2, singularityPre=5)
-- used 16.5631 seconds (m=3, singularityPrc=5)
-- used 49.3681 seconds (m=4, singularityPrc=5)
-- used 119.282 seconds (m=5, singularityPrc=5)


time betti (Iproj = eliminate(I,{x_0,x_1}));
time apply(dIelim = decompose Iproj,degree)
betti time (I:saturate (I,dIelim#0))


-- used 31.2153 seconds (m=3, singularityPrc=10)



interpolProjection.interpolatedIdeals()

interpolProjection.bareIdeals()


P = (e.pointLists())#({2})

time DbbI = unique apply(1,i->time interpolateBB(2,bbI,P#i));#DbbI
-- used 92.6583 seconds (4)

-- used 12.5595 seconds (m=2,d=2)
-- 13:50
-- used 70.4629 seconds (m=2,d=3)

-- used 63.347 seconds (m=2,{2,2,5})

E100 = getEpsRing(K,100)
ideal E100

time jetAt(bbI,P#0,20,1);

