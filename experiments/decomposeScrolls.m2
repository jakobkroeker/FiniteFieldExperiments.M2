-- decompose a union of scrolls using 
-- finite field experiments

quit -- F11 F11 F12

restart
path = append(path,"/Users/bothmer/Desktop/projekte/strudel/Jakob2010/GitHub/padicLiftM2/")

loadPackage "BlackBoxIdeals"
loadPackage "FiniteFieldExperiments"

prime = 11
K =  ZZ/prime
n = 3
R = K[x_0..x_n]

randomScroll = (c) -> minors(2,random(R^{2:0},R^{n-c+2:-1}))

m = 3
betti (I = intersect apply(m,i->randomScroll(2)))
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
e.estimateDecomposition()

interpolProjection = createInterpolatedImage(e)
time interpolProjection.createAllInterpolatedIdeals(2,1)
-- used 41.6433 seconds (m=2)
-- used 114.558 seconds (m=3)

interpolProjection.printInterpolatedIdeals()
interpolProjection.ideal_9
interpolProjection.interpolatedIdealKeys("ideal_9")
interpolProjection



P = (e.pointLists())#({2})

time DbbI = unique apply(10,i->time interpolateBB(2,bbI,P#i));#DbbI
-- used 92.6583 seconds (4)

