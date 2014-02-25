-- does collected count really use the information about
-- number of components of the same dimension?

restart
needsPackage"BlackBoxIdeals"
needsPackage"FiniteFieldExperiments"

K = ZZ/5
R = K[x,y,z]

I = ideal (x*z*(z+x),y*z*(z+x))

-- make a black box from the ideal
bbI = blackBoxIdeal I;
e = new Experiment from bbI
e.run(1000)
e.collectedCount()
-- OK. Experiment collected twice the number of points in codim 1
-- since there are 2 components in codim 1.

assert(2.0 <= ((e.collectedCount())#{1}/e.minPointsPerComponent())*1.0)
assert(3.0 >= ((e.collectedCount())#{1}/e.minPointsPerComponent())*1.0)

