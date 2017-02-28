
         debug BlackBoxIdeals
    idealBlackBoxesProtect()

  break
  restart 
 errorDepth =2
  loadPackage "BlackBoxIdeals"
  check "BlackBoxIdeals"
  uninstallPackage "BlackBoxIdeals"
  installPackage "BlackBoxIdeals"
  debug BlackBoxIdeals

K = ZZ/5
R = K[x,y,z]
I = ideal (x*z,y*z)
bb = blackBoxIdeal I;       
pointOnLine = matrix{{0,0,1_K}}
pointOnPlane = matrix{{0,1,0_K}}

bb.interpolateAt(pointOnLine, 1)
bb.interpolateAt(pointOnPlane, 1)
	   
  
  -- to get meaningful debug output for packages,
  -- we need errorDepth = 2. (default value is 3) => Lower value means better 
  errorDepth = 2
  --kk := ZZ  
  --0_kk - fails
  --while 
  kk=ZZ
  0_kk --is fine
  R := kk[x,y]
  I:= ideal (x*y*(x^2-y^2)*(y^4-3*x-7))
  bb = new BlackBoxIdeal from I
  p1 = matrix {{1,1_kk}}
  p2 = matrix {{1,0_kk}}
  p3 = matrix {{3,2_kk}}
  origin = matrix {{0, 0_kk}}
  singularPoint  = matrix{{0,0_kk}}
  bb.valuesAt p1
  bb.valuesAt p2
  bb.valuesAt p3
  pointList = {p1,singularPoint,p2,p3}
  maxDegree=1
  onComponentPrecision = 1
  iiList = bb.interpolateComponents(pointList,maxDegree,onComponentPrecision)
  maxDegree = 3
  iiList2 = bb.interpolateComponents(pointList,maxDegree,onComponentPrecision)
  maxDegree = 4
  iiList2 = bb.interpolateComponents(pointList,maxDegree,onComponentPrecision)
  ip = bb.interpolator
  II = first  iiList
  II#"jetSet"
  --II#"setName"("c_2")
  II#"setName"("c_5")
  ip.componentByName("c_9")
  ip.renameComponent("c_1","cc");
  ip.renameComponent("c_9","ccc")
  ip.componentNamesAt(p1)
  ip.componentNamesAt(origin)
  
sr = select(  iiList, (val)->("c_5"===val#"name"()))
#sr
first values sr
iiList
  
II  
peek II  
II#"getName"()
netInterpolatedIdeal(II)
size II#"jetset"
l0 = toString "\"bl1222 "  
l1 = {1,223,3}
ls1 = apply(l1, i->toString i)
ls2 = apply(2, i->"=>")


sls1 = stack ls1
sls2 = stack ls2

l0 |sls1 |sls2
