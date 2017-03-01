
         debug BlackBoxIdeals
    idealBlackBoxesProtect()

break
restart
loadPackage "FiniteFieldExperiments"
check "FiniteFieldExperiments"
uninstallPackage "FiniteFieldExperiments"
installPackage "FiniteFieldExperiments"

  break
  restart 
 errorDepth =2
  loadPackage "BlackBoxIdeals"
  check "BlackBoxIdeals"
  uninstallPackage "BlackBoxIdeals"
  installPackage "BlackBoxIdeals"
  viewHelp BlackBoxIdeals
  debug BlackBoxIdeals
  targetJetLength:= 10;
  


	
K = ZZ/5
R = K[x,y,z]
I = ideal (x*z,y*z)
bb = blackBoxIdeal I;       
pointOnLine = matrix{{0,0,1_K}}
pointOnPlane = matrix{{0,1,0_K}}

bb.interpolateAt(pointOnLine, 1,10)
bb.interpolateAt(pointOnPlane, 1,10)
	   
   break
  restart 
 errorDepth =2
  loadPackage "BlackBoxIdeals"
  check "BlackBoxIdeals"
  uninstallPackage "BlackBoxIdeals"
  installPackage "BlackBoxIdeals"
  viewHelp BlackBoxIdeals
  debug BlackBoxIdeals
  -- to get meaningful debug output for packages,
  -- we need errorDepth = 2. (default value is 3) => Lower value means better 
  errorDepth = 2
  --kk := ZZ  
  --0_kk - fails
  --while 
  kk= ZZ
  0_kk --is fine
  R := kk[x,y]
  I= ideal (x*y*(x^2-y^2)*(y^4-3*x-7))
  bb = new BlackBoxIdeal from I
  p1 = matrix {{1,1_kk}}
  p2 = matrix {{1,0_kk}}
  p3 = matrix {{3,2_kk}}
  p4 =  matrix {{2,2_kk}}
  origin = matrix {{0, 0_kk}}
  singularPoint  = matrix{{0,0_kk}}
  bb.valuesAt p1
  bb.valuesAt p2
  bb.valuesAt p3
  monomialDegree = 0;  
  pointList = {p1,singularPoint,p2,p3}
      bb.resetInterpolation()
  bb.interpolateComponents(pointList)
  c3 = first bb.componentsAt p3
  bb.isOnComponent(c3,p3)
  bb.isOnComponent("c3",p3)
  bb.isOnComponent("c1",origin)
  bb.isOnComponent("c1",p3)
  bb.components()
  bb.componentNames()
  bb.componentByName "c1"
  bb.interpolateAt p3  
  bb.interpolateAt (p3,4)    
  bb.interpolateAt p4
  bb.components()
  bb.componentNames()
  bb.refineInterpolation()
  --pointList = {p1,p2,p3}
  --pointList = {p3}
  maxDegree=1
  iiList1 =  bb.interpolateComponents(pointList,maxDegree)
  maxDegree = 3
  iiList3 = bb.interpolateComponents(pointList,maxDegree) 
  bb.resetInterpolation()
  maxDegree = 4
  iiList4 = bb.interpolateComponents(pointList,maxDegree)
  maxDegree = 5
  bb.interpolateComponents(pointList,maxDegree)
  maxDegree = 3
   bb.interpolateComponents(pointList,maxDegree) 
    bb.resetInterpolation()
  bb.interpolateComponents(pointList)
  
  interpolator = bb.interpolator
  II = first  bb.componentsAt(p3)
  first II#"jetSet"#"jets"
  maximalConditions  II#"jetSet"
  II#"setName"("c5")
  II
  II#"name"() 
  interpolator.componentNames()
  interpolator.componentByName("c2")
  interpolator.renameComponent("c1","cc");
  interpolator.renameComponent("c2","ccc")
  interpolator.componentNamesAt(p1)
  interpolator.componentNamesAt(origin)
  interpolator.componentNamesAt(p3)  
  c3 = bb.componentByName("c3")
  peek c3#"jetSet"
  jet = first  c3#"jetSet"#"jets"
  length jet
  bb.valuesAt   jet#"value"
  sort  bb.components()
-- does not work: 
--   bb.valuesAt   jet
II  
peek II  
II#"name"()
II#"setName"("bububu")
II
size II#"jetSet"
