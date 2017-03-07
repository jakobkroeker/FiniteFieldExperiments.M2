
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
  viewHelp FiniteFieldExperiments


	
K = ZZ/5
R = K[x,y,z]
I = ideal (x*z,y*z)
bb = blackBoxIdeal I;       
pointOnLine = matrix{{0,0,1_K}}
pointOnPlane = matrix{{0,1,0_K}}

bb.interpolateAt(pointOnLine, 1)
bb.interpolateAt(pointOnPlane, 1)
	   
   break
  restart 
 errorDepth =2
  loadPackage "BlackBoxIdeals"
--  check "BlackBoxIdeals"
--  uninstallPackage "BlackBoxIdeals"
--  installPackage "BlackBoxIdeals"
--  viewHelp BlackBoxIdeals
--  debug BlackBoxIdeals
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
  maxdeg = 1
  bb.interpolateComponentAt(p1,maxdeg)
-- wenn bei p1 schon mal interpoliert wurde mit kleinerem 
--   Grad, interpolieren mit maxdeg
-- wenn bei p1 schon mal interpoliert wurde mit groesser oder gleichem
--   grad, dann nichts machen
-- wenn bei p1 noch nicht interpoliert wurde, dann neu componente
--   interpolieren mit maxdeg
-- jet zum interpolieren immer laenge = anzahl der Monome+10
-- wenn kein jet dieser laenge gefunden werden kann, dann throw "singulaer"
-- und keine Interpolation

maxdeg = 2
bb.interpolateComponentsAt({p1,p2,p3},maxdeg)
-- ruft interpolateComponentAt(p_i,maxdeg) auf wenn:
--    sub(ideal componente, p_i) == 0  dann nicht interpolieren
--
--  ausser ideal componente == 0
--
--    nicht vorher testen ob glatt, sondern falls kein jet gefunden
--    werden kann => nicht interpolieren, keine Fehlermeldung
-- 
-- possible names:
--   quickIsOnComponentAt
--   isOnComponentAt(...,0)
p1
bb.interpolator.componentNamesAt(p1)
-- erzeugt jet j der laenge onComponenPrecision mit Anfang in point
--    wenn kein j dieser laege gefunden werden kann => "is singulaer"
-- sonst 
--     sub(ideal component, j) == 0
-- fuer alle componenten testen.

bb.interpolatedComponentsAt(p1)
--- {"certainlySingular",{}}
--- {"probablySmooth",{"c1","c2"}}

bb.interpolatedComponentNamesAt(singularPoint)
bb.interpolatedComponentNamesAt(p1)


-- bb.setOnComponentAnswerStrategy(SmoothnessInfoWithAnswerPair)
-- bb.setOnComponentAnswerStrategy(NullIfNotSmoothStrategy)


-- bb.onInterpolatedComponentPrecision()

-- bb.interpolator.onComponentPrecision()
-- bb.interpolator.setOnComponentPrecision(5)

  

-- bb.interpolator.setOnComponentPrecision(5) --- default =2
-- bb.interpolatedComponentNamesAt(p2)
-- bb.interpolatedComponents()

bb.onComponentAvailableAnswerStrategies()
bb.setOnComponentAnswerStrategy(PlainTextSmoothnessInfoWithAnswerPair)


bb.setOnComponentPrecision(5) --- default =2
bb.interpolatedComponentNamesAt(p3)


e.interpolateComponents(maxdeg)
bb.interpolateComponentsAt(e.points(),maxdeg)

e.watchProperty("interpolatedComponentNamesAt")
e.tryProperty("interpolatedComponentNamesAt")




  monomialDegree = 0;  
  pointList = {p1,singularPoint,p2,p3}
      bb.resetInterpolation()
  bb.interpolateComponentsAt(pointList,4)
  c3 = first bb.interpolator.componentsAt p3
  bb.isOnComponent(c3,p3)
  bb.isOnComponent("c3",p3)
  bb.isOnComponent("c1",origin)
  bb.isOnComponent("c1",p3)
  -- uses "isOnComponentPrecision"
  bb.interpolatedComponents()
  bb.interpolatedComponentNames()
  bb.interpolatedComponentByName "c1"
  bb.interpolateComponentAt p3  
  bb.interpolateComponentAt (p3,4)    
  bb.interpolateComponentAt p4
  bb.interpolatedComponents()
  bb.interpolatedComponentNames()
  bb.refineInterpolation()
  --pointList = {p1,p2,p3}
  --pointList = {p3}
  maxDegree=1
  iiList1 =  bb.interpolateComponentsAt(pointList,maxDegree)
  maxDegree = 3
  iiList3 = bb.interpolateComponentsAt(pointList,maxDegree) 
  bb.resetInterpolation()
  maxDegree = 4
  iiList4 = bb.interpolateComponentsAt(pointList,maxDegree)
  maxDegree = 5
  bb.interpolateComponentsAt(pointList,maxDegree)
  maxDegree = 3
   bb.interpolateComponentsAt(pointList,maxDegree) 
    bb.resetInterpolation()
  bb.interpolateComponentsAt(pointList)
  
  interpolator = bb.interpolator
  II = first  bb.interpolatedComponentsAt(p3)
  first II#"jetSet"#"jets"
  maximalConditions  II#"jetSet"
  II#"setName"("c5")
  II
  II#"name"() 
  interpolator.componentNames()
  interpolator.componentByName("c2")
  interpolator.renameComponent("c1","cc");
bb.interpolatedComponents()
bb.renameComponent("linear Factor 1","linearFactor1")
  interpolator.renameComponent("c2","ccc")
  interpolator.componentNamesAt(p1)
  interpolator.componentNamesAt(origin)
  interpolator.componentNamesAt(p3)  
  c3 =bb.interpolatedComponentByName("c5")
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
