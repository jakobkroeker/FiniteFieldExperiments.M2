newPackage(
     "IntervalPkg",
     Version => "0.1", 
     Date => "02.08.2013",
     Authors => {
           { Name => "Jakob Kroeker", 
           Email => "kroeker@uni-math.gwdg.de", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"}
      },
     Configuration => {},
     PackageExports => {},
     Headline => " a rudimentary implementation of an Interval and corresponding operations.",
     DebuggingMode => true
)


export {
  "Interval",
  "IntervalWithCenter"
}


IntervalProtect = ()->
(
  protect center;
);

IntervalExport =()->
(
 exportMutable("center");
);

IntervalExport();


Interval = new Type of HashTable;

roundInterval = method();

-- todo: use interval arithmetic
roundInterval ( HashTable, ZZ) := Interval=>(I,n)->
(
     return new Interval from ( round(n,I.min), round(n,I.max) );
)




new Interval from Thing := (E,ll) -> 
(
   error ("this interface is not implemented");
);

--new Interval from HashTable := (E,ll) -> 
new Interval from List := (E,ll) -> 
( 
  ht := new MutableHashTable from ll;
  --print "new Interval from HashTable";
  assert(ht#?(symbol min) );
  assert(ht#?(symbol max) );
  assert(ht.min<=ht.max);
 
  ht.round = method();
  ht.round(ZZ) := Interval => (n)->
  (
     return roundInterval( ht, n);
  );

  --print( keys ht);
  return new HashTable from ht;
);

-- maybe Interval coould also contain statistical information?

new Interval from Sequence := (E,seq) -> 
(
    --print "new Interval from Sequence";
    assert(#seq==2);
    ht := { (symbol min)=>seq#0*1.0, (symbol max)=>seq#1*1.0};
    return new Interval from ht;
);







formatIntervalBounds = (num)->
(
  f :=(pnum) ->format(10,2,2,2,pnum);
  str:="";
  if (value f num)==0 then return "0.";
  if (abs(value f num)<1) then 
  str = str|"0";
  origstr := format(10,2,2,2,num);
  if origstr#0=="-" then
  (
    origstr = origstr_(1..#origstr-1);
     str = "-"|str;
  );
  str = str|origstr;
  return str;
)


net (Interval) := Net =>(interval)->
(
   str :=  " " | "[" | formatIntervalBounds  interval.min | ", " | formatIntervalBounds interval.max | "]";
   return net str;
)


RR * Interval := (scalingFactor,I) -> (
     new Interval from (scalingFactor*I.min,scalingFactor*I.max)
     )

Interval * RR := (I,scalingFactor) -> (
     scalingFactor*I
     )

ZZ * Interval := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

Interval * ZZ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )

QQ * Interval := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

Interval * QQ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )



-- problem: Sortierung nicht ohne weiteres austauschbar -> Wrappertyp müsste eingeführt werden, der eine andere Sortierung implementiert
-- und alle Objekte müssen in den Wrappertyp konvertiert werden.


Interval ? Interval := (I1,I2) ->
(
   if  ( I1.center == I2.center ) then 
   (
        I2.min ? I1.min 
   )
   else 
   (  
         I1.center ? I2.center
   )
)

Interval == Interval := (I1,I2) -> (
     (I1.min == I2.min) and (I1.max == I2.max)
     )


IntervalWithCenter = new Type of Interval;


intervalCenter = method();


intervalCenter (HashTable) := Interval=>(I)->
(
    return (I.min+I.max)/2.0;
)

new IntervalWithCenter from Thing := (E,ll) -> 
(
   error ("this interface is not implemented");
);

roundIntervalWithCenter = method();

-- todo: use interval arithmetic
roundIntervalWithCenter ( HashTable, ZZ) := IntervalWithCenter=>(I,n)->
(
     return new IntervalWithCenter from ( round(n,I.min), round(n,I.max) );
)




new IntervalWithCenter from List := (E,ll) -> 
( 
  ht := new MutableHashTable from ll;
  --print "new Interval from HashTable";
  assert(ht#?(symbol min) );
  assert(ht#?(symbol max) );
  assert(ht.min<=ht.max);
 
  ht.round = method();
  ht.round(ZZ) := IntervalWithCenter => (n)->
  (
     return roundIntervalWithCenter( ht, n);
  );

  ht.center = intervalCenter(ht);

  --print( keys ht);
  return new HashTable from ht;
);


new IntervalWithCenter from Sequence := (E,seq) -> 
(
    --print "new Interval from Sequence";
    assert(#seq==2);
    ht := { (symbol min)=>seq#0*1.0, (symbol max)=>seq#1*1.0};
    return new IntervalWithCenter from ht;
);



net (IntervalWithCenter) := Net =>(interval)->
(
   str :=  " " | formatIntervalBounds  interval.center  | " " | "[" | formatIntervalBounds  interval.min | ", " | formatIntervalBounds interval.max | "]";
   return net str;
)


RR * IntervalWithCenter := (scalingFactor,I) -> (
     new IntervalWithCenter from (scalingFactor*I.min,scalingFactor*I.max)
     )

IntervalWithCenter * RR := (I,scalingFactor) -> (
     scalingFactor*I
     )

ZZ * IntervalWithCenter := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

IntervalWithCenter * ZZ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )

QQ * IntervalWithCenter := (scalingFactor,I) -> (
     (scalingFactor*1.0)*I
     )

IntervalWithCenter * QQ := (I,scalingFactor) -> (
     (scalingFactor*1.0)*I
     )


-- provide an ordering of intervals. 
-- If they have the same center, then the left bound is compared.

IntervalWithCenter ? IntervalWithCenter := (I1,I2) ->
(
   if  ( I1.center == I2.center ) then 
   (
        I2.min ? I1.min 
   )
   else 
   (  
         I1.center ? I2.center
   )
)






TEST ///
   debug IntervalPkg
   IntervalProtect()
   assert (2*(new IntervalWithCenter from (1,2)) == new IntervalWithCenter from (2,4))
   assert (2.0*(new IntervalWithCenter from (1,2)) == new IntervalWithCenter from (2,4))
///

IntervalWithCenter == IntervalWithCenter := (I1,I2) -> (
     (I1.min == I2.min) and (I1.max == I2.max)
     )

TEST ///
  debug IntervalPkg
  IntervalProtect()
  new IntervalWithCenter from { (symbol min)=>0, (symbol max)=>1};
  i2 := new IntervalWithCenter from (1, 1.04);
  i3 := new IntervalWithCenter from (1, 1.045);
  assert( i3.round(2)== i2);
  i4 :=  new IntervalWithCenter from (2.5,3.5);
  assert(i4.center == 3.0);
  i2<i3;
  i2?i3;
  assert ( (i2?i3) === (1?2));
  assert ( (i2?i2) === (1?1));
  i2?i2;
///




