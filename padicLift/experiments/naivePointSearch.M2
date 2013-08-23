
load "MSchuett.m2";


----------------------- naive Punktsuche --------------------------------------------
----- nur parallelisiert schnell genug .



nextPoint = point -> ( 
 added := false; 
 pos := 0;
 while not added do 
 (  
   point_(0,pos) = point_(0,pos)+1 ;

   if ( point_(0,pos)!=0 or ( numColumns point<=pos+1) ) then
     (    added = true;   )  
   else (     pos=pos+1;   ); 
);

  if (point==0) then return false;  
  return point; 
);


nextPointEx = (point,num) ->
(
	while num>0 do
	(
		point = nextPoint (point);
		if point===false then return point;
		num=num-1;
	);
	return point;
)


varnum := #(gens rng);
solutions := {};

point :=  apply(varnum, i->0) ;
point = apply(point, i->sub( i, ring idFp));
point = mutableMatrix { point };

sub( idFp, matrix point )


 progress := 0;

while ( nextPoint(point) =!= false ) do (  
   progress = progress+1;   print(progress);
   if substitute(id,matrix point)==0 then    solutions =solutions | { point }; 

);




