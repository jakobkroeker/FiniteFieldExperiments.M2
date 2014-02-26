uninstallPackage "BlackBoxIdeals"
installPackage "BlackBoxIdeals"
--loadPackage "BlackBoxIdeals"



   x  := null;
   x  = symbol x;
   rng := ZZ/7[x];
   coeffRng := coefficientRing rng;
   x = (gens rng)#0;

   RP := ZZ/7[x];
   IFP := ideal { 3*x^2+1, 5*x-1 };        
   IFPBlackBox := blackBoxIdeal( IFP );
   point := matrix {{3}};
   rng13 := ZZ/13;
   assert( IFPBlackBox.unknowns=={x} );
   assert( IFPBlackBox.equations==gens IFP);
   assert( IFPBlackBox.jacobian== jacobian IFP);
   

   point = matrix {{3}};

   point = sub( point, coeffRng ) ;
   IFPBlackBox.jacobianAt(point)
 
  rank    IFPBlackBox.jacobianAt(point)

IFPBlackBox.jetAt(point);
IFPBlackBox.jetAt(point, "liftDepth"=>10);

jetAt(IFPBlackBox,point,4,4);


