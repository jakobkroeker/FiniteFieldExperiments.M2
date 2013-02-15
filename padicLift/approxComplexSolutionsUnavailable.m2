
approxComplexSolutions = method (Options=>{});
approxComplexSolutionsOld = method (Options=>{});
approxComplexSolutionsOld (Ideal, Matrix, List) := opts->(inputIdeal, solutionPoint, unknownList)->
(
   return null;
)


-- idealPointsApproxData: creates result data structure (a hash table)
--precondition: root coordinates corresponds to unknown list entries and unknownList is the same as gens ring equationsIdeal
idealPointsApproxData=( inputIdeal, solutionPoint, minPolyTable, approxSolutions, mergedLiftInfo, unknownList )->
(
   return null;
)


approxComplexSolutions (Ideal, Matrix) := opts->(inputIdeal, solutionPoint)->
(   
   return null;
)


doc ///
    Key
        approxComplexSolutions
        (approxComplexSolutions ,Ideal, Matrix)
    Headline
        not available - needs presence of pari/gp 
   Description
       Text
           This method is not available and currently requires presence of the pari/gp binary
///

-- todo: needs a test with more than one solution and one pairing - because I do not trust 'copy'
testApproxComplexSolutions = ()->
(
   return true;
)
