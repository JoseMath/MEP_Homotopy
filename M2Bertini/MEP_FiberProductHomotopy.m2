
newPackage(
    "MEP_FiberProductHomotopy",
    Version => "1.0", 
    Date => "March 2019",
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "Jose@Math.wisc.edu",
       HomePage => "http://www.math.wisc.edu/~jose/"}
    },
    Headline => "Produces equations and solves multiparameter eigenvalue problems. ",
    DebuggingMode => true,
    AuxiliaryFiles => false,
    PackageImports => {"SimpleDoc","Bertini"},
    PackageExports => {"SimpleDoc","Bertini"},
  Configuration => { "RandomCoefficients"=>CC,
      "Continuation"=>Bertini },
  CacheExampleOutput => false
)


--path=prepend("/Users/jo/Documents/GoodGit/EuclideanDistanceDegree",path)
--loadPackage("EuclideanDistanceDegree",Reload=>true)

randomCC=()->random CC
randCC=()->random CC
randomRR=()->((-1)^(random(1,2)) *random RR)
randomZZ=()->random(1,30103)
randomValue=(kk)-> if kk===CC then randomCC() else if kk===RR then randomRR() else randomZZ() 
randomVector=method(Options=>{		})
randomVector(ZZ,Thing):= o->(n,R) ->apply(n,i->randomValue(R))--list of length n of randomValue

load"EDD_Determinantal.m2"

export { 
    "newMultiparameterEigenvalueProblem",
    "MultiparameterEigenvalueProblem"
            }

mepKeys={
    "PolynomialMatrices",
    "StartLinearSpace",
    "TargetLinearSpace",
    "ExtrinsicDimension",
    "IntrinsicDimension",
    "Directory",
    "StartSystemBertiniConfigurations",
    "TargetSystemBertiniConfigurations"    
    }


--###################################
-- TYPE DEFINITIONS
--###################################
MultiparameterEigenvalueProblem=new Type of MutableHashTable


--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
----------------------------------------
parString=(aString)->("("|toString(aString)|")");
addSlash=(aString)->(
    if aString_-1===" " then error (aString|" cannot end with whitespace.");
    if aString_-1=!="/" then aString=aString|"/";
    return aString    )
--newHyperplanes=A->for i to (numColumns A)+1 list randomVector(numRows A)
makeJac=(system,unknowns)->(--it is a list of lists of partial derivatives of a polynomial
         for i in system list for j in unknowns list  diff(j,i))


beginDocumentation()

load "./DOC_EDD.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end
  