
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

randomCCMatrix=(m,n)->matrix apply(m,i->apply(n,j->random CC))

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
IntrinsicDimension="IntrinsicDimension";
ExtrinsicDimension="ExtrinsicDimension";
PolynomialMatrices="PolynomialMatrices";


--###################################
-- TYPE DEFINITIONS
--###################################
MultiparameterEigenvalueProblem=new Type of MutableHashTable
mepTypes:={"GenericExtrinsic",""}
TypeGenericExtrinsic="GenericExtrinsic";
Directory="Directory";
TypeEmpty="";
newMultiparameterEigenvalueProblem=method(Options=>{})
newMultiparameterEigenvalueProblem(List,List,String):= o->(n,d,type)->(
    if not member(type,mepTypes) then error("Last argument needs to be in "|mepTypes);
    dir:=temporaryFileName();
    if not fileExists dir then mkdir dir;
    mepH:=new MultiparameterEigenvalueProblem from {};
    scan(mepKeys,i->mepH#i=null);       
    mepH#Directory=dir;
    mepH#IntrinsicDimension=d;
    mepH#ExtrinsicDimension=n;    
    mepH#PolynomialMatrices=type;
    print n;
    if type==TypeGenericExtrinsic 
    then(
	scan(#n,i->mepH#("L"|i+1)=randomCCMatrix(#n-1,#n));
	mepH#("G")=randomCCMatrix((#n)^2-#n,(#n)^2-#n);
	scan(#n,i->mepH#("H"|i+1)=apply(#n+1,j->randomCCMatrix(n_i,n_i))));
    if type==TypeEmpty 
    then(
	scan(#n,i->mepH#("L"|i+1)=null);
	mepH#("G")=null;
	scan(#n,i->mepH#("H"|i+1)=null));
    return mepH)

writeMultiparameterEigenvalueProblem=method()
writeMultiparameterEigenvalueProblem(MultiparameterEigenvalueProblem):=(mepH)->(
    print mepH#Directory;
    if not (mepH#?"H1") then error" H1 is not a key. ";
    kk:=(first mepH#"H1")_(0,0)//ring;
    kParameters:=#mepH#ExtrinsicDimension;
    R:=kk[lam_0..lam_kParameters];
    print 1;
    scan(kParameters,i->(
	    xv:=apply(mepH#ExtrinsicDimension#i,j->"x"|j+1);
    	    print xv;
	    oneH:=0; 
	    scan(gens R,mepH#("H"|i+1),(l,A)->oneH=oneH+l*sub(A,R));
	    print (numrows oneH,numcols oneH,oneH);
	    bp1:=apply(entries oneH,r->makeB'Section(r,B'NumberCoefficients=>xv));
    	    bp0:=apply(entries mepH#("L"|i+1),r->makeB'Section(r,B'NumberCoefficients=>gens R));		
	    makeB'InputFile(mepH#Directory,NameB'InputFile=>"input_start"|i+1,
		B'Polynomials=>bp0|bp1,
		HomVariableGroup=>{xv},
		AffVariableGroup=>{drop(gens R,1)},
		B'Constants=>{toString first gens R=>1},
		B'Configs=>{"PrintPathProgress"=>100})		
		)));
runStartMultiparameterEigenvalueProblem=method()
runStartMultiparameterEigenvalueProblem(MultiparameterEigenvalueProblem,ZZ):=(mepH,s)->(
    runBertini(mepH#Directory,NameB'InputFile=>"input_start"|s);
    moveB'File(mepH#Directory,"nonsingular_solutions","start"|s))
runStartMultiparameterEigenvalueProblem(MultiparameterEigenvalueProblem):=(mepH)->scan(#mepH#ExtrinsicDimension,i->runStartMultiparameterEigenvalueProblem(mepH,i+1))

--writeStartFiberProductHomotopy







mepH=newMultiparameterEigenvalueProblem({4,5},{2,2},"GenericExtrinsic")    
writeMultiparameterEigenvalueProblem(mepH)
runStartMultiparameterEigenvalueProblem(mepH)

 
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

load "./DOC_MEP.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end
  