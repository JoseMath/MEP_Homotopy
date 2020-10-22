
newPackage(
    "MultiparameterEigenvalueProblemHomotopy",
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
--loadPackage("MultiparameterEigenvalueProblemHomotopy",Reload=>true)

randomCC=()->random CC
randCC=()->random CC
randomRR=()->((-1)^(random(1,2)) *random RR)
randomZZ=()->random(1,30103)
randomValue=(kk)-> if kk===CC then randomCC() else if kk===RR then randomRR() else randomZZ() 
randomVector=method(Options=>{		})
randomVector(ZZ,Thing):= o->(n,R) ->apply(n,i->randomValue(R))--list of length n of randomValue

randomCCMatrix=(m,n)->matrix apply(m,i->apply(n,j->random CC))
randomCCRankOneMatrix=(m,n)->randomCCMatrix(m,1)*randomCCMatrix(1,n)
randomCCRankMatrix=(k,m,n)->(
    M:=randomCCRankOneMatrix(m,n);scan(k-1,i->M=M+randomCCRankOneMatrix(m,n)); return M )

--load"EDD_Determinantal.m2"

export { 
    "newMultiparameterEigenvalueProblem",
    "MultiparameterEigenvalueProblem",
    "writeMultiparameterEigenvalueProblem",
    "runStartMultiparameterEigenvalueProblem",
    "getStartLines",
    "writeStartSolutionsFiberProductHomotopy",
    "mepFiberProductHomotopy",
    "diagonalCoefficientHomotopy"
            }

mepKeys={
    "PolynomialMatrices",
    "StartLinearSpace",
    "TargetLinearSpace",
    "ExtrinsicDimension",
    "IntrinsicDimension",
    "Directory",
    "StartSystemBertiniConfigurations",
    "TargetSystemBertiniConfigurations",
    "TypeStartSolution",    
    "Homogenize"
    }
IntrinsicDimension="IntrinsicDimension";
ExtrinsicDimension="ExtrinsicDimension";
PolynomialMatrices="PolynomialMatrices";
TypeStartSolution="TypeStartSolution";
Homogenize="Homogenize";
startSolutionTypes={"nonsingular_solutions"}
NonsingularSolutions="nonsingular_solutions"

--###################################
-- TYPE DEFINITIONS
--##################################
MultiparameterEigenvalueProblem=new Type of MutableHashTable
StartSolutionsFiberProductHomotopy=new Type of MutableHashTable
mepTypes:={"GenericExtrinsic","","DimensionDeficientRankOne","DimensionDeficientRankNumberOfParameters","DimensionDegree"}
TypeGenericExtrinsic="GenericExtrinsic";
DimensionDeficientRankOne="DimensionDeficientRankOne"
DimensionDeficientRankNumberOfParameters="DimensionDeficientRankNumberOfParameters"
DimensionDegree="DimensionDegree"
Directory="Directory";
TypeEmpty="";
newMultiparameterEigenvalueProblem=method(Options=>{})
newMultiparameterEigenvalueProblem(List,String):= o->(n,type)->(
    if not member(type,mepTypes) then error("Last argument needs to be in "|mepTypes);
    dir:=temporaryFileName();
    if not fileExists dir then mkdir dir;
    mepH:=new MultiparameterEigenvalueProblem from {};
    scan(mepKeys,i->mepH#i=null);       
    mepH#Directory=dir;
    mepH#IntrinsicDimension=n;
    mepH#ExtrinsicDimension=n;    
    mepH#PolynomialMatrices=type;
    mepH#TypeStartSolution=NonsingularSolutions;
    mepH#Homogenize=false;
    mepH#"TargetSystemBertiniConfigurations"={};
    mepH#"StartSystemBertiniConfigurations"={};
--    print n;
    if type==TypeGenericExtrinsic 
    then(
	scan(#n,i->mepH#("L"|i+1)=randomCCMatrix(#n-1,#n+1));
	mepH#("G")=randomCCMatrix((#n)^2-#n,(#n)^2-#n);
	scan(#n,i->mepH#("H"|i+1)=apply(#n+1,j->randomCCMatrix(n_i,n_i))));
    if type==DimensionDeficientRankOne 
    then(
	scan(#n,i->mepH#("L"|i+1)=randomCCMatrix(#n-1,#n+1));
	mepH#("G")=randomCCMatrix((#n)^2-#n,(#n)^2-#n);
	scan(#n,i->mepH#("H"|i+1)=apply(#n+1,j->if j=!=0 then randomCCRankOneMatrix(n_i,n_i) else randomCCMatrix(n_i,n_i) ))
	);
    if type==DimensionDeficientRankNumberOfParameters 
    then(
	scan(#n,i->mepH#("L"|i+1)=randomCCMatrix(#n-1,#n+1));
	mepH#("G")=randomCCMatrix((#n)^2-#n,(#n)^2-#n);
	print DimensionDeficientRankNumberOfParameters;
	scan(#n,i->mepH#("H"|i+1)=apply(#n+1,j->if j=!=0 then randomCCRankMatrix(#n,n_i,n_i) else randomCCMatrix(n_i,n_i) ))
	);
    if type==DimensionDegree 
    then(
	scan(#n,i->mepH#("L"|i+1)=randomCCMatrix(#n-1,#n+1));
	mepH#("G")=randomCCMatrix((#n)^2-#n,(#n)^2-#n);
    	mepH#IntrinsicDimension=apply(#n,i->#n);--(k,k,...,k)
	print DimensionDegree;
--    	if not mepH#?"IntrinsicDimension" then error"IntrinsicDimension not set. ";
	scan(#n,i->(
	    rankOneMatrices:=apply((mepH#IntrinsicDimension)_i,i->randomCCRankOneMatrix(n_i,n_i));
    	    getMatrixInPencil:=()->(M:=0;scan(rankOneMatrices,m->M=m+M);return M);
	    mepH#("H"|i+1)=apply(#n+1,j->if j=!=0 then getMatrixInPencil() else randomCCMatrix(n_i,n_i) )
		)));
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
    R:=kk[apply(kParameters+1,i->"lam"|i)];
--    print 1;
    scan(kParameters,i->(
	    xv:=apply(mepH#ExtrinsicDimension#i,j->"x"|j+1);
    	    --print xv;
	    oneH:=0; 
	    scan(gens R,mepH#("H"|i+1),(l,A)->oneH=oneH+l*sub(A,R));
	    --print (numrows oneH,numcols oneH,oneH);
	    bp1:=apply(entries oneH,r->makeB'Section(r,B'NumberCoefficients=>xv));
    	    bp0:=apply(entries mepH#("L"|i+1),r->makeB'Section(r,B'NumberCoefficients=>gens R));		
	    makeB'InputFile(mepH#Directory,NameB'InputFile=>"input_start"|i+1,
		B'Polynomials=>bp0|bp1,
		HomVariableGroup=>{xv},
		AffVariableGroup=>{drop(gens R,1)},
		B'Constants=>{toString first gens R=>1},
		BertiniInputConfiguration=>{"PrintPathProgress"=>100})		
		));
    ---Now we write the Fiberproduct homotopy part. 	    
    allLam:=apply(kParameters,i->apply(kParameters+1,j->"lam"|i+1|"v"|j));
    S:=kk[flatten allLam];
    lvFix:=(allLam_0)/(i->sub(value i,S));
    fp1:={};    
    fpL:={};
    fpG:={};
    X:={};
    countL:=0;
    countG:=0;    
    --print("Begin scan");
    scan(kParameters,i->(
	xv:=apply(mepH#ExtrinsicDimension#i,j->"x"|i|"v"|j+1);
    	X=append(X,xv);
	lv:=(allLam_i)/(i->sub(value i,S));    	
	--print xv;
	--print lv;
	oneH:=0; 
	scan(lv,mepH#("H"|i+1),(l,A)->oneH=oneH+l*sub(A,S));
	--print (numrows oneH,numcols oneH,oneH);
	fp1=fp1|apply(entries oneH,r->makeB'Section(r,B'NumberCoefficients=>xv));
	scan(entries mepH#("L"|i+1),r->(
	    	countL=countL+1;		
--		print ("countL"=>countL);
		fpL=append(fpL,
		    makeB'Section(r,B'NumberCoefficients=>lv,
		    NameB'Section=>"linearL"|countL))));
    	print("number fpL"=>#fpL)));
    allLam=apply(kParameters,i->apply(kParameters+1,j->sub(value("lam"|i+1|"v"|j),S)));
    baseG:=flatten apply(#allLam-1,i->drop(allLam_0-allLam_(i+1),1));
--    print ("baseG",#baseG,baseG);
    scan(kParameters^2-kParameters,i->(
		countG=countG+1;
		--print ("countG"=>countG);
		gc:=((entries(mepH#"G"))_(countG-1));
		--print (#gc);
		fpG=append(fpG,makeB'Section(baseG,
			B'NumberCoefficients=>gc,
			NameB'Section=>"linearG"|countG))));
--    print("End scan",#fpL==#fpG);    
    fp0:=apply(#fpL,i->"(1-fpT)*linearL"|i+1|"+(fpT)*linearG"|i+1);
--    print fp0;
--    print("write input_FPH");
    writeParameterFile(mepH#Directory,{0},NameParameterFile=>"start_parameters");
    writeParameterFile(mepH#Directory,{1},NameParameterFile=>"final_parameters");
--    print("fpL",fpL);
    bcFP:=mepH#"TargetSystemBertiniConfigurations"|{"ParameterHomotopy"=>2,"PrintPathProgress"=>100};
    if not mepH#Homogenize 
    then (avg,hvg,bConstant,bFunction):=(allLam/(i->drop(i,1))//flatten,X,(allLam/first/toString/(i->i=>1)),fpL|fpG) 
    else  (avg,hvg,bConstant,bFunction)=({},append(X,flatten prepend("lamZero",allLam/(i->drop(i,1)))),{},(allLam/first/toString/(i->i=>"lamZero"))|fpL|fpG);
    makeB'InputFile(mepH#Directory,NameB'InputFile=>"input_FPH",
		B'Polynomials=>fp0|fp1,
		HomVariableGroup=>hvg,
		AffVariableGroup=>avg,
    	    	ParameterGroup=>{"fpT"},
		B'Constants=>bConstant,
    	    	B'Functions=>bFunction,
		BertiniInputConfiguration=>bcFP));

	
runStartMultiparameterEigenvalueProblem=method()
runStartMultiparameterEigenvalueProblem(MultiparameterEigenvalueProblem,ZZ):=(mepH,s)->(
    runBertini(mepH#Directory,NameB'InputFile=>"input_start"|s);
    moveB'File(mepH#Directory,NonsingularSolutions,"start"|s))
runStartMultiparameterEigenvalueProblem(MultiparameterEigenvalueProblem):=(mepH)->scan(#mepH#ExtrinsicDimension,i->runStartMultiparameterEigenvalueProblem(mepH,i+1))

getStartLines=method()
getStartLines(MultiparameterEigenvalueProblem,ZZ):=(mepH,s)->(
    theLines:=lines get(addSlash(mepH#Directory)|"start"|s);
    ns:=value first theLines;
    theLines=drop(theLines,1);  
--    print (theLines_0,theLines_1)  ;
    count:=0;
    allSolutions:={};
    pDim:=#mepH#ExtrinsicDimension;
    iDim:=mepH#IntrinsicDimension#(s-1);    
    eDim:=mepH#ExtrinsicDimension#(s-1);
--    print(pDim,iDim,eDim);
    while count<ns do(
--	print ("count"=>count);
	x:=(apply(eDim,
		i->theLines_(i+1)));
--	print (#x);
    	theLines=drop(theLines,1+eDim);
	lam:=apply(pDim,i->theLines_i);
--    	print (#lam);
    	theLines=drop(theLines,pDim);
	allSolutions=append(allSolutions,(x=>lam));
--	print count;
	count=count+1);
    return (ns,allSolutions)			
	)

writeStartSolutionsFiberProductHomotopy=method()
writeStartSolutionsFiberProductHomotopy(MultiparameterEigenvalueProblem):=(mepH)->(
    if not (mepH#"Homogenize") then ss:={{}=>{}} else ss={{}=>{"1e0 0e0"}};
    scan(#mepH#ExtrinsicDimension,s->(
--	print "getStartLines";    
	(ns,allSolutions):=getStartLines(mepH,s+1);
--	print first allSolutions;
--	print ss;
	ss=flatten apply(ss,oldSol->apply(allSolutions,newSol->(
		x:=((first oldSol)|(first newSol));
		lam:=((last oldSol)|(last newSol));
		x=>lam)))));
    SFile:= openOut(addSlash(mepH#Directory)|"start_FPH"); 
    SFile<< #ss<<endl;
    scan(ss,oneSol->(
	    SFile<<endl;
	    scan(first oneSol,x->SFile<<x<<endl);
	    scan(last oneSol,lam->SFile<<lam<<endl);
	    ));
    close SFile;
    return ss)
---

mepFiberProductHomotopy=method();
mepFiberProductHomotopy(MultiparameterEigenvalueProblem):=mepH->(
    writeMultiparameterEigenvalueProblem(mepH);
    runStartMultiparameterEigenvalueProblem(mepH);
    writeStartSolutionsFiberProductHomotopy(mepH);
    moveB'File(mepH#Directory,"start_FPH","start",CopyB'File=>true);
    runBertini(mepH#Directory,NameB'InputFile=>"input_FPH");
    readFile(mepH#Directory)    
    )



-*
PolynomialMatrix=new Type of MutableHashTable 
new PolynomialMatrix := (X)->new MutableHashTable from {
    "MatrixCoefficients"=>null,
    "SupportPolynomial"=>null,
    "ExtrinsicDimension"=>null,
    "IntrinsicDimension"=>null}
*-

diagonalCoefficientHomotopy = method(Options=>{})
diagonalCoefficientHomotopy(List,List) := o -> (allD,allH) -> (
    extDim := allH/first/numrows;
    k := #allH;
    lam := symbol lam;    
    dhT := symbol dhT;    
    R := (ring allH#0#0)[lam_1..lam_k,dhT];
    allD1:=allD/entries/first/diagonalMatrix/(i->sub(i,R));
    allD2:=allD/entries/last/diagonalMatrix/(i->sub(i,R));
    mons:={1}|drop(gens R,-1);
    startMatrices := apply(#allH, 
	h -> (1-dhT)*(allD1#h) + (1-dhT)*lam_(h+1)*allD2#h);
    targetMatrices := apply(#allH, 
	h -> dhT*sum(allH#h,mons,(a,b)->b*sub(a,R)));
    sys:= startMatrices+targetMatrices;
    xVector :=(i,j)->apply(extDim#i,a->"x"|i|"v"|a);
    win:=apply(#sys,m->apply(entries sys#m, r-> makeB'Section(r,B'NumberCoefficients=>xVector(m,#r))));
--    win = flatten win/(i->i#B'SectionString);
    dir :=temporaryFileName();
    mkdir dir;
    makeB'InputFile(dir,
	AffVariableGroup=>toList (lam_1..lam_k),
	HomVariableGroup=>apply(#extDim,i->xVector(i,extDim#i)),
    	ParameterGroup=>{last gens R},
	B'Polynomials=>flatten win,
	BertiniInputConfiguration=>{"ParameterHomotopy"=>2}
	);
    S := apply(allD,D1->set apply(numcols D1,j->{
		j => -((entries(D1_j))#0)/((entries(D1_j))#1)
		}
	    )
	);
    cartesianProduct := (A,B)->  (toList (A**B))/toList/flatten//set;
    S = toList fold(cartesianProduct,S);
    unitVector := (a,b) -> apply(b,i->if i == a then 1 else 0);
    makeSolution := s->(
    	p:={};
    	scan(s/first,extDim,(a,b)->p=p|unitVector(a,b));
    	p=p|toList(s/last);
    	p);
    S=S/makeSolution;
    writeParameterFile(dir,{0},NameParameterFile=>"start_parameters");
    writeParameterFile(dir,{1},NameParameterFile=>"final_parameters");
    writeStartFile(dir,S);
    runBertini(dir);
    readFile(dir)    
    )

-*
restart
path=prepend("/Users/jo/Documents/GoodGit/MEP_Homotopy/Bertini/M2Bertini",path)
needsPackage"MultiparameterEigenvalueProblemHomotopy"

randMatrix = (i,j)-> matrix for i to i-1 list for j to j-1 list random CC
--Start system
D1 = randMatrix(2,3)
D2 = randMatrix(2,4)

--Target system
H1 = apply(2+1,i->randMatrix(3,3))
H2 = apply(2+1,i->randMatrix(4,4))
allD={D1,D2}--specify a start system by taking the first row to be a diagonal matrix of constants and the second row is diagonal matrix multiplied by a lambda. 
allH={H1,H2}

--Solve
diagonalCoefficientHomotopy(allD,{H1,H2})
*-

-*
Experiment setup
restart
path=prepend("/Users/jo/Documents/GoodGit/MEP_Homotopy/Bertini/M2Bertini",path)
needsPackage"MultiparameterEigenvalueProblemHomotopy"

extDim={5,5,5}
randMatrix = (i,j)-> matrix for i to i-1 list for j to j-1 list random CC
--Start system
allD=apply(#extDim,i->randMatrix(2,extDim#i))--specify a start system by taking the first row to be a diagonal matrix of constants and the second row is diagonal matrix multiplied by a lambda. 

--Target system
allH = apply(extDim,n->apply(#extDim+1,i->randMatrix(n,n)))
#allH
--Solve
diagonalCoefficientHomotopy(allD,allH)

*-

--PolynomialMatrix + PolynomialMatrix := (x,y)->




 
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

--load "./DOC_MEP.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end
  
  
  
  restart
path=prepend("/Users/jo/Documents/GoodGit/MEP_Homotopy/Bertini/M2Bertini",path)
needsPackage"MultiparameterEigenvalueProblemHomotopy"

loadPackage("MultiparameterEigenvalueProblemHomotopy",Reload=>true)
printingPrecision=100
--Simple case.
mepH=newMultiparameterEigenvalueProblem({7,6},"GenericExtrinsic")    
peek mepH
theDir=temporaryFileName()
mkdir theDir
mepH#"Directory"=theDir
sols = mepFiberProductHomotopy(mepH)


"/Users/jo/Desktop/Dump"
"nonsingular_solutions"
--Repeated eigenvalue case. 
mepH=newMultiparameterEigenvalueProblem({3,3},"GenericExtrinsic")    
mepH#"H1"
A10=diagonalMatrix apply(3,i->random CC)
alpha=random CC
A10=diagonalMatrix {alpha,alpha,2_CC}
A11=diagonalMatrix{1,1_CC,2}
A12=diagonalMatrix{1,1_CC,2}
mepH#"H1"={A10,A11,A12}
mepFiberProductHomotopy(mepH)
----
mepH=newMultiparameterEigenvalueProblem({3,3,3,3},"GenericExtrinsic")    
writeMultiparameterEigenvalueProblem(mepH)
runStartMultiparameterEigenvalueProblem(mepH)
getStartLines(mepH,2)
writeStartSolutionsFiberProductHomotopy(mepH);
moveB'File(mepH#"Directory","start_FPH","start",CopyB'File=>true);
runBertini(mepH#"Directory",NameB'InputFile=>"input_FPH");
readFile(mepH#"Directory")
--------

mepH=newMultiparameterEigenvalueProblem({4,5,3},"DimensionDeficientRankOne")    ;
mepFiberProductHomotopy(mepH)
readFile(mepH#"Directory","start1",1000)
readFile(mepH#"Directory","start2",1000)
readFile(mepH#"Directory","start3",1000)
readFile(mepH#"Directory","nonsingular_solutions",1000)
readFile(mepH#"Directory")--2+2+2
readFile(mepH#"Directory","input_FPH",100000)
readFile(mepH#"Directory","input_start1",100000)
readFile(mepH#"Directory","input_start2",100000)

--------

mepH=newMultiparameterEigenvalueProblem({5,5,5},"DimensionDeficientRankNumberOfParameters")    ;
mepFiberProductHomotopy(mepH)
readFile(mepH#"Directory","start1",1000)
readFile(mepH#"Directory","start2",100)
readFile(mepH#"Directory","start3",100)
readFile(mepH#"Directory","nonsingular_solutions",1000)
readFile(mepH#"Directory")--2+2+2
readFile(mepH#"Directory","input_FPH",100000)
readFile(mepH#"Directory","input_start1",100000)
readFile(mepH#"Directory","input_start2",100000)

--------
mepH=newMultiparameterEigenvalueProblem({5,5},"DimensionDegree")    ;
mepFiberProductHomotopy(mepH)
SVD last mepH#"H1"
SVD  (mepH#"H1"#1+mepH#"H1"#2)
readFile(mepH#"Directory","input_FPH",100000)


--Homogenize
printingPrecision=100
mepH=newMultiparameterEigenvalueProblem({4,4},"DimensionDeficientRankOne")    ;
mepH#"TargetSystemBertiniConfigurations"={"MPType"=>1,"SecurityLevel"=>1}
mepH#"Homogenize"=true
mepFiberProductHomotopy(mepH)
readFile(mepH#"Directory","singular_solutions",100000)
readFile(mepH#"Directory","nonsingular_solutions",100000)

readFile(mepH#"Directory","start1",1000)
readFile(mepH#"Directory","start2",1000)

writeMultiparameterEigenvalueProblem(mepH);
runStartMultiparameterEigenvalueProblem(mepH);
getStartLines(mepH,2);
win=writeStartSolutionsFiberProductHomotopy(mepH);
4==#win
#first last win==mepH#"ExtrinsicDimension"//sum
#last last win==(mepH#"ExtrinsicDimension"//length)^2+1
moveB'File(mepH#"Directory","start_FPH","start",CopyB'File=>true);
runBertini(mepH#"Directory",NameB'InputFile=>"input_FPH");
readFile(mepH#"Directory")
readFile(mepH#"Directory","start1",1000)
readFile(mepH#"Directory","start2",1000)
readFile(mepH#"Directory","nonsingular_solutions",100000)


(ns1,ns2)=(5,5)
T=QQ[a,b]**QQ[apply(ns1,i->"x"|i)]**QQ[apply(ns2,i->"y"|i)]
xList=flatten entries basis(degree x1,T)
yList=flatten entries basis(degree y1,T)
randomMatrix=(m,n)->matrix apply(m,i->apply(n,j->random(1,1000)))
randomRankOneMatrix=(m,n)->transpose matrix {apply(m,i->random(1,1000))}* matrix {apply(n,i->random(1,1000))}
oneH=(m,n)->sub(randomMatrix(m,n),T)+a*sub(randomRankOneMatrix(m,n),T)+b*sub(randomRankOneMatrix(m,n),T)
I1=ideal ((oneH(ns1,ns1))*transpose basis(degree x1,T))
I2=ideal ((oneH(ns2,ns2))*transpose basis(degree y1,T))
I=I1+I2
decI=decompose I
netList oo
decI/degree

decompose 
J1=I1+ideal(10124+random({1,0,0},T))
decJ1=decompose J1
decJ1/(i->eliminate(xList,i))
oo/degree



eliminate(xList,J1)
oo/degree

netList decI
degree(I)


I=ideal(    det oneH(4,4),    det oneH(4,4))
degree I
radical I==I


R=QQ[x1,x2,x3,lam]
primaryDecomposition ideal(diagonalMatrix{lam,lam}*transpose matrix{{x1,x2}})
--slice x and tracking will work. 

primaryDecomposition ideal(matrix{{lam,1},{0,lam}}*transpose matrix{{x1,x2}})

loadPackage"Bertini"
(ideal(diagonalMatrix{lam,lam}*transpose matrix{{x1,x2}}))_*//bertiniPosDimSolve
(ideal(matrix{{lam,1},{0,lam}}*transpose matrix{{x1,x2}}))_*//bertiniPosDimSolve
nv=oo
peek nv#2#1


R=QQ[a1,a2,b1,b2,t]
H1= (1-t)*(2*a1+3*a2+7)    + 19*t*(a1-b1)+26*t*(a2-b2)
H2= (1-t)*(11*a1+13*a2+17) + 23*t*(a2-b2)+123*t*(a1-b1)
eliminate({a2},ideal(H1,H2))
