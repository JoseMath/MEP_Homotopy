-- Add path
path = append(path, "/home/dujinhong/M2/")
loadPackage("Bertini",Reload=>true,Configuration=>{"BERTINIexecutable"=>"/home/dujinhong/bertini/bin/bertini"})
loadPackage("MultiparameterEigenvalueProblemHomotopy",Reload=>true)
printingPrecision=100

-- Table 4
n = 5; -- 5, 15, 30, 50
extDim = {n,n};
randMatrix = (i,j)-> matrix for i to i-1 list for j to j-1 list random CC;

theDirFiber = concatenate {"/home/dujinhong/result/table4/fiber/", toString(n), "/"};
theDirDiag = concatenate {"/home/dujinhong/result/table4/diag/", toString(n), "/"};
mkdir theDirFiber;
mkdir theDirDiag;
ComputeTime = for i when i < 10 list (
    -- fiber product
    subDir = concatenate {theDirFiber, toString(i), "/"};
    mkdir subDir;
    mepH = newMultiparameterEigenvalueProblem(extDim,"GenericExtrinsic");
    mepH#"Directory" = subDir;

    t1 = elapsedTiming(
        writeMultiparameterEigenvalueProblem(mepH);
        runStartMultiparameterEigenvalueProblem(mepH);
        writeStartSolutionsFiberProductHomotopy(mepH);
        moveB'File(mepH#"Directory","start_FPH","start",CopyB'File=>true);
        runBertini(mepH#"Directory",NameB'InputFile=>"input_FPH");
        readFile(mepH#"Directory");
    );


    -- diagonal coefficients
    subDir = concatenate {theDirDiag, toString(i), "/"};
    mkdir subDir;

    allD = apply(#extDim,i->randMatrix(2,extDim#i));
    allH = {mepH#"H1", mepH#"H2"};
    t2 = elapsedTiming(
    diagonalCoefficientHomotopy(allD, allH, subDir);    
    );    

    (t1#0,t2#0)
    )

theDir = "/home/dujinhong/result/table4/";
SFile = openOut(concatenate {theDir, "time_", toString(n), ".txt"}); 
scan(ComputeTime,x->(
    SFile<<x#0<<","<<x#1<<endl;	    
    ));
close SFile;



-- Table 5 QMEP
n = 2; -- 2, 5, 10, 20, 40
Tmpdim = {n,n}
extDim = {3*n,3*n};
randMatrix = (i,j)-> matrix for i to i-1 list for j to j-1 list random CC;

theDirFiber = concatenate {"/home/dujinhong/result/table5/fiber/", toString(n), "/"};
theDirDiag = concatenate {"/home/dujinhong/result/table5/diag/", toString(n), "/"};
mkdir theDirFiber;
mkdir theDirDiag;
ComputeTime = for i when i < 10 list (
    
    B = apply(Tmpdim,n->apply(#Tmpdim+1,i->randMatrix(n,n)));
    C = apply(Tmpdim,n->apply(#Tmpdim+1,i->randMatrix(n,n)));

    H1 = {(B#0#0 | B#0#1 | B#0#2) || (id_(CC^n)*0 | -id_(CC^n) | id_(CC^n)*0) || (id_(CC^n)*0 | id_(CC^n)*0 | -id_(CC^n)),
            (id_(CC^n)*0 | B#1#0 | B#1#1) || (id_(CC^n) | id_(CC^n)*0 | id_(CC^n)*0) || (id_(CC^n)*0 | id_(CC^n)*0 | id_(CC^n)*0),
            (id_(CC^n)*0 | id_(CC^n)*0 | B#1#2) || (id_(CC^n)*0 | id_(CC^n)*0 | id_(CC^n)*0) || (id_(CC^n) | id_(CC^n)*0 | id_(CC^n)*0)};
    H2 = {(C#0#0 | C#0#1 | C#0#2) || (id_(CC^n)*0 | -id_(CC^n) | id_(CC^n)*0) || (id_(CC^n)*0 | id_(CC^n)*0 | -id_(CC^n)),
            (id_(CC^n)*0 | C#1#0 | C#1#1) || (id_(CC^n) | id_(CC^n)*0 | id_(CC^n)*0) || (id_(CC^n)*0 | id_(CC^n)*0 | id_(CC^n)*0),
            (id_(CC^n)*0 | id_(CC^n)*0 | C#1#2) || (id_(CC^n)*0 | id_(CC^n)*0 | id_(CC^n)*0) || (id_(CC^n) | id_(CC^n)*0 | id_(CC^n)*0)};

    -- fiber product
    subDir = concatenate {theDirFiber, toString(i), "/"};
    mkdir subDir;
    mepH=newMultiparameterEigenvalueProblem(extDim,"GenericExtrinsic");
    mepH#"H1" = H1;
    mepH#"H2" = H2;
    mepH#"Directory"=subDir;

    t1 = elapsedTiming(
        writeMultiparameterEigenvalueProblem(mepH);
        runStartMultiparameterEigenvalueProblem(mepH);
        writeStartSolutionsFiberProductHomotopy(mepH);
        moveB'File(mepH#"Directory","start_FPH","start",CopyB'File=>true);
        runBertini(mepH#"Directory",NameB'InputFile=>"input_FPH");
        readFile(mepH#"Directory");
    );


    -- diagonal coefficients
    subDir = concatenate {theDirDiag, toString(i), "/"};
    mkdir subDir;

    allD = apply(#extDim,i->randMatrix(2,extDim#i));
    allH = {H1, H2};
    t2 = elapsedTiming(
    diagonalCoefficientHomotopy(allD, allH, subDir);    
    );    

    (t1#0,t2#0)
    )

theDir = "/home/dujinhong/result/table5/";
SFile = openOut(concatenate {theDir, "time_", toString(n), ".txt"}); 
scan(ComputeTime,x->(
    SFile<<x#0<<","<<x#1<<endl;	    
    ));
close SFile;
