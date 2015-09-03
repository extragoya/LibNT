function tests = sparse_mia_mult_test
tests = functiontests(localfunctions);
end

function testSimpleMult1(testCase)    
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[1 2 3],[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[2 3 1],[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[1 2 3],[3 2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[2 3 1],[3 2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[1 2 3],[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[2 3 1],[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[1 2 3],[3 2 1],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6],[7 8 5],'!i,j,k','j,l,!i',[2 3 1],[3 2 1],0,0);
    verifyTrue(testCase,C==DenseC);
end


function testMult1(testCase)        
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[1 2 3 4],[1 2 3 4],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[4 3 2 1],[1 2 3 4],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[1 2 3 4],[3 4 2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[2 3 1 4],[3 1 4 2],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[1 2 3 4],[1 2 3 4],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[4 3 2 1],[1 2 3 4],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[1 2 3 4],[3 4 2 1],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!i,j,k,l','m,j,l,!i',[2 3 1 4],[3 1 4 2],0,0);
    verifyTrue(testCase,C==DenseC);
end

function testMult2(testCase)    
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[1 2 3 4],[1 2 3 4],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[4 3 2 1],[1 2 3 4],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[1 2 3 4],[3 4 2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[2 3 1 4],[3 1 4 2],1,1);
    verifyTrue(testCase,C==DenseC);    
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[1 2 3 4],[1 2 3 4],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[4 3 2 1],[1 2 3 4],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[1 2 3 4],[3 4 2 1],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!i,l,k,j','m,j,l,!i',[2 3 1 4],[3 1 4 2],0,1);
    verifyTrue(testCase,C==DenseC);    
end

function testMultSmallDim1(testCase)    
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[1 2],[1 2],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[2 1],[1 2],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[1 2],[2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[2 1],[2 1],1,1);
    verifyTrue(testCase,C==DenseC);    
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[1 2],[1 2],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[2 1],[1 2],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[1 2],[2 1],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 6 ],[4 5],'!i,k','m,!i',[2 1],[2 1],0,0);
    verifyTrue(testCase,C==DenseC);    
    
end

function testMultSmallDim2(testCase)    
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[1 2 3],[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[2 3 1],[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[1 2 3],[3 2 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[2 1 3],[1 3 2],1,1);
    verifyTrue(testCase,C==DenseC);    
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[1 2 3],[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[2 3 1],[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[1 2 3],[3 2 1],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([5 7 8],[7 8 5],'!i,j,l','j,l,!i',[2 1 3],[1 3 2],0,0);
    verifyTrue(testCase,C==DenseC);    
    
end

function testMultSmallDim3(testCase)    
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[1 2 3],[1 2 3],1,1);    
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[2 3 1],[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[1 2 3],[3 2 1],1,1);
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[2 1 3],[1 3 2],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[1 2 3],[1 2 3],0,1);    
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[2 3 1],[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[1 2 3],[3 2 1],0,1);
    verifyTrue(testCase,C==DenseC);   
    [C, DenseC]=sparse_mia_mult_test_work([7 6 8],[4 7 8],'j,k,l','m,j,l',[2 1 3],[1 3 2],0,0);
    verifyTrue(testCase,C==DenseC);
end