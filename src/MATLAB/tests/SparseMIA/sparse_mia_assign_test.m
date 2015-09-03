function tests = sparse_mia_assign_test
tests = functiontests(localfunctions);
end

function testSimpleAssign1(testCase)    
    dims=[5 7 6];
    lexOrder=[1 2 3];
    [C, DenseC]=sparse_mia_assign_test_work(dims,'i,j,k','i,j,k',lexOrder);
    verifyTrue(testCase,C.isSorted);    
    verifyEqual(testCase,C.dims,dims);
    verifyEqual(testCase,C.lexOrder,lexOrder);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleAssign1b(testCase)    
    dims=[5 7 6];
    lexOrder=[2 1 3];
    [C, DenseC]=sparse_mia_assign_test_work(dims,'i,j,k','i,j,k',lexOrder);    
    verifyTrue(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims);
    verifyEqual(testCase,C.lexOrder,lexOrder);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleAssign2(testCase)        
    dims=[5 7 6];
    lexOrder=[1 2 3];
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'i,j,k','i,k,j',lexOrder);    
    verifyFalse(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims([1,3,2]));
    verifyEqual(testCase,C.lexOrder,lexOrder([1,3,2]));
    verifyTrue(testCase,C==DenseC);    
end

function testSimpleAssign2b(testCase)        
    dims=[5 7 6];
    lexOrder=[3 2 1];
    reverse_permute_idx=[1 3 2];
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'i,j,k','i,k,j',lexOrder);    
    verifyFalse(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims([1,3,2]));
    verifyEqual(testCase,C.lexOrder,reverse_permute_idx(lexOrder));
    verifyTrue(testCase,C==DenseC);    
end

function testSimpleAssign3(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6 8 9],'j,k,i,l,m','m,i,l,j,k');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign1(testCase)       
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'!i,j,k','i,j,k');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign2(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'i,j,k','!i,j,k');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign3(testCase)       
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'i,!j,k','i,k,j');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign4(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'!i,!j,!k','i,k,j');
    verifyTrue(testCase,C==DenseC);
end