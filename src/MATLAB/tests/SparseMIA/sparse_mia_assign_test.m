function tests = sparse_mia_assign_test
tests = functiontests(localfunctions);
end

function testSimpleAssign1(testCase)    
    dims=[5 7 6];
    lexOrder=[1 2 3];
    [C, DenseC]=sparse_mia_assign_test_work(dims,'ijk','ijk',lexOrder);
    verifyTrue(testCase,C.isSorted);    
    verifyEqual(testCase,C.dims,dims);
    verifyEqual(testCase,C.lexOrder,lexOrder);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleAssign1b(testCase)    
    dims=[5 7 6];
    lexOrder=[2 1 3];
    [C, DenseC]=sparse_mia_assign_test_work(dims,'ijk','ijk',lexOrder);    
    verifyTrue(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims);
    verifyEqual(testCase,C.lexOrder,lexOrder);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleAssign2(testCase)        
    dims=[5 7 6];
    lexOrder=[1 2 3];
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'ijk','ikj',lexOrder);    
    verifyFalse(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims([1,3,2]));
    verifyEqual(testCase,C.lexOrder,lexOrder([1,3,2]));
    verifyTrue(testCase,C==DenseC);    
end

function testSimpleAssign2b(testCase)        
    dims=[5 7 6];
    lexOrder=[3 2 1];
    reverse_permute_idx=[1 3 2];
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'ijk','ikj',lexOrder);    
    verifyFalse(testCase, C.isSorted);
    verifyEqual(testCase,C.dims,dims([1,3,2]));
    verifyEqual(testCase,C.lexOrder,reverse_permute_idx(lexOrder));
    verifyTrue(testCase,C==DenseC);    
end

function testSimpleAssign3(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6 8 9],'jkilm','miljk');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign1(testCase)       
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'!ijk','ijk');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign2(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'ijk','!ijk');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign3(testCase)       
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'i!jk','ikj');
    verifyTrue(testCase,C==DenseC);
end

function testTrickyAssign4(testCase)        
    [C, DenseC]=sparse_mia_assign_test_work([5 7 6],'!i!j!k','ikj');
    verifyTrue(testCase,C==DenseC);
end