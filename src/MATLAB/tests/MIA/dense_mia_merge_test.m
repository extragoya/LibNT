function tests = dense_mia_merge_test
tests = functiontests(localfunctions);
end

function testSimpleMerge1(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','i,j,k',[1 2 3],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','i,j,k',[1 2 3],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleMerge2(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','i,k,j',[1 3 2],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','i,k,j',[1 3 2],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testSimpleMerge3(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','k,i,j',[2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','k,i,j',[2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleMerge4(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m','m,i,l,j,k',[4 5 2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m','m,i,l,j,k',[4 5 2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testTrickyMerge1(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,!j,k','i,j,k',[1 2 3],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i,j,k','i,!j,k',[1 2 3],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyMerge2(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'!i,j,k','i,k,!j',[1 3 2],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'!i,j,k','i,k,!j',[1 3 2],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testTrickyMerge3(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'!!i,!j,!k','k,i,j',[2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'!!i,!j,!k','k,i,j',[2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end