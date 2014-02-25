function tests = dense_mia_merge_test
tests = functiontests(localfunctions);
end

function testSimpleMerge1(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','ijk',[1 2 3],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','ijk',[1 2 3],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleMerge2(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','ikj',[1 3 2],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','ikj',[1 3 2],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testSimpleMerge3(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','kij',[2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','kij',[2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleMerge4(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6 8 9],'jkilm','miljk',[4 5 2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6 8 9],'jkilm','miljk',[4 5 2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testTrickyMerge1(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'i!jk','ijk',[1 2 3],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'ijk','i!jk',[1 2 3],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyMerge2(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'!ijk','ik!j',[1 3 2],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'!ijk','ik!j',[1 3 2],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end
function testTrickyMerge3(testCase)    
    [C, c_data]=dense_mia_merge_test_work([5 7 6],'!!i!j!k','kij',[2 3 1],1);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
     [C, c_data]=dense_mia_merge_test_work([5 7 6],'!!i!j!k','kij',[2 3 1],0);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end