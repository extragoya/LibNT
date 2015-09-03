function tests = dense_mia_assign_test
tests = functiontests(localfunctions);
end

function testSimpleAssign1(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'i,j,k','i,j,k',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleAssign2(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'i,j,k','i,k,j',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleAssign3(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6 8 9],'j,k,i,l,m','m,i,l,j,k',[5 3 4 1 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign1(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!i,j,k','i,j,k',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign2(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!i,j,k','!i,j,k',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign3(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'i,!j,k','i,k,j',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign4(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!i,!j,!k','i,k,j',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end