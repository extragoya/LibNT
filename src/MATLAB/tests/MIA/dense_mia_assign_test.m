function tests = dense_mia_assign_test
tests = functiontests(localfunctions);
end

function testSimpleAssign1(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'ijk','ijk',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleAssign2(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'ijk','ikj',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testSimpleAssign3(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6 8 9],'jkilm','miljk',[5 3 4 1 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign1(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!ijk','ijk',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign2(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!ijk','!ijk',[1 2 3]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign3(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'i!jk','ikj',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testTrickyAssign4(testCase)    
    [C, c_data]=dense_mia_assign_test_work([5 7 6],'!i!j!k','ikj',[1 3 2]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end