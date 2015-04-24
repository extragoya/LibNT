function tests = dense_mia_mult_test
tests = functiontests(localfunctions);
end

function testSimpleMult1(testCase)    
    [C, c_data]=dense_mia_mult_test_work([5 7 6],[7 8 5],'!ijk','jl!i',3,2,1,2,1,3);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end


function testMult1(testCase)    
    [C, c_data]=dense_mia_mult_test_work([5 7 6 8],[4 7 8 5],'!ijkl','mjl!i',3,[2 4],1,1,[2 3],4);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testMult2(testCase)    
    [C, c_data]=dense_mia_mult_test_work([5 7 6 8],[4 8 7 5],'!ilkj','mjl!i',3,[2 4],1,1,[3 2],4);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testMultSmallDim1(testCase)    
    [C, c_data]=dense_mia_mult_test_work([5 6 ],[4 5],'!ik','m!i',2,[],1,1,[],2);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testMultSmallDim2(testCase)    
    [C, c_data]=dense_mia_mult_test_work([5 7 8],[7 8 5],'!ijl','jl!i',[],[2 3],1,[],[1 2],3);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end

function testMultSmallDim3(testCase)    
    [C, c_data]=dense_mia_mult_test_work([7 6 8],[4 7 8],'jkl','mjl',2,[1 3],[],1,[2 3],[]);
    verifyEqual(testCase,C.data,c_data,'AbsTol',100*eps );
end