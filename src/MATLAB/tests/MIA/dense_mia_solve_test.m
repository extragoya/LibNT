function tests = dense_mia_solve_test
tests = functiontests(localfunctions);
end

function testSimpleSolve1(testCase)    
    [B, B_test]=dense_mia_solve_test_work([5 10 10],[10 10 5],'!i,j,k','j,l,!i',0);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end

function testSimpleSolve2(testCase)    
    [B, B_test]=dense_mia_solve_test_work([10 5 10 ],[5 10 10],'j,!i,k','!i,l,j',0);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end

function testSimpleSolve3(testCase)    
    [B, B_test]=dense_mia_solve_test_work([4 4 5 4 4],[4 5 4 4 4],'m,j,!i,n,k','n,!i,l,o,j',0);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end

function testLSQRSolve1(testCase)    
    [B, B_test]=dense_mia_solve_test_work([5 20 10],[20 10 5],'!i,j,k','j,l,!i',1);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end

function testLSQRSolve2(testCase)    
    [B, B_test]=dense_mia_solve_test_work([20 5 10],[5 10 20],'j,!i,k','!i,l,j',1);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end

function testLSQRSolve3(testCase)    
    [B, B_test]=dense_mia_solve_test_work([4 8 5 8 4],[8 5 4 4 8],'m,j,!i,n,k','n,!i,l,o,j',1);
    verifyEqual(testCase,B.data,B_test.data,'AbsTol',1000*eps );
end