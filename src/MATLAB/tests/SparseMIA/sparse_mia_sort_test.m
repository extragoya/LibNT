function tests = sparse_mia_sort_test
tests = functiontests(localfunctions);
end

function testSort1(testCase)

    dims=[5 5 5];
    [A, A2]=sparse_mia_sort_test_work(dims);
    verifyEqual(testCase,A.data,A2.data);
    verifyEqual(testCase,A.indices,A2.indices);    
    
end


function testSort2(testCase)

    dims=[1 1 1];
    [A, A2]=sparse_mia_sort_test_work(dims);
   
       
    verifyEqual(testCase,A.data,A2.data);
    verifyEqual(testCase,A.indices,A2.indices);
    
end