function tests = sparse_mia_permute_test
tests = functiontests(localfunctions);
end

function testLinIdx3_1(testCase)

    
    newIdx=[3 2 1];
    dims=[5 5 5];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);    
    verifyEqual(testCase,A.data,A2.data);    
    
end

function testLinIdx3_2(testCase)

    
    newIdx=[3 1 2];
    dims=[5 5 5];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);    
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx3_3(testCase)

    
    newIdx=[2 1 3];
    dims=[5 5 5];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx3_4(testCase)

    
    newIdx=[2 1 3];
    dims=[6 8 9];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx3_5(testCase)

    
    newIdx=[1 2 3];
    dims=[6 8 9];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx4_1(testCase)

    
    newIdx=[2 1 3 4];
    dims=[5 5 5 5];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);  
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx4_2(testCase)

    
    newIdx=[2 1 3 4];
    dims=[10 4 5 17];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    verifyEqual(testCase,A.indices,A2.indices);
    verifyEqual(testCase,A.data,A2.data); 
    
end

function testLinIdx4_3(testCase)

    
    newIdx=[2 1 3 4];
    dims=[10 4 5 17];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims);   
    newIdx=[1 2 3 4];
    [A, A2]=sparse_mia_permute_test_work(newIdx,dims,A,A2);   
    verifyEqual(testCase,A.indices,A2.indices); 
    verifyEqual(testCase,A.data,A2.data); 
    
end
