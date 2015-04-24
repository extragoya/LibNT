function tests = sparse_lattice_merge_test
tests = functiontests(localfunctions);
end

function testPlus1(testCase)

    m=5;n=5;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testPlus2(testCase)

    m=1;n=5;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testPlus3(testCase)

    m=5;n=1;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testPlus4(testCase)

    m=5;n=5;p=1;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testPlus5(testCase)

    m=1;n=1;p=1;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,1);    
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testMinus1(testCase)

    m=5;n=5;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,2);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testMinus2(testCase)

    m=1;n=5;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,2);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testMinus3(testCase)

    m=5;n=1;p=5;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,2);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testMinus4(testCase)

    m=5;n=5;p=1;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,2);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end

function testMinus5(testCase)

    m=1;n=1;p=1;
    [C, denseC]=sparse_lattice_merge_test_work(m,n,p,2);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    
    
end