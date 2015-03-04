function tests = sparse_dense_lattice_mult_test
tests = functiontests(localfunctions);
end

function testMult1(testCase)
    m1=5;n1=6;p=2;m2=6;n2=10;
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,0);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
end


function testMult2(testCase)
    m1=5;n1=6;p=1;m2=6;n2=10;
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,0);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
end

function testMult3(testCase)
    m1=5;n1=1;p=2;m2=1;n2=10;
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,0);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
end


function testMult4(testCase)
    m1=1;n1=6;p=2;m2=6;n2=1;
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,0);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
end

function testMult5(testCase)
    m1=1;n1=25;p=1;m2=25;n2=1;
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,0);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
    [C, denseC]=sparse_dense_lattice_mult_test_work(m1,n1,m2,n2,p,1);
    verifyTrue(testCase,fuzzy_eq(C,denseC,0));
end
