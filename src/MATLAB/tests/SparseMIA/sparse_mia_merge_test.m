function tests = sparse_mia_merge_test
tests = functiontests(localfunctions);
end

function testSimpleMerge1(testCase)        
        
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 7 6],'i,j,k',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 7 6],'i,j,k',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[5 7 6],'i,j,k',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[5 7 6],'i,j,k',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 7 6],'i,j,k',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 7 6],'i,j,k',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[5 7 6],'i,j,k',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[5 7 6],'i,j,k',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    
end


function testSimpleMerge2(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 6 7],'i,k,j',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 6 7],'i,k,j',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[5 6 7],'i,k,j',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[5 6 7],'i,k,j',[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 6 7],'i,k,j',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[5 6 7],'i,k,j',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[5 6 7],'i,k,j',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[5 6 7],'i,k,j',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleMerge3(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[6 5 7],'k,i,j',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[6 5 7],'k,i,j',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[6 5 7],'k,i,j',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[3 2 1],0,[6 5 7],'k,i,j',[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[6 5 7],'k,i,j',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],1,[6 5 7],'k,i,j',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[6 5 7],'k,i,j',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'i,j,k',[1 2 3],0,[6 5 7],'k,i,j',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleMerge4(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[1 2 3 4 5],1,[9 6 8 5 7],'m,i,l,j,k',[1 2 3 4 5],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[5 4 3 2 1],1,[9 6 8 5 7],'m,i,l,j,k',[5 4 3 2 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[1 5 4 2 3],0,[9 6 8 5 7],'m,i,l,j,k',[2 5 4 3 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[2 5 4 3 1],0,[9 6 8 5 7],'m,i,l,j,k',[1 5 4 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[4 3 2 1 5],1,[9 6 8 5 7],'m,i,l,j,k',[4 5 3 2 1],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[1 2 3 5 4],1,[9 6 8 5 7],'m,i,l,j,k',[1 2 3 5 4],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[4 5 3 2 1],0,[9 6 8 5 7],'m,i,l,j,k',[4 3 2 1 5],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'j,k,i,l,m',[5 1 2 3 4],0,[9 6 8 5 7],'m,i,l,j,k',[5 1 2 3 4],0,0);
    verifyTrue(testCase,C==DenseC);
end
