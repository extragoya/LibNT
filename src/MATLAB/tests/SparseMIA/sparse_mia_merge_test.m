function tests = sparse_mia_merge_test
tests = functiontests(localfunctions);
end

function testSimpleMerge1(testCase)        
        
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 7 6],'ijk',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 7 6],'ijk',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[5 7 6],'ijk',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[5 7 6],'ijk',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 7 6],'ijk',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 7 6],'ijk',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[5 7 6],'ijk',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[5 7 6],'ijk',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    
end


function testSimpleMerge2(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 6 7],'ikj',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 6 7],'ikj',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[5 6 7],'ikj',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[5 6 7],'ikj',[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 6 7],'ikj',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[5 6 7],'ikj',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[5 6 7],'ikj',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[5 6 7],'ikj',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleMerge3(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[6 5 7],'kij',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[6 5 7],'kij',[2 3 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[6 5 7],'kij',[1 2 3],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[3 2 1],0,[6 5 7],'kij',[1 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[6 5 7],'kij',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],1,[6 5 7],'kij',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[6 5 7],'kij',[1 2 3],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6],'ijk',[1 2 3],0,[6 5 7],'kij',[1 2 3],0,0);
    verifyTrue(testCase,C==DenseC);
end

function testSimpleMerge4(testCase)    
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[1 2 3 4 5],1,[9 6 8 5 7],'miljk',[1 2 3 4 5],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[5 4 3 2 1],1,[9 6 8 5 7],'miljk',[5 4 3 2 1],0,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[1 5 4 2 3],0,[9 6 8 5 7],'miljk',[2 5 4 3 1],1,1);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[2 5 4 3 1],0,[9 6 8 5 7],'miljk',[1 5 4 2 3],0,1);
    verifyTrue(testCase,C==DenseC);

    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[4 3 2 1 5],1,[9 6 8 5 7],'miljk',[4 5 3 2 1],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[1 2 3 5 4],1,[9 6 8 5 7],'miljk',[1 2 3 5 4],0,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[4 5 3 2 1],0,[9 6 8 5 7],'miljk',[4 3 2 1 5],1,0);
    verifyTrue(testCase,C==DenseC);
    [C, DenseC]=sparse_mia_merge_test_work([5 7 6 8 9],'jkilm',[5 1 2 3 4],0,[9 6 8 5 7],'miljk',[5 1 2 3 4],0,0);
    verifyTrue(testCase,C==DenseC);
end
