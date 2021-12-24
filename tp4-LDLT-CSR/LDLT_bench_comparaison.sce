function [] = LDLT_bench_comparaison(matrix_size_range, runs)
    filename1 = "data/LDLT3b.dat";
    filename2 = "data/LDLT1b.dat";
    [f1, mode_f1] = mopen(filename1, "w");
    [f2, mode_f2] = mopen(filename2, "w");
    
    for i = matrix_size_range
        printf("Running matrix size %d*%d\n", i,i);
        t = zeros(2,runs)
        for j = 1 : runs
            x = rand(i,i);
            A = x*x';

            tic();
            [L,D] = myLDLT3b(A); 
            t(1,j) = toc();
            n(1,j) = norm(A-L*D*L');
            
            tic();
            [L,D] = myLDLT1b(A); 
            t(2,j) = toc();
            n(2,j) = norm(A-L*D*L');
        end
        mean1 = mean(t(1,:));
        mean2 = mean(t(2,:));
        meannorm1 = mean(n(1,:));
        meannorm2 = mean(n(2,:));
        mfprintf(f1, "%d %.17f %.17f %.17f %.17f %.17f %.17f\n", i, mean1, min(t(1,:)), max(t(1,:)), meannorm1, min(n(1,:)), max(n(1,:)));
        mfprintf(f2, "%d %.17f %.17f %.17f %.17f %.17f %.17f\n", i, mean2, min(t(2,:)), max(t(2,:)), meannorm2, min(n(2,:)), max(n(2,:)));
    end
    mclose(f1);
    mclose(f2);
endfunction