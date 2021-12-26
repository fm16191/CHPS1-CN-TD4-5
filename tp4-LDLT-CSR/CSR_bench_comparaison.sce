function [] = CSR_bench_comparaison(m, matrix_size_range, runs, density)
    filename1 = "data/CSR.dat";
    [f1, mode_f1] = mopen(filename1, "w");

    // Variables d'environnement
    p = 1;
    
    for i = matrix_size_range
        printf("Running matrix size %d*%d\n", m,i);
        t = zeros(2,runs)
        for j = 1 : runs

            n = i;
            
            // Mise en place
            sp = sprand(m,n,density);
            A = full(sp);
            v = grand(n,1, "bin", 1, p);
            
            tic();
            x = (A*v)';
            t(1,j) = toc();
            
            tic();
            [AX, AI, AJ] = csmtCSR(A);
            t(2,j) = toc();
            
            tic();
            xex = mCSRv(AX, AI, AJ, v);
            t(3,j) = toc();

            np(j) = norm(x - xex);

        end
        mean_t1 = mean(t(1,:));
        mean_t2 = mean(t(2,:));
        mean_t3 = mean(t(3,:));
        mean_np = mean(np(1,:));

        mfprintf(f1, "%d %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f\n", i, mean_t1, min(t(1,:)), max(t(1,:)), mean_t2, min(t(2,:)), max(t(2,:)), mean_t3, min(t(3,:)), max(t(3,:)), mean_np, min(np(1,:)), max(np(1,:)));
    end
    mclose(f1);
endfunction

function [] = CSR_bench_comparaison_density(m, n, runs, density_range)
    filename1 = "data/CSR_density.dat";
    [f1, mode_f1] = mopen(filename1, "w");

    // Variables d'environnement
    p = 1;
    density = 0.1;

    printf("= Matrix size : %d*%d =\n \n", m, n)

    it = 1

    for density=density_range
        printf("Running density : %f\n", density);
        t = zeros(2,runs)
        for j = 1 : runs

            // Mise en place
            sp = sprand(m,n,density);
            A = full(sp);
            v = grand(n,1, "bin", 1, p);
            
            tic();
            x = (A*v)';
            t(1,j) = toc();
            
            tic();
            [AX, AI, AJ] = csmtCSR(A);
            t(2,j) = toc();
            
            tic();
            xex = mCSRv(AX, AI, AJ, v);
            t(3,j) = toc();

            np(j) = norm(x - xex);

        end
        mean_t1 = mean(t(1,:));
        mean_t2 = mean(t(2,:));
        mean_t3 = mean(t(3,:));
        mean_np = mean(np(1,:));

        mfprintf(f1, "%f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f %.17f\n", density, mean_t1, min(t(1,:)), max(t(1,:)), mean_t2, min(t(2,:)), max(t(2,:)), mean_t3, min(t(3,:)), max(t(3,:)), mean_np, min(np(1,:)), max(np(1,:)));
    end
    mclose(f1);
endfunction