function [a, b, a2, b2] = bench_test(max_size, runs)
a = 0
b = 0
a2 = 0
b2 = 0

    for i = 0: runs
        step = i * max_size/runs
        x = rand(step, step)
        A = x * x'
        tic();
        [L, D] = ldltt(A);
        a = a + toc();
        a2 = a2 + norm(L*D*L' - A)
        tic()
        [L, D] = myldlt(A);
        b = b + toc();
        b2 = b2 + norm(L*D*L' - A)
    end
endfunction

[a, b, a2, b2] = bench_test(100, 50)