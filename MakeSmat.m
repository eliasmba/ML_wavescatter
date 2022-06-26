function S = MakeSmat(R,r,  centers, k, N_multi)
N = length(centers);
N_oneblock = (N_multi+1); 
S = zeros(N_oneblock*N);

for i=1:N
    for j=1:N
        if i==j
            S((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= makeS_diag(R, r,k,N_multi);
        else
            S((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= makeS_offdiag(centers(i),centers(j),R, k,N_multi);
        end
    end
end 