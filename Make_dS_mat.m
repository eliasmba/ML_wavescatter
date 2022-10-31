function dS = Make_dS_mat(R, centers, k, N_multi)
N_oneblock = (N_multi+1);
N = length(centers(1,:));
dS = zeros(N_oneblock*N);

for i=1:N
    for j=1:N
        if i==j
            dS((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= make_dS_diag(R,R,k,N_multi);
        else
            dS((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= make_dS_offdiag(centers(:,i),centers(:,j),R, k,N_multi);
        end
    end
end 