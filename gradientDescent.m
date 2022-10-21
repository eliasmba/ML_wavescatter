function[y, nor] = gradientDescent(U, Y, z, x, k, N)

YtUT = (pinv(Y)*U);

