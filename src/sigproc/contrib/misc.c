/*
 * Functions for 4x4 matrix inversion by Gaussian elimination. The operation
 * is usually run only once during the whole process execution. Therefore,
 * there is no need for optimizations for speed here. The code is adopted
 * from the gamesxmath's blog:
 * http://gamesxmath.tumblr.com/post/86495837013/inverting-a-4x4-matrix-using-c-code-double
 */

static float
invmtx_invf(int i, int j, const float *m)
{
    int o = 2 + (j - i);

    i += 4 + o;
    j += 4 - o;

#define e(a,b) m[((j + b) % 4) * 4 + ((i + a) % 4)]

    float inv =
        + e(+1,-1)*e(+0,+0)*e(-1,+1)
        + e(+1,+1)*e(+0,-1)*e(-1,+0)
        + e(-1,-1)*e(+1,+0)*e(+0,+1)
        - e(-1,-1)*e(+0,+0)*e(+1,+1)
        - e(-1,+1)*e(+0,-1)*e(+1,+0)
        - e(+1,-1)*e(-1,+0)*e(+0,+1);

    return (o % 2) ? inv : -inv;

#undef e
}


int
inverse_4x4_matrix(const float *m, float *out)
{
    float inv[16];
    int i, j, k;
    float D;
	
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            inv[j * 4 + i] = invmtx_invf(i, j, m);

    D = 0.;

    for (k = 0; k < 4; k++)
        D += m[k] * inv[k * 4];

    if (D == 0)
        return -1;
	
    D = 1. / D;
	
    for (i = 0; i < 16; i++)
        out[i] = inv[i] * D;

    return 0;
}
