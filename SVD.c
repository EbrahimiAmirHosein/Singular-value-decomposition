#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define N 3

/**
 * @authors : AmirHosein ebrahimi - amirh.ebrahimi1377@gmail.com
 *            Hoda vafaei sefat   - h.vafaeisefat@gmail.com
 * 
 * @brief calculate SVD Singular value decomposition in c without any library:
 * 
 *  Singular value decomposition is a matrix factorization method that
 *  generalizes the eigendecomposition of a square matrix (n x n) to any matrix (n x m)
 *  SVD is similar to Principal Component Analysis (PCA)
 * 
 * ormula of SVD is:
 *                  A=UΣVᵗ, where :
 * 
 *                      * M-is original matrix we want to decompose
 *                      * U-is left singular matrix (columns are left singular vectors). U columns contain eigenvectors of matrix MMᵗ
 *                      * Σ-is a diagonal matrix containing singular (eigen)values  
 *                      * V-is right singular matrix (columns are right singular vectors). V columns contain eigenvectors of matrix MᵗM
 * 
 *                  Power iteration :
 *                      Power iteration starts with b₀ which might be a random vector. At every iteration this vector is updated using following rule:
 *                              b[k+1] = A*b[k] / || A*b[k] ||
 *                     First we multiply b₀ with original matrix A (Abₖ) and divide result with the norm (||Abₖ||).
 *                       We’ll continue until result has converged (updates are less than threshold).
 * 
 *                      Power method has few assumptions:
 * 
 *                          1- b₀ has a nonzero component in the direction of an eigenvector associated with the dominant eigenvalue. 
 *                             Initializing b₀ randomly minimizes possibility that this assumption is not fulfilled
 *                          2- matrix A has dominant eigenvalue which has strictly greater magnitude than other eigenvalues 
 *                              source : https://en.wikipedia.org/wiki/Power_iteration
 * 
 *                      These assumptions guarantee that algorithm converges to a reasonable result. 
 *                      The smaller is difference between dominant eigenvalue and second eigenvalue,
 *                      the longer it might take to converge
 * 
 * 
 *                      * To get more than just most dominant singular value from matrix, we could still use power iteration: 
 *                          We should remove dominant direction from the matrix and repeat finding most dominant singular value 
 *                              source : https://jeremykun.com/2016/05/16/singular-value-decomposition-part-2-theorem-proof-algorithm/
 * 
 *                      To do that we could subtract previous eigenvector(s) component(s) from the original matrix :
 *                          
 *                                      A_next = A-(singular_value₁)(u₁)(v₁)ᵗ
 * 
 * @param A matrix
 * @return ** SVD matrix 
 * @source https://towardsdatascience.com 
 */

void SVD(double A[][N])
{
    double svd_res[2*N+1][N];
    double eps = 0.99999;
    double ACopy[N][N];
    double AT[N][N];
    double AC[N][N];
    double singularValues[N];
    double u[N][N], vT[N][N], v[N][N], v_unnorm[N];
    double sigma = 0;
    double SVD_arr[2 * N + 1][N];

    double x_p[N];
    double lastV[N];
    double currentV[N];
    double A_T[N][N];
    double B[N][N];
    int iterations = 0;
    double norm = 0;
    double currentV_tmp[N];
    double VdotVold = 0;
    double unnormalized[10];
    double sumation = 0;

    // A.T
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            AT[j][i] = A[i][j];
        }
    }

// main loop
    for (size_t iter_i = 0; iter_i < N; iter_i++)
    {

        for (size_t i = 0; i < N; i++)
        {
            v_unnorm[i] = 0;
        }

        // make a copy of A
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                ACopy[i][j] = A[i][j];
            }
        }

        // Updating Acopy ->  Acopy = Acopy - ( ) singularValue * outer_product of(u, v) )
        for (size_t z = 0; z < iter_i; z++)
        {
            for (size_t i = 0; i < N; i++)
            {

                for (size_t j = 0; j < N; j++)
                {
                    ACopy[i][j] -= u[i][z] * v[j][z] * singularValues[z];
                }
            }
        }

    // calculating first svd1D # next singular vector using power method


        // randomUnitVector -> unnormalized
        for (size_t i = 0; i < N; i++)
        {
            unnormalized[i] = sqrt(log(((double)(rand()) + 1.) / ((double)(RAND_MAX) + 1.)) * -2) * cos(((double)(rand()) + 1.) / ((double)(RAND_MAX) + 1.) * 2 * 3.14);
        }

        for (size_t i = 0; i < N; i++)
        {
            sumation += unnormalized[i] * unnormalized[i];
        }

        for (size_t i = 0; i < N; i++)
        {
            unnormalized[i] /= sqrt(sumation);
        }
        sumation = 0;

        // x_p = unnormalized
        for (size_t i = 0; i < N; i++)
        {
            x_p[i] = unnormalized[i];
        }
        // currentV = x_p
        for (size_t i = 0; i < N; i++)
        {
            currentV[i] = x_p[i];
        }


    // calculating A Transpose and B = A @ A.T
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {

                A_T[j][i] = ACopy[i][j];
            }
        }
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                B[i][j] = 0;
            }
        }
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                for (size_t k = 0; k < N; k++)
                {

                    B[i][j] += ACopy[i][k] * A_T[k][j];
                }
            }
        }

    // main while until convergence
        while (1)
        {
            iterations++;
            VdotVold = 0;
            // lastV = currentV
            for (size_t i = 0; i < N; i++)
            {
                lastV[i] = currentV[i];
            }
            // reseting  currentV_tmp for later use
            for (size_t i = 0; i < N; i++)
            {
                currentV_tmp[i] = 0;
            }
            // currentV = B @ lastV
            for (size_t j = 0; j < N; j++)
            {
                for (size_t k = 0; k < N; k++)
                {
                    currentV_tmp[j] += B[j][k] * lastV[k];
                }
            }
            // currentV = currentV / norm_currentV
            // norm_currentV = sqrt (sum ( norm_currentV^2 ))
            norm = 0;
            for (size_t i = 0; i < N; i++)
            {
                norm += (currentV_tmp[i] * currentV_tmp[i]);
            }
            for (size_t i = 0; i < N; i++)
            {
                currentV_tmp[i] /= sqrt(norm);
            }

            for (size_t i = 0; i < N; i++)
            {
                currentV[i] = currentV_tmp[i];
            }

            // VdotVold = currentV * lastV -> convergence
            for (size_t i = 0; i < N; i++)
            {
                VdotVold += currentV[i] * lastV[i];
            }

            VdotVold = VdotVold > 0 ? VdotVold : -1 * VdotVold;

            if (VdotVold > eps) break; //convergence checking 
        }

        // current_V = first svd1D of A , now calculating other dimensions of SVD
        for (size_t i = 0; i < N; i++)
        {
            u[i][iter_i] = currentV[i];
        }

        //calculating v_unnormalized = np.dot(A.T, u)
        for (size_t i = 0; i < N; i++)
        {
            for (size_t k = 0; k < N; k++)
            {
                v_unnorm[i] += AT[i][k] * u[k][iter_i];
            }
        }
        // calculating norm(v_unnormalized)
        sigma = 0; //# next singular value

        for (size_t i = 0; i < N; i++)
        {
            sigma += (v_unnorm[i] * v_unnorm[i]);
        }
        for (size_t i = 0; i < N; i++)
        {
            v_unnorm[i] /= sqrt(sigma);
        }
        for (size_t i = 0; i < N; i++)
        {
            v[i][iter_i] = v_unnorm[i];
        }

        singularValues[iter_i] = sqrt(sigma);
    }

    //calculating V Transpose
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            vT[j][i] = v[i][j];
        }
    }

    /**
     * @brief storing singularValues in svd_res 
     * 
     * svd_res is (2N+1 , N) array which :
     * first row indicates all singular values of A or (S)
     * next N rows inidicates all u matrices of A or (V)
     * last N rows indicates all v matrices of A or (D)
     * 
     *                         svd_res
     *              ------------------------------- 
     *            0 |     All singular values     |     S
     *              -------------------------------
     *            1 |                             |
     *            : |                             |     V
     *            : |         u matrices          |
     *            N |                             |
     *              ------------------------------- 
    *           N+1 |                             |
     *            : |                             |     D
     *            : |         v matrices          |
     *           2N |                             |
     *              -------------------------------    
     * 
     */
    for (size_t i = 0; i < N; i++)
    {
        svd_res[0][i] = singularValues[i];
    }
    int row = 1;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            svd_res[row][j] = u[i][j];
        }
        row++;
    }
    row = N+1;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            svd_res[row][j] = vT[i][j];
        }
        row++;
    }
    for (size_t i = 0; i < 2*N+1; i++)
    {
        printf("\n");
        for (size_t j = 0; j < N; j++)
        {
            printf("svd : %f ", svd_res[i][j]);
        }
    }
}


// int main()
// {
//     double B[N][N] = {
//         {7, 2, 9},
//         {1, 2, 3},
//         {5, 3, 1}};
//     SVD(B);

//     return 0;
// }
