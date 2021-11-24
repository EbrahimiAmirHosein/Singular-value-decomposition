# -Singular-value-decomposition-in-c
SVD algotihm in C without any library

calculate SVD Singular value decomposition in c without any library:
  
   Singular value decomposition is a matrix factorization method that
   generalizes the eigendecomposition of a square matrix (n x n) to any matrix (n x m)
   SVD is similar to Principal Component Analysis (PCA)
  
  formula of SVD is:
                   
                   A=UΣVᵗ, where :
  
                       * M-is original matrix we want to decompose
                       * U-is left singular matrix (columns are left singular vectors). U columns contain eigenvectors of matrix MMᵗ
                       * Σ-is a diagonal matrix containing singular (eigen)values  
                       * V-is right singular matrix (columns are right singular vectors). V columns contain eigenvectors of matrix MᵗM
  
                   Power iteration :
                       Power iteration starts with b₀ which might be a random vector. At every iteration this vector is updated using following rule:
                               b[k+1] = A*b[k] / || A*b[k] ||
                      First we multiply b₀ with original matrix A (Abₖ) and divide result with the norm (||Abₖ||).
                        We’ll continue until result has converged (updates are less than threshold).
  
                       Power method has few assumptions:
  
                           1- b₀ has a nonzero component in the direction of an eigenvector associated with the dominant eigenvalue. 
                              Initializing b₀ randomly minimizes possibility that this assumption is not fulfilled
                           2- matrix A has dominant eigenvalue which has strictly greater magnitude than other eigenvalues 
                               source : https://en.wikipedia.org/wiki/Power_iteration
  
                       These assumptions guarantee that algorithm converges to a reasonable result. 
                       The smaller is difference between dominant eigenvalue and second eigenvalue,
                       the longer it might take to converge
  
  
                       * To get more than just most dominant singular value from matrix, we could still use power iteration: 
                           We should remove dominant direction from the matrix and repeat finding most dominant singular value 
                               source : https://jeremykun.com/2016/05/16/singular-value-decomposition-part-2-theorem-proof-algorithm/
  
                       To do that we could subtract previous eigenvector(s) component(s) from the original matrix :
                           
                                       A_next = A-(singular_value₁)(u₁)(v₁)ᵗ
                                       
  This algorithm was created for personal purpose , Implementing SFA and NICA algorith , all other explanation needed are described in code.
