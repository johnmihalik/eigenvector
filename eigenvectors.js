/**
 *  Power Iteration algorithm for eigenvector decomposition.  Computes the greatest eigenvalue of
 *  the input matrix and its corresponding eigenvector.  Input matrix MUST be diagonalizable.
 *
 *  See: https://en.wikipedia.org/wiki/Power_iteration
 *
 *
 * @param  {[][]} M   An m x n matrix
 * @param  {Number}   The convergence tolerance
 * @param  {Number}   The maximum iteration to allow before convergence
 * @return {[]}       An array with the eigenvalue and the eigenvector
 */
function power_iteration(M, tolerance=0.00001, max_iterations=1000) {

    let rank = M.length;

    // Initialize the first guess pf the eigenvector to a row vector of the sqrt of the rank
    let eigenvector = new Array(rank).fill(0).map(() => new Array(1).fill(Math.sqrt(rank)));

    // Compute the corresponding eigenvalue
    let eigenvalue = eigenvalue_of_vector(M, eigenvector);

    let epsilon = 1.0;
    let iterations = 0;
    do {

        let old_eigenvalue = deep_copy(eigenvalue);

        // Multiply the Matrix M by the guessed eigenveector
        let Mv = multiply_matrices(M,eigenvector);

        // Normalize the eigenvector to unit length
        eigenvector = normalize_vector(Mv);

        // Calculate the associated eigenvalue with the eigenvector (transpose(v) * M * v)
        eigenvalue = eigenvalue_of_vector(M, eigenvector);

        // Calculate the epsilon of the differences
        epsilon = Math.abs( eigenvalue - old_eigenvalue);
        iterations++;

    } while (epsilon > tolerance && iterations < max_iterations);

    return [eigenvalue, eigenvector];
}

//
// Matrix M MUST be a square and symmetric!

/**
 *  Compute all of the eigenvectors/values for a given matrix.  Uses shifting deflation
 *  with power iteration.
 *
 *  See: http://mlwiki.org/index.php/Power_Iteration#Finding_Other_Eigenvectors
 *
 *
 * @param  {[][]} M   An m x n matrix
 * @param  {Number}   The convergence tolerance
 * @param  {Number}   The maximum iteration to allow before convergence
 * @return {[]}       An array with the eigenvalues and the eigenvectors
 */
function eigens(M, tolerance=0.0001, max_iterations=1000) {

    let eigenvalues = [];
    let eigenvectors = [];

    for (let i = 0; i < M.length; i++ ) {

        // Compute the remaining most prominent eigenvector of the matrix M
        let result = power_iteration(M, tolerance, max_iterations);

        // Separate the eigenvalue and vector from the return array
        let eigenvalue = result[0];
        let eigenvector = result[1];

        eigenvalues[i] = eigenvalue;
        eigenvectors[i] = flatten_vector(eigenvector);

        // Now remove or peel off the last eigenvector
        M = shift_deflate(M, eigenvalue, eigenvector);
    }

    return [eigenvalues, eigenvectors];
}


/**
 * Print a matrix nicely for the console.
 *
 * @param  {Array} matrix  The m x n matrix to output.  Matrices are simple 2 dimensional arrays.
 */
// Print a matrix formatted nicely
function pretty_print_matrix(mat) {
    for ( let i = 0; i < mat.length; i++) {
        process.stdout.write("[ ");
        for ( let j = 0; j < mat[i].length; j ++) {
            if ( j != 0) {
                process.stdout.write(", ");
            }
            process.stdout.write( " " + mat[i][j] + " ");
        }
        process.stdout.write("] \n");
    }
}

/**
 *  Extract a column vector at column index x from the matrix.
 *
 * @param  {[][]} matrix  The m x n matrix to output.  Matrices are simple 2 dimensional arrays.
 * @param  {Number} x     The column index to extract
 * @return {[][]}         The column vector of the matrix (one element per row)
 */
function column(mat, x) {
    let col = new Array(mat.length).fill(0).map(() => new Array(1).fill(0));
    for (let i = 0; i < mat.length; i ++) {
        col[i][0] = mat[i][x];
    }
    return col;
}

/**
 *  Multiply 2 matrices together.  Note matrices must be compatible, n == r.  Output is a m x s matrix.
 *
 * @param  {[][]} m1  An m x n matrix
 * @param  {[][]} m2  An r x s matrix
 * @return {[][]}     The output product
 */
function multiply_matrices(m1, m2) {
    var result = [];
    for (var i = 0; i < m1.length; i++) {
        result[i] = [];
        for (var j = 0; j < m2[0].length; j++) {
            var sum = 0;
            for (var k = 0; k < m1[0].length; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }

    // If the result is a scalar dot product, just return the scalar value
    // This simplifies the results with vector dot products
    if (result.length == 1 && result[0].length == 1) {
        return result[0][0];
    }
    return result;
}

/**
 *  Multiple a matrix and a simple scalar value
 *
 * @param  {[][]} m    An m x n matrix
 * @param  {Number} s  The scalar to multiply the entries by
 * @return {[][]}      The output product
 */
function multiply_matrix_scalar(m, s) {
    let result = [];
    for (var i = 0; i < m.length; i++) {
        result[i] = [];
        for (let j = 0; j < m[0].length; j++) {
            result[i][j] = m[i][j] * s;
        }
    }
    return result;
}

/**
 *  Transpose a matrix
 *
 * @param  {[][]} matrix   An m x n matrix
 * @return {[][]}          The output matrix with rows and columns transposed
 */
function transpose(matrix) {
    return matrix[0].map((col, i) => matrix.map(row => row[i]));
}


/**
 *  deep copy an array object with nested elements
 *
 * @param  {Object}   The object to copy
 * @return {Object}   The copied output object
 */
function deep_copy(a) {
    return JSON.parse(JSON.stringify(a));
}


/**
 *  Sum the squared differences of 2 vectors
 *
 * @param  {[][]} v1   An m x 1 vector
 * @param  {[][]} v2   An m x 1 vector
 * @return {Number}    The sum of the squared differences
 */
function squared_difference(v1, v2) {
    let sum = 0.0;
    for (let i = 0; i < v1.length; i ++) {
        sum = sum + Math.pow( v1[i] - v2[i], 2 );
    }
    return sum;
}


/**
 *  Subtract matrices A - B  Note matrices must be compatible equal sizes
 *
 * @param  {[][]} m1  An m x n matrix
 * @param  {[][]} m2  An r x s matrix
 * @return {[][]}     The output difference
 */
function subtract_matrices(A,B) {
    let result = [];
    for (let i = 0; i < A.length; i++) {
        result[i] = [];
        for (var j = 0; j < A[0].length; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

//
/**
 *  Utility function to simplify a vector into a standard array
 *
 * @param  {[][]} v1   An m x 1 vector
 * @return {Array}    A simple array of the vector values
 */
function flatten_vector(v) {
    let v_new = [];
    for (let i = 0; i < v.length; i++) {
        v_new[i] = v[i][0];
    }
    return v_new;
}


/**
 *   Uses deflation to compute remaining eigenvectors.  At each iteration the most
 *   prominent eigenvector is computed.  This eigenvectors contribution is then
 *   removed from the original input matrix iteratively until all eigenvectors have
 *   been computed.
 *
 *   See: https://math.stackexchange.com/questions/768882/power-method-for-finding-all-eigenvectors
 *
 * @param  {[][]} m1  An m x n matrix
 * @param  {[][]} m2  An r x s matrix
 * @return {[][]}     The output difference
 */
function shift_deflate(M, eigenvalue, eigenvector)  {
    let len = Math.sqrt( multiply_matrices(transpose(eigenvector), eigenvector)  );
    let U = multiply_matrix_scalar(eigenvector, 1.0/len);
    let delta = multiply_matrix_scalar( multiply_matrices(U, transpose(U)) , eigenvalue);
    let M_new = subtract_matrices(M, delta);
    return M_new;
}

/**
 *  Computes the eigenvalue for the input eigenvector.
 *
 *     lambda = transpose(v) * M * v
 *
 * @param M
 * @param eigenvector
 * @return {*[][]}
 */
function eigenvalue_of_vector(M, eigenvector) {
    // Xt * M * x
    ev = multiply_matrices( multiply_matrices(transpose(eigenvector), M ), eigenvector);
    return ev;
}


/**
 *  Normalize a vector, divide by the length
 *
 * @param v
 * @return {*[][]}
 */
function normalize_vector(v) {
    let len = Math.sqrt( multiply_matrices(transpose(v), v));
    let unit_vector = multiply_matrix_scalar(v, 1.0 / len);
    return unit_vector;
}


function test() {

    console.log("Testing eigenvector decomposition...");

    const M = [[1,2,3], [2,7,9], [3,9,8]];
    console.log("Input Matrix: ");
    console.log(M);

    /*

     [ 17.290099760810172, -1.6986800597414065, 0.40857840843684384 ],
      [
        [ 0.21334698825933474, 0.6662542466119157, 0.7145546455466228 ],
        [ 0.2901637561729478, 0.653947732066262, -0.6986825876812721 ],
        [ 0.9334786124556919, -0.35486388173886213, 0.051858514489167264 ]
    ]

    np.linalg.eig(M)
    (array([17.29010017,  0.40857622, -1.69867639]),
     matrix([[ 0.21338482,  0.93292845, -0.29001971],
             [ 0.66635256, -0.35607071, -0.65512435],
             [ 0.71445166,  0.05346178,  0.69763935]]))

     */

    let E = eigens(M);
    console.log("Calculated Eigens:");
    console.log( E );

    let eigenvalues = E[0];
    let eigenvectors = E[1];

    for (let i = 0; i < eigenvalues.length; i++) {
        console.log( "Eigenvalue " + i + " =  " + eigenvalues[i] + ",  Eigenvector = " + eigenvectors[i] );
    }

}

test();