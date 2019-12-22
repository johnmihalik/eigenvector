
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

// Extract a column vector x from the matrix
function column(mat, x) {
//    let col = new Array(mat.length);
    let col = new Array(mat.length).fill(0).map(() => new Array(1).fill(0));
    for (let i = 0; i < mat.length; i ++) {
        col[i][0] = mat[i][x];
    }
    return col;
}

// Multiply 2 matrices together
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

    // If the result is a scalar dor product, just return the scalar value
    if (result.length == 1 && result[0].length == 1) {
        return result[0][0];
    }
    return result;
}

// Multiple a matrix and a scalar value
function multiply_matrix_scalar(m, v) {
    let result = [];
    for (var i = 0; i < m.length; i++) {
        result[i] = [];
        for (let j = 0; j < m[0].length; j++) {
            result[i][j] = m[i][j] * v;
        }
    }
    return result;
}
// Transpose a matrix by exchanging rows and columns
function transpose(matrix) {
    return matrix[0].map((col, i) => matrix.map(row => row[i]));
}

// Deep copy an array object
function deep_copy(a) {
    return JSON.parse(JSON.stringify(a));
}

// Sum the squared differences of 2 vectors
function squared_difference(v1, v2) {
    let sum = 0.0;
    for (let i = 0; i < v1.length; i ++) {
        sum = sum + Math.pow( v1[i] - v2[i], 2 );
    }
    return sum;
}

// Subtract matrices A - B
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


function shift_deflate(M, eigenvalue, eigenvector)  {
    let len = Math.sqrt( multiply_matrices(transpose(eigenvector), eigenvector)  );
    let U = multiply_matrix_scalar(eigenvector, 1.0/len);
    let delta = multiply_matrix_scalar( multiply_matrices(U, transpose(U)) , eigenvalue);
    let M_new = subtract_matrices(M, delta);
    return M_new;
}

function eigenvalue_of_vector(M, eigenvector) {
    // Xt * M * x
    ev = multiply_matrices( multiply_matrices(transpose(eigenvector), M ), eigenvector);
    return ev;
}

// Normalize a vector, divide by the length
function normalize_vector(v) {
    let len = Math.sqrt( multiply_matrices(transpose(v), v));
    let unit_vector = multiply_matrix_scalar(v, 1.0 / len);
    return unit_vector;
}

// Compute the most prominent eigenvector and eigenvalue
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

        let Mv = multiply_matrices(M,eigenvector);
        eigenvector = normalize_vector(Mv);
        eigenvalue = eigenvalue_of_vector(M, eigenvector);

        epsilon = Math.abs( eigenvalue - old_eigenvalue);
        iterations++;

    } while (epsilon > tolerance && iterations < max_iterations);

    return [eigenvalue, eigenvector];
}

// Compute all of the eigenvectors for a given matrix
// Matrix M MUST be a square and symmetric!
function eigens(M, tolerance=0.0001, max_iterations=1000) {

    let eigenvalues = [];
    let eigenvectors = [];

    for (let i = 0; i < M.length; i++ ) {

        // Compute the remaining most prominent eigenvector of the matrix M
        let result = power_iteration(M, tolerance, max_iterations);
        eigenvalues[i] = result[0];
        eigenvectors[i] = result[1];

        // Now remove or peel off the last eigenvector
        M = shift_deflate(M, eigenvalues[i], eigenvectors[i]);
    }

    return [eigenvalues, eigenvectors];
}

console.log("Starting Eigenvectors...");

const M = [[1,2,3], [2,7,9], [3,9,8]];
/*
np.linalg.eig(M)
(array([17.29010017,  0.40857622, -1.69867639]),
 matrix([[ 0.21338482,  0.93292845, -0.29001971],
         [ 0.66635256, -0.35607071, -0.65512435],
         [ 0.71445166,  0.05346178,  0.69763935]]))
 */

//let eig = power_iteration(mat);
let E = eigens(M);
console.log( E );

for (let i = 0; i <E[0].length; i++) {
    console.log("Eigenvalue " + i + " =  " + E[0][i]);
    console.log("Eigenvector:");
    for (let j = 0; j < E[1][0].length; j++) {
       console.log(" " + E[1][i][j]);
   }
}


/*

let X = deep_copy(mat);

console.log("Original");
pretty_print_matrix(X);

console.log("Transposed:");
let X_t = transpose(X);
pretty_print_matrix(X_t);

console.log("Vector T:");
let t = column(X,0);
pretty_print_matrix(t);

console.log("P Product:");
let p = multiply_matrices(X_t, t);
pretty_print_matrix(p);

p_length = Math.sqrt(multiply_matrices(transpose(p), p));
console.log("P Length = " + p_length);

console.log("Normalized P");
p  = multiply_matrix_scalar(p,1.0 / p_length);
pretty_print_matrix(p);
*/



