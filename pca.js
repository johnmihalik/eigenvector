
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
    return result;
}

// Multiple a matrix and a scalar value
function multiply_matrix_scalar(m, v) {
    var result = [];
    for (var i = 0; i < m.length; i++) {
        result[i] = [];
        for (var j = 0; j < m[0].length; j++) {
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


function eigen(mat) {

    let X = deep_copy(mat);
    pretty_print_matrix(X);

    const rank = X.length;
    console.log(" Matrix Rank = " + rank);

    // Set the tolerance for convergence
    let tollerance = 0.00001;

    // Initialize the principal component vectors, each PC is an array
    //let t = new Array(rank).fill(0).map(() => new Array(rank).fill(0));

    let t = new Array(rank);

    let p = new Array(rank);

    let X_t = transpose(X);

    // Set the first principal component to the first column vector of the input matrix
    t[0] = column(X, 0);

    let epsilon = 1.0;

    let iterations = 0;
    do {

        p[0] = multiply_matrices(X_t, t[0]);
        let tp = multiply_matrices(transpose(t[0]), t[0]);
        p[0] = multiply_matrix_scalar(p[0], 1.0 / tp);

        // Normalize p
        let p_length = Math.sqrt(multiply_matrices(transpose(p[0]), p[0]));
        p[0] = multiply_matrix_scalar(p[0], 1.0 / p_length);

        let t_new = multiply_matrices(X, p[0]);
        let pp = multiply_matrices(transpose(p[0]), p[0]);
        t_new = multiply_matrix_scalar(t_new, 1.0 / pp);

        epsilon = squared_difference(t[0], t_new);

        t[0] = deep_copy(t_new);

        iterations++;
    } while ( epsilon > tollerance);

    console.log("Iterations = " + iterations);

    console.log("T: ")
    pretty_print_matrix(t[0]);

    console.log("P");
    pretty_print_matrix(p[0]);

    console.log("PC?");
    let pc = multiply_matrices(transpose(t[0]), t[0]);
    pretty_print_matrix(pc)

}

console.log("Starting Eigenvectors...");

// example:
// m = numpy.mat([[4, -6, 5], [-6, 3, 4], [5, 4, -3]])
// numpy.linalg.eigvals(m) #=> array([-9.12030391,  9.62192181,  3.4983821 ])
const mat = [
    [4, -6, 5],
    [-6, 3, 4],
    [5, 4, -3]
];

eigen(mat);

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



