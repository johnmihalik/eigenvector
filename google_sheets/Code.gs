/**
 * A self contained implementation of eigen-decomposition using the power iteration method and
 * shifting by deflation.  This implementation uses a simple array based matrix implementation
 * to allow this code to be embedded without external dependencies.
 *
 *   See:  https://en.wikipedia.org/wiki/Power_iteration
 */

/**
 * This is the button callback to invoke the eigens calculations.
 * This is a workaround due to the fact that custom functions must have deterministic inputs.
 * This prevents us passing references to functions like NOW() and Google Finance queries
 * as inputs to a custom function.
 *
 * This simply copies the values of the covariance matrix to a target range as
 * a copy paste values.  Due to the complexity of the cell functions, it appears
 * that Google ApScript could not interpret or keep up with the cell calculations.
 */
function copyMatrixValues() {
  
  var ss = SpreadsheetApp.getActiveSpreadsheet();
  
  // Copy the values only 
  ss.getRange('Portfolio!C23:L32').copyTo(ss.getRange('PCA!C101:C111'), SpreadsheetApp.CopyPasteType.PASTE_VALUES, false);
  
}


/**
 *  Compute the eigenvalue and eigenvector decomposition of a cell range.  
 *  The cell range must be square and symetric.
 *  
 *  @return {Array[]} A 2 dimensional array.  The first row are the eigenvalues followed by colummns of the cooresponding eigenvectors.
 *
 *  @customfunction
 */
function eigens(input) {

  try {

    var M = filterMatrixBlanks(input);
    
    var E = __eigens(M, 0.001, 100);

    var eigenvalues = E[0];
    var eigenvectors = E[1];

    // Format the output as a M+1 x M multi-dimensional array
    var out = [];

    // The eigenvalues form the first row.
    out[0] = E[0];

    // The subsequent columns under each eigenvalue are the eigenvectors
    eigenvectors = transpose(eigenvectors);
    for (var i = 0; i < eigenvectors.length; i++ ) {
        out[i+1] = eigenvectors[i];
    }

    return out;
    
  } catch( ex ) {
     return "Error: " + ex;
  }
  
}


/**
 *  Compute the eigenvalues of a cell range as input.
 *  The cell range must be square and symetric.
 *  
 *  @return {Array} The eigenvalues of the matrix.
 *
 *  @customfunction
 */
function eigenvalues(input) {

  try {
    
    var E = eigens(input)
    return E[0];
    
  } catch( ex ) {
     return "Error: " + ex;
  }
  
}


/**
 *  Compute the eigenvectors of a cell range as imput.  
 *  The cell range must be square and symetric.
 *  
 *  @return {Array[]} The eigenvectors of the matrix as column vectors.
 *
 *  @customfunction
 */
function eigenvectors(input) {

  try {
    
    var E = eigens(input)
    return E.slice(1);
    
  } catch( ex ) {
     return "Error: " + ex;
  }
  
}



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
function power_iteration(M, tolerance, max_iterations) {

    var rank = M.length;

    // Initialize the first guess pf the eigenvector to a row vector of the sqrt of the rank
    //var eigenvector = new Array(rank).fill(0).map( function() {new Array(1).fill(Math.sqrt(rank)) });
    var eigenvector = init_array(rank,1,Math.sqrt(rank));

    // Compute the corresponding eigenvalue
    var eigenvalue = eigenvalue_of_vector(M, eigenvector);

    var epsilon = 1.0;
    var iterations = 0;
    do {

        var old_eigenvalue = deep_copy(eigenvalue);

        // Multiply the Matrix M by the guessed eigenveector
        var Mv = multiply_matrices(M,eigenvector);

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
function __eigens(M, tolerance, max_iterations) {

    var eigenvalues = [];
    var eigenvectors = [];

    for (var i = 0; i < M.length; i++ ) {

        // Compute the remaining most prominent eigenvector of the matrix M
        var result = power_iteration(M, tolerance, max_iterations);

        // Separate the eigenvalue and vector from the return array
        var eigenvalue = result[0];
        var eigenvector = result[1];

        eigenvalues[i] = eigenvalue;

        // Compress the eigenvector to a simple array for eas of manipulation
        eigenvectors[i] = flatten_vector(eigenvector);

        // Now remove or peel off the last eigenvector
        M = shift_deflate(M, eigenvalue, eigenvector);
    }

    return [eigenvalues, eigenvectors];
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
    var len = Math.sqrt( multiply_matrices(transpose(eigenvector), eigenvector)  );
    var U = multiply_matrix_scalar(eigenvector, 1.0/len);
    var delta = multiply_matrix_scalar( multiply_matrices(U, transpose(U)) , eigenvalue);
    var M_new = subtract_matrices(M, delta);
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


/*
 * Matrix Functions
 *
 * Created this implementation for Google ApScript compatibility 
 * and to have a fully encapsultated implementation with no other
 * dependencies.
 */


/**
 * Filters a M x N set of input cell ranges of values to eliminate blank rows and columns.
 * Creates a matrix of size of rows x columns to be input into the eigen functions.
 * This is neccesary since the input of stocks is variable from 1 to 10.  The input
 * size of the range is 10 x 10 cells fixed and this method reduces that to eliminate 
 * possible blanks.
 */
function filterMatrixBlanks(M) {
    
    var filtered = [];
    for( var i = 0; i < M.length; i++){
        for( var j = 0; j < M[i].length; j++ ) {
          if (typeof(M[i][j]) == "number") {
              if ( j == 0 ) {
                 filtered[i] = [];
              }
              filtered[i][j] = M[i][j];
          }
        }
     }

    //return "Filtered rows = " + filtered.length + ", cols = " + filtered[0].length;   //filtered;
    return filtered;

}


/**
 *  Helper function to initialize arrays without using Arrays.fill due to Google AppScript compatibility.
 *
 */
function init_array(rows, cols, value) {

    result = [];

    //init the grid matrix
    for ( var i = 0; i < rows; i++ ) {
        result[i] = [];
        for (var j = 0; j < cols; j++) {
            result[i][j] = value;
        }
    }
    return result;
}



/**
 * Print a matrix nicely for the console.
 *
 * @param  {Array} matrix  The m x n matrix to output.  Matrices are simple 2 dimensional arrays.
 */
// Print a matrix formatted nicely
function pretty_print_matrix(mat) {
    for ( var i = 0; i < mat.length; i++) {
        process.stdout.write("[ ");
        for ( var j = 0; j < mat[i].length; j ++) {
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
    //var col = new Array(mat.length).fill(0).map( function() { new Array(1).fill(0) });
    var col = init_array(mat.length, 1, 1    );

    for (var i = 0; i < mat.length; i ++) {
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
    var result = [];
    for (var i = 0; i < m.length; i++) {
        result[i] = [];
        for (var j = 0; j < m[0].length; j++) {
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
//    return matrix[0].map(  function(col, i) { matrix.map(function(row) { row[i] } ) });
//    return matrix[0].map((col, i) => matrix.map(row => row[i]));
    var newArray = [];
    for(var i = 0; i < matrix[0].length; i++){
        newArray.push([]);
    };

    for(var i = 0; i < matrix.length; i++){
        for(var j = 0; j < matrix[0].length; j++){
            newArray[j].push(matrix[i][j]);
        };
    };

    return newArray;
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
    var sum = 0.0;
    for (var i = 0; i < v1.length; i ++) {
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
    var result = [];
    for (var i = 0; i < A.length; i++) {
        result[i] = [];
        for (var j = 0; j < A[0].length; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

/**
 *  Utility function to simplify a vector into a standard array
 *
 * @param  {[][]} v1   An m x 1 vector
 * @return {Array}    A simple array of the vector values
 */
function flatten_vector(v) {
    var v_new = [];
    for (var i = 0; i < v.length; i++) {
        v_new[i] = v[i][0];
    }
    return v_new;
}

/**
 *  Normalize a vector, divide by the length
 *
 * @param v
 * @return {*[][]}
 */
function normalize_vector(v) {
    var len = Math.sqrt( multiply_matrices(transpose(v), v));
    var unit_vector = multiply_matrix_scalar(v, 1.0 / len);
    return unit_vector;
}
