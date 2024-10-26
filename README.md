#  Collaborators

### Al Mubtasim Preom ([nashokota](https://www.github.com/nashokota))- 2107014

### Shakhoyat Rahman Shujon ([Shakhoyat](https://www.github.com/Shakhoyat))- 2107104

### Sree. Shuvo Kumar Joy ([commoner02](https://www.github.com/commoner02)) - 2107116



<br/>


# Numerical Methods Console Application

## Overview

Welcome to the Numerical Methods Console Application! This project provides various numerical methods for solving mathematical problems, making it a versatile tool for anyone interested in numerical analysis.

## Table of Contents

[Application Structure](#application-structure)

- [1. Solution of Linear Equations](#solution-of-linear-equations)
  - [a. Jacobi Iterative Method](#jacobi-iterative-method)
  - [b. Gauss-Seidel Iterative Method](#gauss-seidel-iterative-method)
  - [c. Gauss Elimination](#gauss-elimination)
  - [d. Gauss-Jordan Elimination](#gauss-jordan-elimination)
  - [e. LU Factorization](#lu-factorization)
- [2. Solution of Non-linear Equations](#solution-of-non-linear-equations)
  - [a. Bisection Method](#bisection-method)
  - [b. False Position Method](#false-position-method)
  - [c. Secant Method](#secant-method)
  - [d. Newton-Raphson Method](#newton-raphson-method)
- [3. Solution of Differential Equations](#solution-of-differential-equations)
  - [a. Runge-Kutta Method](#runge-kutta-method)
- [4. Matrix Inversion](#matrix-inversion)

---

## Application Structure

### Key Features

1. **Solve Linear Equations**

   - Several techniques to handle systems of linear equations:
     - **Jacobi Method**: An iterative approach where we make guesses until the variables converge to a solution.
     - **Gauss-Seidel Method**: Similar to Jacobi but it updates the variables immediately for faster convergence than Jacobi Method.
     - **Gauss Elimination**: Transforms the system into an upper triangular form for solving.
     - **Gauss-Jordan Elimination**: A more modified version of Gauss Elimination that directly reveals the solutions.
     - **LU Factorization**: Decomposes a matrix into two simpler matrices for efficient computation.

2. **Solve Non-Linear Equations**

   - Methods designed to find roots of non-linear equations:
     - **Bisection Method**: Narrows down an interval to find a root by checking the function’s values at both ends.
     - **False Position Method**: Similar to Bisection but employs linear interpolation for quicker convergence.
     - **Secant Method**: Utilizes two initial guesses to find the root without needing derivatives.
     - **Newton-Raphson Method**: Uses a function’s derivative to accelerate root finding.

3. **Solve Differential Equations**

   - **Runge-Kutta Method (RK4)**: A step-by-step approach for approximating solutions to ordinary differential equations.

4. **Matrix Inversion**
   - **Gauss-Jordan Elimination**: Finds the inverse of a matrix by transforming it into a simpler form.

### Solution of Linear Equations

---

### Jacobi Iterative Method

---

### Description

The Jacobi method is an iterative algorithm used to solve a system of linear equations of the form **Ax = b**. It is particularly useful for large systems and can converge particularly when the coefficient matrix \( A \) is diagonally dominant.

### Algorithm Steps

1. **Initialization**:
   - Start with an initial guess for the solution vector \( x \) (commonly \( x = 0 \)).
2. **Iteration**:

   - For each variable x(i) update x(i)^(k+1) from the k-th iteration using formula : x(i) = 1/a(ii)*{b(i)- sumof(a(i,j)*x(j)^k)} where i!=j.
   - Repeat the process until convergence (i.e., until the change in \( x \) is below a specified tolerance).

3. **Convergence Check**:
   - After each iteration, the difference between the current and previous estimates of \( x \) is checked to be less than the defined tolerance value.

### Input

- **Matrix \( A \)**: Coefficient matrix (must be square and diagonally dominant for guaranteed convergence).
- **Vector \( b \)**: Constants vector (the same length as the number of rows in \( A \)).
- **max_iterations**: Maximum number of iterations (default is 1000).
- **tolerance**: Convergence tolerance is 0.0001 .

### Output

- The method outputs the solution vector \( x \) and the number of iterations taken to converge.

### Edge Cases

- The algorithm may fail to converge if:
  - The matrix \( A \) is not diagonally dominant.
  - The maximum number of iterations is reached without convergence.
- It is advisable to check the diagonal dominance of \( A \) before using the Jacobi method.

---

### Gauss-Seidel Iterative Method

---

## Overview

The Gauss-Seidel method is an iterative technique used for solving a system of linear equations of the form **Ax = b**. It improves upon the Jacobi method by using the most recently updated values of the solution vector in subsequent calculations, which can lead to faster co

## Algorithm Steps

1. **Input**:
   - Read the coefficient matrix \( A \) and the constants vector \( b \).
2. **Diagonal Dominance Check**:

   - Check if the matrix \( A \) is diagonally dominant. If it is not, attempt to modify it to ensure diagonal dominance. If the modification fails, warn the user that the method may not converge.

3. **Initialization**:

   - Initialize the solution vector \( x \) with zeros.

4. **Iteration**:

   - For a maximum number of iterations:
     - Set a convergence flag to true.
     - For each variable \( x_i \):
       - Calculate the sum:
         $$\text{sum} = b[i] - \sum\_{j \neq i} A[i][j] \cdot x[j]$$
       - Check for zero or near-zero diagonal elements and handle appropriately.
       - Update \($x[i]$ \):
         $$x[i] = \frac{\text{sum}}{A[i][i]}$$
       - Check if the change in \( $x[i]$ \) is above the tolerance level. If so, set the convergence flag to false.
     - If all updates are within the tolerance, exit the loop.

5. **Convergence Check**:
   - If the method converges, display the number of iterations taken and the final solution vector. If it does not converge, notify the user.

## Input

- **Matrix \( A \)**: Coefficient matrix (must be square and ideally diagonally dominant).
- **Vector \( b \)**: Constants vector (the same length as the number of rows in \( A \)).
- **highest_iter**: Maximum number of iterations (default is 1000).
- **tol**: Convergence tolerance 0.0001 .

## Output

- The method outputs the solution vector \( x \) and the number of iterations taken to converge.

## Edge Cases

- The algorithm may fail to converge if:
  - The matrix \( A \) is not diagonally dominant and cannot be modified to achieve dominance.
  - The maximum number of iterations is reached without convergence.
- A warning is issued if a zero or near-zero diagonal element is encountered during computation.

---

### Gauss Elimination

---

### Description

The Gauss Elimination method is a direct algorithm used to solve systems of linear equations. It transforms the system into an upper triangular form and then performs back substitution to find the solution. This method is effective for small to moderate-sized systems.

### Algorithm Steps

1. **Input**:

   - Read the augmented matrix \( A | b \) where \( A \) is the coefficient matrix and \( b \) is the constants vector.

2. **Forward Elimination**:

   - For each variable \( x_i \):
     - Eliminate the variables below \($x_i$\) in the column.
     - Calculate the elimination factor for each row below the pivot row.
     - Update the rows of the matrix \( A \) and vector \( b \).

3. **Back Substitution**:
   - Starting from the last variable, compute the values of \($x$\) by substituting back into the equations.

### Input

- **Matrix \( A \)**: Coefficient matrix (must be square).
- **Vector \( b \)**: Constants vector (the same length as the number of rows in \( A \)).

### Output

- The method outputs the solution vector \( x \) derived from the back substitution.

---

### Gauss-Jordan Elimination

---

### Description

The Gauss-Jordan elimination method is an extension of the Gauss elimination method. It reduces the augmented matrix to reduced row echelon form (RREF). This method directly provides the solutions without the need for back substitution.

### Algorithm Steps

1. **Input**:

   - Read the augmented matrix \($A | b$\).

2. **Normalization**:
   - For each variable \( $x_i$ \):
     - Normalize the pivot row to make the pivot element equal to 1.
3. **Elimination**:
   - Eliminate all other entries in the pivot column by subtracting appropriate multiples of the pivot row from the other rows.

### Input

- **Matrix \( A \)**: Coefficient matrix (must be square).
- **Vector \( b \)**: Constants vector (the same length as the number of rows in \( A \)).

### Output

- The method outputs the solution vector \( b \) which now contains the values of the variables after reaching RREF.

---

### LU Factorization

---

## Overview

This README provides a detailed explanation of the LU Decomposition algorithm and the functions implemented for solving systems of linear equations. The LU Decomposition method is used to factor a matrix \( A \) into two matrices \( L \) (lower triangular) and \( U \) (upper triangular), allowing for efficient solving of the equation \( Ax = b \).

## LU Decomposition Algorithm

LU Decomposition decomposes a given square matrix \( A \) into the product of a lower triangular matrix \( L \) and an upper triangular matrix \( U \):

\[A = LU\]

### Steps of the Algorithm

1. **Pivoting**:
   - Check for zero pivots and perform row swapping if necessary to ensure stability.
2. **Decomposing the Matrix**:
   - The upper triangular matrix \( U \) is filled by subtracting the contributions from the previously computed values in \( L \) and \( U \).
   - The lower triangular matrix \( L \) is filled with the appropriate coefficients, ensuring the diagonal elements are set to 1.

### Functions

#### 1. LU_Decomposition

- **Purpose**: To perform LU decomposition of matrix \( A \).
- **Parameters**:
  - `vector<vector<double>>& A`: The input matrix \( A \) to be decomposed.
  - `vector<vector<double>>& L`: The output lower triangular matrix.
  - `vector<vector<double>>& U`: The output upper triangular matrix.
  - `int n`: The number of equations (size of the matrix).
- **Returns**: None. The function modifies \( L \) and \( U \) in place.

#### 2. `forwardSubstitution`

- **Purpose**: To solve the equation \( Ly = b \) using forward substitution.
- **Parameters**:
  - `const vector<vector<double>>& L`: The lower triangular matrix.
  - `const vector<double>& b`: The constants vector.
  - `int n`: The number of equations.
- **Returns**: The solution vector \( y \).

#### 3. `backwardSubstitution`

- **Purpose**: To solve the equation \( Ux = y \) using backward substitution.
- **Parameters**:
  - `const vector<vector<double>>& U`: The upper triangular matrix.
  - `const vector<double>& y`: The solution vector from the forward substitution.
  - `int n`: The number of equations.
- **Returns**: The solution vector \( x \).

#### 4. `L_U_factorization`

- **Purpose**: To handle user input, perform LU decomposition, and print the results.
- **Parameters**: None.
- **Returns**: None. This function manages the entire workflow from input to output.

#### 5. `printMatrix`

- **Purpose**: To print the augmented matrix \( [A | b] \).
- **Parameters**:
  - `const vector<vector<double>>& A`: The coefficient matrix.
  - `const vector<double>& b`: The constants vector.
- **Returns**: None. This function outputs the augmented matrix to the console.

### Solution of Non-linear Equations

---

### Bisection Method

---

This method finds a root of the function by repeatedly narrowing down an interval where a sign change occurs.

1. **Check Interval**: Ensure that the function values at the endpoints `f(a)` and `f(b)` have opposite signs (`f(a) * f(b) < 0`). This indicates that a root lies between `a` and `b`.
2. **Calculate Midpoint**: Compute the midpoint of the interval:
   $$\text{midpoint} = \frac{a + b}{2}$$
3. **Check Midpoint**: If $f(\text{midpoint})$ is zero, you've found the root!
4. **Narrow the Interval**: Determine which half of the interval contains the root:
   - If \( $f(a)$ \) and \($f(\text{midpoint})$ \) have opposite signs, set \($b = \text{midpoint}$\).
   - Otherwise, set \($a = \text{midpoint}$\).
5. **Repeat**: Continue this process until the interval is sufficiently small or the maximum number of iterations is reached

---

### False Position Method

---

This method is a faster alternative to the Bisection Method, using linear interpolation to improve convergence speed.

1. **Check Interval**: Verify that `f(a)` and `f(b)` have opposite signs.
2. **Linear Interpolation**: Calculate the new guess for the root:
   $$\text{guess} = \frac{a \cdot f(b) - b \cdot f(a)}{f(b) - f(a)}$$
3. **Update Interval**: Check where the sign change occurs:
   - If $f(\text{guess})$ has the same sign as $f(a)$, set a = guess
   - Otherwise, set $b = \text{guess}$.
4. **Convergence Check**: Stop if $|f(\text{guess})| < kTolerance$.

---

### Secant Method

---

The Secant Method approximates roots using secant lines instead of derivatives.

1. **Initial Guesses**: Start with two initial values, $x_0$ and $x_1$.
2. **Calculate Next Approximation**:
   $$x_2 = x_1 - \frac{f(x_1) \cdot (x_1 - x_0)}{f(x_1) - f(x_0)}$$
3. **Update**: Move \($x_0$\) to \($x_1$\) and \($x_1$\) to \($x_2$\).
4. **Convergence Check**: Stop if \($|x_2 - x_1| < kTolerance$\).

---

### Newton-Raphson Method

---

This method is powerful for finding roots but requires the derivative of the function.

1. **Derivative Requirement**: Ensure you have the function's derivative $f'(x)$
2. **Iteration Formula**:
   $$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$
3. **Update**: Set
   $$x_n = x_{n+1}$$
4. **Convergence Check**: Stop if $$|x\_{n+1} - x_n| < kTolerance$$

### Special Case handling for Non-Linear Equations

- **Bisection and False Position Methods**: If $$f(a) \cdot f(b) \geq 0$$ the program outputs a message stating no root exists in that interval and terminates.
- **Secant Method**: Ensures that the initial guesses $x_0$ and $x_1$ are not the same to prevent errors.
- **Newton-Raphson Method**: Checks if $$f'(x_n)$$ is too close to zero to avoid division errors.

### Solution of Differential Equations

---

### Runge-Kutta Method

---

Here is implementation of the classic fourth-order Runge-Kutta method (RK4) for solving ordinary differential equations. This method uses a weighted average of four increments to estimate the solution with higher accuracy.

The main steps of the RK4 method are:

**a. Calculate the initial conditions**: Determine the starting values for \( x \) and \( y \).<br/>
**b. Estimate the four increments (slopes)**:

- \( k1 \): Calculate the slope at the initial point.
- \( k2 \): Calculate the slope at the midpoint using \( k1 \).
- \( k3 \): Calculate another midpoint slope using \( k2 \).
- \( k4 \): Calculate the slope at the end of the interval using \( k3 \).

**c. Combine the slopes**: Use a weighted average of these four slopes to estimate the value of \( y \) at the next step.<br/>
**d. Repeat the process**: Move to the next interval and repeat the calculations to trace out the entire path of the solution.

This iterative process provides a highly accurate approximation of the solution to the differential equation.

---

### Matrix Inversion

---

For matrix inversion, the implemententation was done by the Gauss-Jordan Elimination method, which is straightforward and effective.

Here's the process:

**a. Forming an Augmented Matrix**: Start by placing the original matrix side-by-side with an identity matrix.<br/>
**b. Applying Gauss-Jordan Elimination**: This step involves transforming the original matrix into an identity matrix through a series of row operations.<br/>
**c. Simultaneous Transformation**: While the original matrix is being transformed into an identity matrix, the identity matrix undergoes the same row operations, gradually turning into the inverse of the original matrix.<br/>

While this method has a time complexity of \( O(n^3) \) (with \( n \) being the order of the square matrix), it is intuitive and highly effective for many practical applications, even if it's not the fastest option for very large matrices.

### So, We can Conclude that

We organized our application in a suitable manner as we preferred. For the linear,non-linear and differential equations section.Our implementation of code is working fine in solving a Matrix Inversion in step by step manner ;)
