#include<bits/stdc++.h>

using namespace std;

const double tol = 0.0001;
const int highest_iter = 1000;

void jacobiMethod();
void gaussSeidelMethod();
void gaussElimination();
void gaussJordanElimination();
void L_U_factorization();
void printMatrix(const vector<vector<double>>& A, const vector<double>& b);

int linear_main() {
    int choice;
    while (true) {
        cout << "\nNumerical Methods Application\n";
        cout << "1. Solution of Linear Equations\n";
        cout << "2. Exit\n";
        cout << "Choose a method : ";
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "1. Jacobi Iterative Method\n"
                     << "2. Gauss-Seidel Iterative Method\n"
                     << "3. Gauss Elimination\n"
                     << "4. Gauss-Jordan Elimination\n"
                     << "5. LU Factorization\n";
                int linearChoice;
                cin >> linearChoice;
                switch (linearChoice) {
                    case 1: jacobiMethod(); break;
                    case 2: gaussSeidelMethod(); break;
                    case 3: gaussElimination(); break;
                    case 4: gaussJordanElimination(); break;
                    case 5: L_U_factorization(); break;
                    default: cout << "Invalid choice.\n";
                }
                break;
            case 2:
                cout << "Exiting program.\n";
                return 0;
            default:
                cout << "Invalid choice.\n";
        }
    }

    return 0;
}

// Function to input the augmented matrix
void inputMatrix(int& n, vector<vector<double>>& A, vector<double>& b) {
    cout << "Input variables number : ";
    cin >> n;
    if(n<1)
    {
        cout << "Insufficient variables , check again and input ";
        return ;
    }
    A.resize(n, vector<double>(n));
    b.resize(n);

    cout << "Enter the coefficients of matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
         {
            cin >> A[i][j];
        }
    }
    cout << "Enter the constants vector b:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }
}

// Check if matrix A is diagonally dominant
bool is_diag_dominant(const vector<vector<double>>& A) {
    for (int i = 0; i < A.size(); i++) {
        double sum = 0;
        for (int j = 0; j < A.size(); j++) {
            if (i != j) {
                sum += fabs(A[i][j]);
            }
        }
        if (fabs(A[i][i]) < sum) {
            return false;
        }
    }
    return true;
}

// Attempt to make matrix A diagonally dominant by row swapping
bool make_diag_dominant(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        double max_val = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > max_val) {
                maxRow = k;
                max_val = fabs(A[k][i]);
            }
        }
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]);
        }
    }
    return is_diag_dominant(A);
}


// Jacobi method
void jacobiMethod() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;
    if (n != 5) {  // Ensure it's exactly 5 variables
        cout << "This method requires exactly 5 variables." << endl;
        return;
    }

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    cout << "Enter the coefficients of matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
    cout << "Enter the constants vector b:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    if (!is_diag_dominant(A)) {
        cout << "Warning: The matrix is not diagonally dominant. Attempting to make it diagonally dominant..." << endl;
        if (!make_diag_dominant(A, b)) {
            cout << "Error: Unable to make the matrix diagonally dominant. The Jacobi method may not converge." << endl;
            return;
        }
    }

    vector<double> x(n, 0), x_old(n, 0);
    bool converged = false;

    cout << setw(10) << "Iteration";
    for (int i = 0; i < n; i++) {
        cout << setw(12) << "x" + to_string(i + 1);
    }
    cout << endl;

    for (int iter = 0; iter < highest_iter && !converged; iter++) {
        x_old = x;
        converged = true;

        cout << setw(10) << iter + 1;
        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i][j] * x_old[j];
                }
            }

            if (fabs(A[i][i]) < 1e-9) {
                cout << "\nError: Zero or near-zero diagonal element at row " << i + 1 << "." << endl;
                return;
            }

            x[i] = sum / A[i][i]; 

            if (fabs(x[i] - x_old[i]) > tol) {
                converged = false;
            }

            cout << setw(12) << fixed << setprecision(6) << x[i];
        }
        cout << endl;

        if (converged) {
            cout << "Converged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }

    if (!converged) {
        cout << "Warning: Solution did not converge within " << highest_iter << " iterations." << endl;
    }

    cout << "\nFinal Solution: ";
    for (double xi : x) {
        cout << fixed << setprecision(6) << xi << " ";
    }
    cout << endl;
}

// Gauss-Seidel Method
void gaussSeidelMethod() {
    int n;
    vector<vector<double>> A;
    vector<double> b;
    inputMatrix(n, A, b);

    if (n != 5) {  // Ensure it's exactly 5 variables
        cout << "This method requires exactly 5 variables." << endl;
        return;
    }

    if (!is_diag_dominant(A)) {
        cout << "Warning: The matrix is not diagonally dominant. Attempting to make it diagonally dominant..." << endl;
        if (!make_diag_dominant(A, b)) {
            cout << "Error: Unable to make the matrix diagonally dominant. The Gauss-Seidel method may not converge." << endl;
            return;
        }
    }

    vector<double> x(n, 0);  
    bool converged = false;

    cout << setw(10) << "Iteration";
    for (int i = 0; i < n; i++) {
        cout << setw(12) << "x" + to_string(i + 1);
    }
    cout << endl;
    cout << string(10 + n * 12, '-') << endl;

    for (int iter = 0; iter < highest_iter && !converged; iter++) {
        converged = true;
        cout << setw(10) << iter + 1;

        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i][j] * x[j];
                }
            }

            if (fabs(A[i][i]) < 1e-9) {
                cout << "\nError: Zero or near-zero diagonal element encountered at row " << i + 1 << "." << endl;
                return;
            }

            double x_new = sum / A[i][i];
            if (fabs(x_new - x[i]) > tol) {
                converged = false;
            }
            x[i] = x_new;

            cout << setw(12) << fixed << setprecision(6) << x[i];
        }
        cout << endl;

        if (converged) {
            cout << "\nConverged after " << iter + 1 << " iterations." << endl;
            break;
        }
    }

    if (!converged) {
        cout << "\nSolution did not converge within " << highest_iter << " iterations." << endl;
    }

    cout << "\nFinal Solution: ";
    for (double xi : x) {
        cout << fixed << setprecision(6) << xi << " ";
    }
    cout << endl;
}


// Gauss Elimination Method
void gaussElimination() {
    int nn;
    vector<vector<double>> A;
    vector<double> b;
    inputMatrix(nn, A, b);

    int n = A.size();
    cout << "Initial Augmented Matrix:" << endl;
    printMatrix(A, b);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            cout << "Eliminating variable x" << i + 1 << " from row " << j + 1 << " using factor " << factor << endl;
            for (int k = i; k < n; k++)
                A[j][k] -= factor * A[i][k];
            b[j] -= factor * b[i];

            cout << "Updated Augmented Matrix after elimination from row " << j + 1 << ":" << endl;
            printMatrix(A, b);
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }

    cout << "Solution: ";
    for (double xi : x) cout << xi << " ";
    cout << endl;
}

// Gauss-Jordan Elimination Method
void gaussJordanElimination() {
    int nn;
    vector<vector<double>> A;
    vector<double> b;
    inputMatrix(nn, A, b);

    int n = A.size();
    cout << "Initial Augmented Matrix:" << endl;
    printMatrix(A, b);

    for (int i = 0; i < n; i++) {
        double pivot = A[i][i];
        cout << "Normalizing row " << i + 1 << " with pivot " << pivot << endl;
        for (int j = 0; j < n; j++)
            A[i][j] /= pivot;
        b[i] /= pivot;

        cout << "Updated Augmented Matrix after normalization:" << endl;
        printMatrix(A, b);

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double factor = A[j][i];
                cout << "Eliminating variable x" << i + 1 << " from row " << j + 1 << " using factor " << factor << endl;
                for (int k = 0; k < n; k++)
                    A[j][k] -= factor * A[i][k];
                b[j] -= factor * b[i];

                cout << "Updated Augmented Matrix after elimination from row " << j + 1 << ":" << endl;
                printMatrix(A, b);
            }
        }
    }

    cout << "Solution: ";
    for (double bi : b) cout << bi << " ";
    cout << endl;
}

// Function to perform LU decomposition of matrix A
void LU_Decomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int n) {
    // Decompose matrix into L and U
    for (int i = 0; i < n; i++) {
        // Check for zero pivot
        if (fabs(A[i][i]) < 1e-9) {
            // Try to find a non-zero pivot by row swapping
            bool swapped = false;
            for (int k = i + 1; k < n; k++) {
                if (fabs(A[k][i]) > 1e-9) {
                    swap(A[i], A[k]);
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                cout << "Error: Matrix is singular; LU Factorization cannot be performed." << endl;
                return;
            }
        }

        // Upper triangular matrix U
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        // Lower triangular matrix L
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal elements of L are 1
            else {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Function to perform forward substitution
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b, int n) {
    vector<double> y(n, 0);
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
    return y;
}

// Function to perform backward substitution
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y, int n) {
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}
void L_U_factorization(){
    int n;

    cout << "Enter the number of equations: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "Enter the coefficients of the matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the constants of the equation (vector b):" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    // Perform LU decomposition
    LU_Decomposition(A, L, U, n);

    // Print L and U matrices
    cout << "Lower Triangular Matrix L:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(3) << L[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nUpper Triangular Matrix U:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(3) << U[i][j] << " ";
        }
        cout << endl;
    }

    // Solve Ly = b using forward substitution
    vector<double> y = forwardSubstitution(L, b, n);

    // Solve Ux = y using backward substitution
    vector<double> x = backwardSubstitution(U, y, n);

    // Print the solution
    cout << "\nSolution (x):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
    }

}
// Function to print matrix
void printMatrix(const vector<vector<double>>& A, const vector<double>& b) {
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++)
            cout << fixed << setprecision(2) << A[i][j] << "\t";
        cout << "| " << fixed << setprecision(2) << b[i] << endl;
    }
}
