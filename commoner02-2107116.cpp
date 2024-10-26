#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

using namespace std;

// Function to perform RK4 method
void rungeKutta() {
    cout<<"Equation List"<<endl;
    cout << "1. a*sin(bx)" << endl;
    cout << "3. a*sin(bx) + a*cos(bx)" << endl;
    cout << "4. a*e^x" << endl;
    cout << "5. a*x^n + b*y^n" << endl<<endl;;
    double x0, y0, x, h, a, b, p;
    //Enter the coefficients and initial values
    cout << "Enter the coefficient a: ";
    cin >> a;
    cout << "Enter the coefficient b: ";
    cin >> b;
    cout << "Enter the initial value of x (x0): ";
    cin >> x0;
    cout << "Enter the initial value of y (y0): ";
    cin >> y0;
    cout << "Enter the value of x at which you want to find y: ";
    cin >> x;
    cout << "Enter the step size (h): ";
    cin >> h;
    //Choose the differential equation form
    cout << "Choose the differential equation form: " << endl;
    cout << "1. a*sin(bx)" << endl;
    cout << "3. a*sin(bx) + a*cos(bx)" << endl;
    cout << "4. a*e^x" << endl;
    cout << "5. a*x^n + b*y^n" << endl;
    int choice;
    cin >> choice;
    function<double(double, double)> f;
    switch (choice) {
        case 1:
            f = [a, b](double x, double y) {
                return a * sin(b * x);
            };
            break;
        case 3:
            f = [a, b](double x, double y) {
                return a * sin(b * x) + a * cos(b * x);
            };
            break;
        case 4:
            f = [a](double x, double y) {
                return a * exp(x);
            };
            break;
        case 5:
            cout << "Enter the power n: ";
            cin >> p;
            f = [a, b, p](double x, double y) {
                return a * pow(x, p) + b * pow(y, p);
            };
            break;
        default:
            cout << "Invalid choice. Using a*sin(bx) by default." << endl;
            f = [a, b](double x, double y) {
                return a * sin(b * x);
            };
            break;
    }
    int n = (int)((x - x0) / h);
    double k1, k2, k3, k4;
    double y = y0;
    cout << "x\ty" << endl;
    cout << x0 << "\t" << y0 << endl; // Initial values
    for (int i = 1; i <= n; i++) {
        k1 = h * f(x0, y);
        k2 = h * f(x0 + 0.5 * h, y + 0.5 * k1);
        k3 = h * f(x0 + 0.5 * h, y + 0.5 * k2);
        k4 = h * f(x0 + h, y + k3);
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        x0 = x0 + h;
        cout << x0 << "\t" << y << endl; // Show x and y at every step
    }
}


//Function for Gauss-Jordan elimination
void matrixInverse() {
    int n;
    cout << "Enter the size of the matrix: ";
    cin >> n;
    vector<vector<double>> matrix(n, vector<double>(n));
    vector<vector<double>> matrix2(n, vector<double>(2 * n));
    cout << "Enter the elements of the matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> matrix[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix2[i][j] = matrix[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            if (j == n + i) {
                matrix2[i][j] = 1;
            }
            else {
                matrix2[i][j] = 0;
            }
        }
    }
    cout<< "The Augumented Matrix is :" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2 * n; j++) {
            if (j == n-1) {
                cout << matrix2[i][j] << " | ";
            }
            else {
                cout << matrix2[i][j] << " ";
            }
        }
        cout << endl;
    }
    cout << endl;
    cout<<"Applying Gauss-Jordan Elimination : " << endl;
    for (int i = 0; i < n; i++) {
        double diag = matrix2[i][i];
        for (int j = 0; j < 2 * n; j++){
            matrix2[i][j] /= diag;
            if (abs(matrix2[i][j]) < 0.0001) {
                matrix2[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++){
            for (int j = 0; j < 2 * n; j++){
                if (j == n-1) {
                    cout << matrix2[i][j] << " | ";
                }
                else {
                    cout << matrix2[i][j] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = matrix2[k][i];
                for (int j = 0; j < 2 * n; j++)
                {
                    matrix2[k][j] -= factor * matrix2[i][j];
                    if (abs(matrix2[i][j]) < 0.0001)
                    {
                        matrix2[i][j] = 0;
                    }
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 2 * n; j++) {
                if (j == n-1) {
                    cout << matrix2[i][j] << " | ";
                }
                else {
                    cout << matrix2[i][j] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << "The Inverse Matrix is:\n\n";
    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            cout << matrix2[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    int c;
    cin>>c;
    if(c==3){
       rungeKutta();
    }
    if(c==4){
        matrixInverse();
    }
    return 0;
}
