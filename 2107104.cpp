#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
using namespace std;

// Tolerance for convergence and maximum number of iterations limit
const double kTolerance = 1e-6;
const int kMaxIterations = 1000;

// Function prototypes for different non-linear methods
void ApplyBisection(function<double(double)> eqn, double start, double end);
void ApplyFalsePosition(function<double(double)> eqn, double start, double end);
void ApplySecant(function<double(double)> eqn, double x0, double x1);
void ApplyNewtonRaphson(function<double(double)> eqn, function<double(double)> derivative, double x0);

// Evaluating the polynomial at a given point
double EvaluatePolynomial(const vector<double>& coefficients, double x) {
    double sum = 0;
    int poly_degree = coefficients.size() - 1;
    for (int i = 0; i <= poly_degree; i++) {
        sum += coefficients[i] * pow(x, poly_degree - i);
    }
    return sum;
}

// Evaluating the d/dx of the polynomial at a given point
double EvaluateDerivative(const vector<double>& coefficients, double x) {
    double sum = 0;
    int poly_degree = coefficients.size() - 1;
    for (int i = 0; i < poly_degree; i++) {
        sum += (poly_degree - i) * coefficients[i] * pow(x, poly_degree - i - 1);
    }
    return sum;
}

int Non_lin_main() {
    int poly_degree;
    vector<double> coefficients;

    int choice;
    while (true) {
        cout << "\nNumerical Methods Application\n";
        cout << "1. Solution of Non-Linear eqns\n";
        cout << "2. Exit\n";
        cout << "Choose an option: ";
        cin >> choice;

        if (choice == 2) {
            cout << "Exiting program.\n";
            break;
        } else if (choice != 1) {
            cout << "Invalid choice.\n";
            continue;
        }
 cout << "\n1. Bisection Method\n"
             << "2. False Position Method\n"
             << "3. Secant Method\n"
             << "4. Newton-Raphson Method\n";
        cout << "Select a root-finding method: ";
        
        int nonlinearChoice;
        cin >> nonlinearChoice;
        cout << "Enter the polynomial degree : ";
        cin >> poly_degree;

        if (poly_degree < 1 || poly_degree > 5) {
            cout << "Degree should be between 1 and 5.\n";
            continue;
        }

        coefficients.resize(poly_degree + 1);
        cout << "Enter the polynomial coefficients from highest to lowest degree:\n";
        for (int i = 0; i <= poly_degree; i++) {
            cout << "Coefficient for x^" << (poly_degree - i) << ": ";
            cin >> coefficients[i];
        }

        // Lambda function for the polynomial and its derivative
        auto eqn = [coefficients](double x) {
            return EvaluatePolynomial(coefficients, x);
        };
        auto derivative = [coefficients](double x) {
            return EvaluateDerivative(coefficients, x);
        };

       

        double start, end, x0, x1;
        switch (nonlinearChoice) {
            case 1:
                cout << "Enter the interval [start, end] for Bisection method:\n";
                cout << "start: "; cin >> start;
                cout << "end: "; cin >> end;
                ApplyBisection(eqn, start, end);
                break;
            case 2:
                cout << "Enter the interval [start, end] for False Position method:\n";
                cout << "start: "; cin >> start;
                cout << "end: "; cin >> end;
                ApplyFalsePosition(eqn, start, end);
                break;
            case 3:
                cout << "Enter initial guesses x0 and x1 for Secant method:\n";
                cout << "x0: "; cin >> x0;
                cout << "x1: "; cin >> x1;
                ApplySecant(eqn, x0, x1);
                break;
            case 4:
                cout << "Enter initial guess x0 for Newton-Raphson method:\n";
                cout << "x0: "; cin >> x0;
                ApplyNewtonRaphson(eqn, derivative, x0);
                break;
            default:
                cout << "Invalid choice.\n";
        }
    }
    return 0;
}

//  Implementation of Bisection Method
void ApplyBisection(function<double(double)> eqn, double start, double end) {
    if (eqn(start) * eqn(end) >= 0) {
        cout << "No root in this interval.\n";
        return;
    }

    double midpoint = start;
    for (int i = 0; i < kMaxIterations; i++) {
        midpoint = (start + end) / 2;
        if (eqn(midpoint) == 0 || fabs(end - start) / 2 < kTolerance) {
            cout << "Root: " << midpoint << ", Iterations: " << i + 1 << endl;
            return;
        }
        (eqn(midpoint) * eqn(start) < 0) ? end = midpoint : start = midpoint;
    }
    cout << "Approximate root after " << kMaxIterations << " iterations: " << midpoint << endl;
}

// False Position Method Implementation
void ApplyFalsePosition(function<double(double)> eqn, double start, double end) {
    if (eqn(start) * eqn(end) >= 0) {
        cout << "No root in this interval.\n";
        return;
    }

    double guess = start;
    for (int i = 0; i < kMaxIterations; i++) {
        double fa = eqn(start), fb = eqn(end);
        guess = (start * fb - end * fa) / (fb - fa);
        if (eqn(guess) == 0 || fabs(eqn(guess)) < kTolerance) {
            cout << "Root: " << guess << ", Iterations: " << i + 1 << endl;
            return;
        }
        (eqn(guess) * fa < 0) ? end = guess : start = guess;
    }
    cout << "Approximate root after " << kMaxIterations << " iterations: " << guess << endl;
}

// Secant Method Implementation
void ApplySecant(function<double(double)> eqn, double x0, double x1) {
    double x2;
    for (int i = 0; i < kMaxIterations; i++) {
        double f0 = eqn(x0), f1 = eqn(x1);
        if (fabs(f1 - f0) < kTolerance) {
            cout << "Close to division by zero. No solution found.\n";
            return;
        }
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        if (fabs(x2 - x1) < kTolerance) {
            cout << "Root: " << x2 << ", Iterations: " << i + 1 << endl;
            return;
        }
        x0 = x1;
        x1 = x2;
    }
    cout << "Approximate root after " << kMaxIterations << " iterations: " << x2 << endl;
}

// Newton-Raphson Method Implementation
void ApplyNewtonRaphson(function<double(double)> eqn, function<double(double)> derivative, double x0) {
    double x1;
    for (int i = 0; i < kMaxIterations; i++) {
        double df0 = derivative(x0);
        if (fabs(df0) < kTolerance) {
            cout << "Derivative near zero. No solution found.\n";
            return;
        }
        x1 = x0 - eqn(x0) / df0;
        if (fabs(x1 - x0) < kTolerance) {
            cout << "Root: " << x1 << ", Iterations: " << i + 1 << endl;
            return;
        }
        x0 = x1;
    }
    cout << "Approximate root after " << kMaxIterations << " iterations: " << x1 << endl;
}
