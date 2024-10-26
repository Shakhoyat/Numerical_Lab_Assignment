#include <bits/stdc++.h>
#include "2107014.cpp"
#include "commoner02-2107116.cpp"
#include "2107104.cpp"
using namespace std;



int main()
{
    int choice;
    while (true) {
        cout << "\nNumerical Methods Application\n";
        cout << "1. Solution of Linear Equations\n";
        cout << "2. Solution of Non-linear Equations\n";
        cout << "3. Solution of Differential Equations\n";
        cout << "4. Matrix Inversion\n";
        cout << "5. Exit\n";
        cout << "Choose a method: ";
        cin >> choice;

        switch (choice) {
            case 1:
                linear_main();
                break;
            case 2:
              Non_lin_main();
                break;
            case 3:
               shuvo_main(3);
                break;
            case 4:
                shuvo_main(4);
                break;
            case 5:
                cout << "Exiting program.\n";
                return 0;
            default:
                cout << "Invalid choice.\n";
        }
    }
}
