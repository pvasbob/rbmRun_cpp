#include <iostream>
#include <fstream>
#include <complex>

int main()
{
    // Open the file for reading
    std::ifstream inputFile("complex_output.txt");

    // Check if the file is open
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Variables to store real and imaginary parts
    double realPart, imagPart;

    // Read the complex number from the file
    char leftParenthesis, comma;
    inputFile >> leftParenthesis >> realPart >> comma >> imagPart;

    // Check if the read was successful
    if (inputFile.fail())
    {
        std::cerr << "Error reading complex number from file!" << std::endl;
        return 1;
    }

    // Create a complex number using the read parts
    std::complex<double> readComplexNumber(realPart, imagPart);

    // Close the file
    inputFile.close();

    // Display the read complex number
    std::cout << "Read complex number from file: " << readComplexNumber << std::endl;

    return 0;
}
