#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

int main()
{
    // Open the file for reading
    std::ifstream inputFile("output_file_complex_numbers.txt");

    // Check if the file is open
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Declare variables
    double value1, realPart, imagPart;
    int value2, value3, value4, value5;
    std::string str1, str2, str3;
    bool boolValue, boolValues[8];
    // std::complex<double> complexNumbers[2];
    std::complex<double> complexNumber1;
    std::complex<double> complexNumber2;

    // Read values from the file
    inputFile >> value1 >> value2 >> value3 >> str1 >> str2 >> str3 >> value4 >> std::boolalpha >> boolValue;

    // Read complex numbers
    // for (int i = 0; i < 2; ++i)
    //{
    //    char discard;
    //    inputFile >> discard >> realPart >> discard >> imagPart >> discard;
    //    complexNumbers[i] = std::complex<double>(realPart, imagPart);
    //}
    char discard;
    inputFile >> discard >> realPart >> discard >> imagPart >> discard;
    complexNumber1 = std::complex<double>(realPart, imagPart);
    inputFile >> discard >> realPart >> discard >> imagPart >> discard;
    complexNumber2 = std::complex<double>(realPart, imagPart);

    // Read the remaining values
    inputFile >> value5;
    for (int i = 0; i < 8; ++i)
    {
        inputFile >> std::boolalpha >> boolValues[i];
    }

    // Close the file
    inputFile.close();

    // Display the read values
    std::cout << "Read values from the file:" << std::endl;
    std::cout << value1 << std::endl;
    std::cout << value2 << std::endl;
    std::cout << value3 << std::endl;
    std::cout << str1 << std::endl;
    std::cout << str2 << std::endl;
    std::cout << str3 << std::endl;
    std::cout << value4 << std::endl;
    std::cout << boolValue << std::endl;

    // Display complex numbers
    std::cout << "Complex numbers:" << std::endl;
    // for (const auto &complexNumber : complexNumbers)
    //{
    //     std::cout << "(" << complexNumber.real() << ", " << complexNumber.imag() << ")" << std::endl;
    // }

    std::cout << "(" << complexNumber1.real() << ", " << complexNumber1.imag() << ")" << std::endl;
    std::cout << "(" << complexNumber2.real() << ", " << complexNumber2.imag() << ")" << std::endl;
    std::cout << value5 << std::endl;

    // Display bool values
    std::cout << "Bool values:" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << std::boolalpha << boolValues[i] << std::endl;
    }

    return 0;
}
