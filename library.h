#ifndef NUMERICALMETHODSLIB_LIBRARY_H
#define NUMERICALMETHODSLIB_LIBRARY_H

#include "iostream"
#include "vector"
#include "string"
#include "random"

double randomDouble(double start, double end);

class Dichotomy {
public:
    explicit Dichotomy(double start, double end, double epsilon);

    virtual ~Dichotomy() = default;

public:
    virtual double calculateRoot(double function(double x), int *div_counter);

private:
    double a_;
    double b_;
    double c_;
    double eps_;

};

class Secant : protected Dichotomy {
public:

    explicit Secant(double start,
                    double end,
                    double epsilon) : Dichotomy(start, end, epsilon) {};

    ~Secant() override = default;

public:
    double calculateRoot(double function(double x), int *div_counter) override;

private:
    double a_;
    double b_;
    double c_;
    double eps_;

};

class Gauss {
public:
    explicit Gauss(int size,
                   const std::string &matrix_filepath,
                   const std::string &vector_filepath,
                   bool create = true);

    explicit Gauss(short power,
                   std::vector<double> &matrix,
                   std::vector<double> &free_vector_);

    virtual ~Gauss() = default;

public:
    std::vector<double> calculateGauss();

private:
    std::vector<double> prototype_;
    std::vector<double> vec_;
    std::vector<double> result_vec_;

    short int size_;
};

class Newton {
public:

    explicit Newton(double function(double x),
                    const std::vector<double> &nodes,
                    double exp_point,
                    short int power);

    virtual ~Newton() = default;

public:
    double calculateNewton();

private:
    std::vector<double> nodes_;
    std::vector<double> fvalue_;;
    std::vector<double> difference_;
    std::vector<double> diag_vec_;

    short int power_;

    double exp_point_;
};

class IntegralMethod {
public:
    explicit IntegralMethod(double start_interval, double end_interval, double eps, int step);

    virtual ~IntegralMethod() = default;

public:
    virtual double calculateSimpson(double function(double x));

private:
    double start_interval_;
    double end_interval_;
    double width_;
    double eps_;

    int step_;

};

class CubeSpline {
public:
    explicit CubeSpline(short power,
                        double function(double x),
                        std::vector<double> &nodes,
                        double alpha0,
                        double alpha1,
                        double beta0,
                        double beta1,
                        double gamma0,
                        double gamma1,
                        double exp_point);

    virtual ~CubeSpline() = default;

public:
    double calculateCubeSpline();

private:
    void sweepMethod();

private:
    std::vector<double> nodes_;
    std::vector<double> fvalue_;;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
    std::vector<double> m_;
    std::vector<double> q_;
    std::vector<double> p_;

    double exp_point_;
    double alpha0_;
    double alpha1_;
    double beta0_;
    double beta1_;
    double gamma0_;
    double gamma1_;

    short power_;

};

class MonteCarlo {
public:
    explicit MonteCarlo(long power);

    virtual ~MonteCarlo() = default;

public:
    double calculateMonteCarlo1(double function(double x, double y));

    double calculateMonteCarlo2(double function(double x, double y));

private:
    double random_x_;
    double random_y_;
    double random_z_;
    double summary_;

    long power_;

};

#endif //NUMERICALMETHODSLIB_LIBRARY_H
