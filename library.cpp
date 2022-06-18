#include "library.h"
#include "cmath"
#include <fstream>

double randomDouble(double start, double end) {
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distribution(start, end);
    return distribution(eng);
}

Dichotomy::Dichotomy(double start,
                     double end,
                     double epsilon) :
        a_(start),
        b_(end),
        eps_(epsilon) {

}

double Dichotomy::calculateRoot(double (*function)(double),
                                int *div_counter) {

    double function_a_ = function(a_);
    double function_b_ = function(b_);
    double function_c_ = 0;

    if (function_a_ * function_b_ > 0) {
        std::cout << "Function has no root on this section." << "\n";
        return 0;
    }

    while (fabs(b_ - a_) >= eps_) {
        c_ = (a_ + b_) * 0.5;
        function_c_ = function(c_);
        *div_counter += 1;

        if (function_a_ * function_c_ <= 0) {
            b_ = c_;
            function_b_ = function_c_;
        } else {
            a_ = c_;
            function_a_ = function_c_;
        }

        if (fabs(b_ - a_) < eps_) {
            return c_;
        }
    }

    return 0;
}

double Secant::calculateRoot(double (*function)(double), int *div_counter) {
    double temp;

    double function_a_ = function(a_);
    double function_b_ = function(b_);
    double function_c_ = 0;

    while (function_a_ * function_b_ < 0) {
        c_ = (-function_a_ * b_ + function_b_ * a_) / (-function_a_ + function_b_);
        function_c_ = function(c_);
        *div_counter += 1;

        if (function_a_ * function_c_ <= 0) {
            b_ = c_;
            function_b_ = function_c_;
        } else {
            a_ = c_;
            function_a_ = function_c_;
        }

        if (fabs(function_c_) < eps_ && fabs(c_ - temp) < eps_) {
            return c_;
        }

        temp = c_;
    }

    return 0;
}

Gauss::Gauss(short size,
             const std::string &matrix_filepath,
             const std::string &vector_filepath,
             bool create) : size_(size) {

    std::ifstream numbers;
    numbers.open(matrix_filepath);

    prototype_.resize(size_ * size_);
    vec_.resize(size_);

    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            numbers >> prototype_[i * size_ + j];
        }
    }
    numbers.close();

    std::ifstream free_row;
    free_row.open(vector_filepath);

    if (create) {
        for (int i = 0; i < size_; i++) {
            free_row >> vec_[i];
        }
    }

    free_row.close();
}

Gauss::Gauss(short power,
             std::vector<double> &matrix,
             std::vector<double> &free_vector) :
        size_(power),
        prototype_(matrix),
        vec_(free_vector) {
}

std::vector<double> Gauss::calculateGauss() {
    int index_;

    double maximum_;
    double temp_;


    result_vec_.resize(size_);

    result_vec_ = vec_;

    for (int rows = 0; rows < size_; rows++) {
        maximum_ = fabs(prototype_[rows * size_ + rows]);
        index_ = rows;

        for (int i = rows + 1; i < size_; i++) {
            if (fabs(prototype_[rows * size_ + rows]) > maximum_) {
                maximum_ = fabs(prototype_[rows * size_ + rows]);
                index_ = i;
            }
        }

        if (maximum_ < 1e-15) {
            std::cout << "Matrix has no solution" << "\n";
            return {};
        }

        for (int j = 0; j < size_; j++) {
            std::swap(prototype_[rows * size_ + j], prototype_[index_ * size_ + j]);
        }

        std::swap(vec_[rows], vec_[index_]);

        for (int i = rows + 1; i < size_; i++) {
            temp_ = prototype_[i * size_ + rows] / prototype_[rows * size_ + rows];
            for (int j = rows; j < size_; j++) {
                prototype_[i * size_ + j] -= temp_ * prototype_[rows * size_ + j];
            }
            vec_[i] -= temp_ * vec_[rows];
        }
    }

    temp_ = 0.0;

    if (prototype_[(size_ - 1) * size_ + (size_ - 1)] == 0) {
        if (vec_[size_ - 1] == 0) {
            std::cout << "Matrix has infinite solutions." << "\n";
            return {};
        } else {
            std::cout << "Matrix has no solution." << "\n";
            return {};
        }
    } else {
        for (int i = size_ - 1; i >= 0; i--) {
            temp_ = 0.0;
            for (int j = i + 1; j < size_; j++) {
                temp_ += prototype_[i * size_ + j] * result_vec_[j];
            }
            result_vec_[i] = (vec_[i] - temp_) / prototype_[i * size_ + i];
        };
    }

    return result_vec_;
}

Newton::Newton(double (*function)(double),
               const std::vector<double> &nodes,
               double exp_point,
               short power) :
        exp_point_(exp_point),
        power_(power + 1) {

    fvalue_.resize(power_);

    nodes_.resize(power_);
    nodes_ = nodes;

    for (int i = 0; i <= nodes_.size(); i++) {
        fvalue_[i] = function(nodes[i]);
    }

    difference_.resize(power_ * power_);

    for (int i = 0; i < power_; i++) {
        difference_[i * power_ + 0] = fvalue_[i];
    }

}

double Newton::calculateNewton() {
    double mul = 1.0;
    double result = 0.0;

    diag_vec_.resize(power_ * power_);

    for (int i = 1; i < power_; i++) {
        for (int j = 1; j < power_; j++) {
            if (i < j) {
                difference_[i * power_ + j] = 0;
            } else {
                difference_[i * power_ + j] =
                        (difference_[i * power_ + (j - 1)] -
                         difference_[(j - 1) * power_ + (j - 1)]) /
                        (nodes_[i] - nodes_[j - 1]);
            }
        }
    }

    for (int i = 0; i < power_; i++) {
        diag_vec_[i] = difference_[i * power_ + i];
    }

    for (int i = 0; i < power_; i++) {
        mul = 1.0;
        for (int j = 0; j < i; j++) {
            mul *= (exp_point_ - nodes_[j]);
        }
        result += diag_vec_[i] * mul;
        std::cout << "Interval value " << i << ": " << result << "\n";
    }

    return result;
}

IntegralMethod::IntegralMethod(double start, double end, double eps, int step) :
        start_(start),
        end_(end),
        width_((end - start) / step),
        eps_(eps),
        step_(step) {

}

double IntegralMethod::calculateSimpson(double function(double x)) {
    const double temp = 15. / 16;

    double integral1;
    double integral2;

    double s0 = function(start_) + function(end_);
    double s1 = 0.0;
    double s2 = 0.0;

    eps_ *= temp;

    for (int i = 1; i < step_; i += 2) {
        s1 += function(start_ + (i * width_));
    }

    for (int i = 0; i < step_; i += 2) {
        s2 += function(start_ + (i * width_));
    }

    integral1 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);

    step_ *= 2;
    width_ /= 2;

    s2 += s1;
    s1 = 0;

    for (int i = 1; i < step_; i += 2) {
        s1 += function(start_ + (i * width_));
    }

    integral2 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);

    while (fabs(integral2 - integral1) >= eps_) {
        step_ *= 2;
        width_ /= 2;

        integral1 = integral2;
        s2 += s1;
        s1 = 0;

        for (int i = 1; i < step_; i += 2) {
            s1 += function(start_ + (i * width_));
        }

        integral2 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);
    }

    return integral2;
}

CubeSpline::CubeSpline(short power,
                       double (*function)(double),
                       std::vector<double> &nodes,
                       double alpha0,
                       double alpha1,
                       double beta0,
                       double beta1,
                       double gamma0,
                       double gamma1,
                       double exp_point) :
        alpha0_(alpha0),
        alpha1_(alpha1),
        beta0_(beta0),
        beta1_(beta1),
        gamma0_(gamma0),
        gamma1_(gamma1),
        exp_point_(exp_point),
        nodes_(nodes),
        power_(power - 1) {

    a_.resize(power_ + 1);
    b_.resize(power_ + 1);
    c_.resize(power_ + 1);
    d_.resize(power_ + 1);

    p_.resize(power_);
    q_.resize(power_);
    m_.resize(power_ + 1);

    fvalue_.resize(power_ + 1);

    for (int i = 0; i < power_ + 1; i++) {
        fvalue_[i] = function(nodes_[i]);
    }

}

void CubeSpline::sweepMethod() {
    p_[0] = -c_[0] / b_[0];
    q_[0] = d_[0] / b_[0];

    for (int i = 1; i < power_; i++) {
        p_[i] = -c_[i] / (a_[i] * p_[i - 1] + b_[i]);
        q_[i] = (d_[i] - a_[i] * q_[i - 1]) / (a_[i] * p_[i - 1] + b_[i]);
    }

    m_[power_] = (d_[power_] - a_[power_] * q_[power_ - 1]) / (a_[power_] * p_[power_ - 1] + b_[power_]);

    for (int i = power_; i > 0; i--) {
        m_[i - 1] = p_[i - 1] * m_[i] + q_[i - 1];
    }

}

double CubeSpline::calculateCubeSpline() {
    double result_ = 0;

    int index = 0;

    a_[0] = 0;
    b_[0] = -alpha0_ * (nodes_[1] - nodes_[0]) / 3 + beta0_;
    c_[0] = alpha0_ * (nodes_[1] - nodes_[0]) / 6;
    d_[0] = gamma0_ - alpha0_ * (fvalue_[1] - fvalue_[0]) / (nodes_[1] - nodes_[0]);

    for (int i = 1; i < power_; i++) {
        a_[i] = (nodes_[i] - nodes_[i - 1]) / 6;
        b_[i] = (nodes_[i + 1] - nodes_[i - 1]) / 3;
        c_[i] = (nodes_[i + 1] - nodes_[i]) / 6;
        d_[i] = (fvalue_[i + 1] - fvalue_[i]) / (nodes_[i + 1] - nodes_[i]) -
                (fvalue_[i] - fvalue_[i - 1]) / (nodes_[i] - nodes_[i - 1]);
    }

    a_[power_] = alpha1_ * (nodes_[power_] - nodes_[power_ - 1]) / 6;
    b_[power_] = alpha1_ * (nodes_[power_] - nodes_[power_ - 1]) / 3 + beta1_;
    c_[power_] = 0;
    d_[power_] = gamma1_ - alpha1_ * (fvalue_[power_] - fvalue_[power_ - 1])
                           / (nodes_[power_] - nodes_[power_ - 1]);

    sweepMethod();

    for (int i = 1; i <= power_; i++) {
        if (exp_point_ >= nodes_[i - 1] && exp_point_ <= nodes_[i]) {
            index = i;
            break;
        }
    }

    result_ = m_[index - 1] * pow(nodes_[index] - exp_point_, 3) / (6 * (nodes_[index] - nodes_[index - 1])) +
              m_[index] * pow(exp_point_ - nodes_[index - 1], 3) / (6 * (nodes_[index] - nodes_[index - 1])) +
              (fvalue_[index - 1] - (m_[index - 1] * pow(nodes_[index] - nodes_[index - 1], 2)) / 6) *
              (nodes_[index] - exp_point_) / (nodes_[index] - nodes_[index - 1]) +
              (fvalue_[index] - (m_[index] * pow(nodes_[index] - nodes_[index - 1], 2)) / 6) *
              (exp_point_ - nodes_[index - 1]) / (nodes_[index] - nodes_[index - 1]);

    return result_;
}

MonteCarlo::MonteCarlo(long power) :
        power_(power),
        random_x_(0),
        random_y_(0),
        random_z_(0) {

}

double MonteCarlo::calculateMonteCarlo1(double function(double x, double y)) {
    for (int i = 0; i <= power_; i++) {
        random_x_ = randomDouble(0, 1);
        random_y_ = randomDouble(0, 1);

        summary_ += function(random_x_, random_y_);
    }

    return summary_ / power_;
}

double MonteCarlo::calculateMonteCarlo2(double (*function)(double, double)) {
    double counter = 0.0;

    for (int i = 0; i <= power_; i++) {
        random_x_ = randomDouble(0, 1);
        random_y_ = randomDouble(0, 1);
        random_z_ = randomDouble(0, 1);

        if (function(random_x_, random_y_) >= random_z_) {
            counter += 1.0;
        }
    }

    return double(counter) / double(power_);
}

