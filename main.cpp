#include <iostream>
#include <complex>
#include <cmath>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include<opencv2/opencv.hpp>

typedef double precision_t;
typedef std::complex<precision_t> complex_t;

void create_complex_plane(const Eigen::Matrix<complex_t, Eigen::Dynamic, 1>& x_axis,
                          const Eigen::Matrix<complex_t, Eigen::Dynamic, 1>& y_axis,
                          Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>& complex_plane) {
    
    for(int r = 0; r < x_axis.rows(); r++) {
        for(int c = 0; c < y_axis.rows(); c++) {
            complex_plane(c, r) = x_axis(r, 0) + y_axis(c, 0);
        }
    }
}

void thread_mandelbrot_map(const Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>& c0,
                          Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>& z0,
                          int maxiter, int row_start, int row_end, int col_start, int col_end) {

    for(int r = row_start; r < row_end; r++) {
        for(int c = col_start; c < col_end; c++) {
            for(int iters = 0; iters < maxiter; iters++) {
                z0(r, c) = std::pow(z0(r, c), 2) + c0(r, c);

                if(std::pow(z0(r, c).real(), 2) + std::pow(z0(r, c).imag(), 2) > 4) {
                    z0(r, c) = 4 + iters;
                    break;
                }
            }

        }
    }

    std::cout << "Finished: " << "(" << row_start << ", " << row_end << ") x (" << col_start << ", " << col_end << ")" <<"\n";
}

void apply_mandelbrot_map(const Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>& c0, Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>& z0, int maxiter) {
    std::vector<std::thread> thread_vector;

    int max_rows = z0.rows();
    int max_cols = z0.cols();
    int step_size = 1000;
    
    for(int row_start = 0; row_start < max_rows; row_start += step_size){
        for(int col_start = 0; col_start < max_cols; col_start += step_size) {
            int row_end = row_start + step_size < max_rows ? row_start + step_size : max_rows;
            int col_end = col_start + step_size < max_cols ? col_start + step_size : max_cols;

            // std::cout << "(" << row_start << ", " << row_end << ") x (" << col_start << ", " << col_end << ")" <<"\n";
            std::thread t(thread_mandelbrot_map, std::ref(c0), std::ref(z0), maxiter, row_start, row_end, col_start, col_end);
            thread_vector.push_back(move(t));            
        }
    }

    std::cout << "Joining" << std::endl;
    for(auto& t: thread_vector) {
        t.join();
    }
    std::cout << "Done" << std::endl;
}

int main(int argc, char* argv[]) {
    const int RESOLUTION = 10000;
    
    using namespace std::complex_literals;
    
    Eigen::Matrix<complex_t, Eigen::Dynamic, 1> x_axis;
    Eigen::Matrix<complex_t, Eigen::Dynamic, 1> y_axis;
    Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic> complex_plane;
    Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic> z0;
    
    x_axis = Eigen::RowVectorXd::LinSpaced(RESOLUTION, -1.5, 0.5);
    y_axis = 1i * Eigen::RowVectorXd::LinSpaced(RESOLUTION, -1.25, 1.25);

    complex_plane.resize(RESOLUTION, RESOLUTION);
    z0.resize(RESOLUTION, RESOLUTION);
    
    create_complex_plane(x_axis, y_axis, complex_plane);

    apply_mandelbrot_map(complex_plane, z0, 100);

    Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic> mandelbrot_plot;    
    mandelbrot_plot = z0.array().abs().matrix();

    mandelbrot_plot = (mandelbrot_plot.array().abs() < 2).select(0, mandelbrot_plot);
    mandelbrot_plot = (mandelbrot_plot.array().abs() >= 2).select(255, mandelbrot_plot);
    
    cv::Mat_<precision_t> mandelbrot_cv = cv::Mat_<precision_t>::zeros(RESOLUTION, RESOLUTION);
    cv::eigen2cv(mandelbrot_plot, mandelbrot_cv);
    cv::imwrite("mandelbrot.png", mandelbrot_cv);
    
    return 0;
}
