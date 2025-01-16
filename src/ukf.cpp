#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);
  
  n_x_ = 5;
  n_aug_ = 7;

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  

  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  lambda_ = 3 - n_aug_;
  
  MatrixXd Xsig_out = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  

  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(n_x_) = x_;
  
  for (int i = n_x_; i < n_aug_; i++){
  	x_aug(i) = 0;
  }

  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  *Xsig_out = Xsig_aug;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

VectorXd UKF::motion_model(VectorXd SigmaVal, double delta_t){
  VectorXd x_dot =  x_;
  
  double x = SigmaVal(0);
  double y = SigmaVal(1);
  double v = SigmaVal(2);
  double omega = SigmaVal(3);
  double omega_dot = SigmaVal(4);
  double tau_a = SigmaVal(5);
  double tau_w = SigmaVal(6);
  
  x_dot(0) += v / omega * (sin(omega + omega_dot * delta_t) - sin(omega)) + 1 / 2 * delta_t * delta_t * cos(omega) * tau_a;
  x_dot(1) += v / omega * (-cos(omega + omega_dot * delta_t) + cos(omega)) + 1 / 2 * delta_t * delta_t * sin(omega) * tau_a;
  x_dot(2) += delta_t * tau_a;
  x_dot(3) += omega_dot * delta_t + 1 / 2 * delta_t * delta_t * tau_w;
  x_dot(4) += delta_t * tau_w;
  
  return x_dot;
                           
  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  AugmentedSigmaPoints(&Xsig_out);
  
  MatrixXd transformed_points = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd SigmaVal = Xsig_out.col(i);
    motion_model(SigmaVal, delta_t);
  }
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}