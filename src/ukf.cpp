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

  // state vector dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.75;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  lambda_ = 3 - n_x_;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  // Initialize weights.
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Initialize measurement noice covarieance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_,0,
              0,std_laspy_ * std_laspy_;

}


MatrixXd UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda, int n_sig) {
  int n = x.size();
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd( n, n_sig );

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  Xsig.col(0) = x;

  double lambda_plue_n_x_sqrt = sqrt(lambda + n);
  for (int i = 0; i < n; i++){
      Xsig.col( i + 1 ) = x + lambda_plue_n_x_sqrt * A.col(i);
      Xsig.col( i + 1 + n ) = x - lambda_plue_n_x_sqrt * A.col(i);
  }
  return Xsig;
}



UKF::~UKF() {}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig, double delta_t, int n_x, int n_sig, double nu_am, double nu_yawdd) {
  MatrixXd Xsig_pred = MatrixXd(n_x, n_sig);

  for (int i = 0; i< n_sig; i++)
  {

    double p_x = Xsig(0,i);
    double p_y = Xsig(1,i);
    double v = Xsig(2,i);
    double yaw = Xsig(3,i);
    double yawd = Xsig(4,i);
    double nu_a = Xsig(5,i);
    double nu_yawdd = Xsig(6,i);

    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a* delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  return Xsig_pred;
}



void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if ( !is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      double rho = meas_package.raw_measurements_[0]; 
      double phi = meas_package.raw_measurements_[1]; 
      double rho_dot = meas_package.raw_measurements_[2]; 

      double x = rho * cos(phi);
      double y = rho * sin(phi);

      double vx = rho_dot * cos(phi);
  	  double vy = rho_dot * sin(phi);

      double v = sqrt(vx * vx + vy * vy);

      x_ << x, y, v, 0, 0;

      double time_us_ = meas_package.timestamp_ ;


    } else {;
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_ ;
    is_initialized_ = true;
    return;

  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

}

void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // 1. Generate sigma points.
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;


  MatrixXd X_sig_aug = GenerateSigmaPoints(x_aug, P_aug, lambda_, 2 * n_aug_ + 1);
  Xsig_pred_ = PredictSigmaPoints(X_sig_aug, delta_t, n_x_, 2 * n_aug_ + 1, std_a_, std_yawdd_);

  x_ = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngleOnComponent(x_diff, 3);
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

void UKF::NormalizeAngleOnComponent(VectorXd vector, int index) {
  while (vector(index)> M_PI) vector(index)-=2.*M_PI;
  while (vector(index)<-M_PI) vector(index)+=2.*M_PI;
}



void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_laser = 2;
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_laser, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_laser);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_laser,n_laser);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 

    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_lidar_;


  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_laser);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;


  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();




  }

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    int n_radar = 3;
    MatrixXd Zsig = MatrixXd(n_radar, 2 * n_aug_ + 1);

    for (int i = 0; i < 2 * n_aug_ + 1; i++){

      double x = Xsig_pred_(0,i);
      double y = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double omega = Xsig_pred_(3,i);

      Zsig(0, i) = sqrt(x * x + y * y);
      Zsig(1, i) = atan2(y, x);
      Zsig(2, i) = (x * cos(omega) * v + y * sin(omega) * v) / sqrt(x * x + y * y);

    }

  VectorXd z_pred = VectorXd(n_radar);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_radar,n_radar);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar_;

  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_radar);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;


  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}
