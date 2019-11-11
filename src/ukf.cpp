#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  time_us_ = 0.0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   cov_radar = MatrixXd(3,3);

   cov_radar << std_radr_*std_radr_,0,0,
                0,std_radphi_*std_radphi_,0,
                0,0,std_radrd_*std_radrd_;

   cov_lidar = MatrixXd(2,2);

  cov_lidar << std_laspx_*std_laspx_,0,
               0,std_laspy_*std_laspy_;
   n_x_ = 5;
   n_aug_ = n_x_ + 2;

   weights_ = VectorXd(2 * n_aug_ + 1);
   weights_.fill(0.5/(lambda_ + n_aug_));
   weights_(0) = lambda_/(lambda_ + n_aug_);

  P_ = MatrixXd(n_x_, n_x_);
  P_ << 0.1, 0, 0, 0, 0,
        0, 0.2, 0, 0, 0,
        0, 0, 0.4, 0, 0,
        0, 0, 0, 0.2, 0,
        0, 0, 0, 0, 0.1;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
}

UKF::~UKF() {}



void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //std::cout<<"started///////"<<std::endl;
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   if(!is_initialized_) {
     if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0]; // range
      double phi = meas_package.raw_measurements_[1]; // bearing
      double rho_dot = meas_package.raw_measurements_[2]; //velocity in radial direction
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      x_ << x , y, v, 0, 0;
     }
     else {
       //take sensor measurements from lidar
       x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
     }

     is_initialized_ = true;
     return;
   }

   // calculate delta t
   double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
   time_us_ = meas_package.timestamp_;

   // Prediction step
   Prediction(delta_t);


   // Update step
   if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
     UpdateRadar(meas_package);
   }

   if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
     UpdateLidar(meas_package);
   }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
   double lambda_ = 3 - n_aug_;
   VectorXd x_aug = VectorXd(n_aug_);
   MatrixXd p_aug = MatrixXd(n_aug_, n_aug_);
   MatrixXd x_aug_temp = MatrixXd(n_aug_, 2*n_aug_+1);

   x_aug.head(5) = x_;
   x_aug(5) = 0.0;
   x_aug(6) = 0.0;

   p_aug.fill(0);
   p_aug.topLeftCorner(5,5) = P_;
   p_aug(5,5) = std_a_*std_a_;
   p_aug(6,6) = std_yawdd_*std_yawdd_;
   //calculate square root of p_aug
   MatrixXd L = p_aug.llt().matrixL();
   x_aug_temp.col(0) = x_aug;


   for(size_t i{0}; i<n_aug_; ++i) {
     x_aug_temp.col(i+1)  = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
     x_aug_temp.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
   }


   //predict sigma points generated above

   for(size_t i{0}; i<2*n_aug_+1; ++i) {
     double p_x = x_aug_temp(0,i);
     double p_y = x_aug_temp(1,i);
     double v = x_aug_temp(2,i);
     double yaw = x_aug_temp(3,i);
     double yawd = x_aug_temp(4,i);
     double nu_a = x_aug_temp(5,i);
     double nu_yawd = x_aug_temp(6,i);

     double px_p, py_p, v_p, yaw_p, yawd_p;

     // catch division by zero error
     if(fabs(yawd) > 0.001) {
       px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
       py_p = p_y + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
     } else {
       px_p = p_x + v * delta_t * cos(yaw);
       py_p = p_y + v * delta_t * sin(yaw);
     }
     v_p = v;
     yaw_p = yaw + yawd * delta_t;
     yawd_p = yawd;

     // add noise
     px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
     py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
     v_p = v_p + nu_a * delta_t;
     yaw_p = yaw_p + 0.5 * nu_yawd * delta_t * delta_t;
     yawd_p = yawd_p + nu_yawd * delta_t;

     Xsig_pred_(0, i) = px_p;
     Xsig_pred_(1, i) = py_p;
     Xsig_pred_(2, i) = v_p;
     Xsig_pred_(3, i) = yaw_p;
     Xsig_pred_(4, i) = yawd_p;
   }


   //estimate mean of our state vector
   x_.fill(0);
   for (size_t i{0}; i<2*n_aug_+1; i++) {
     x_ = x_ + weights_(i)*Xsig_pred_.col(i);
   }


   //estimate covariance of our state vector
   P_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     // normalize angle
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
     P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
   }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
   VectorXd z_ = meas_package.raw_measurements_;
   uint dim_z = 2;
   MatrixXd z_sig = MatrixXd(dim_z, 2*n_aug_+1);
   for(size_t i{0}; i< 2*n_aug_+1; ++i) {
     z_sig(0,i) = Xsig_pred_(0,i);
     z_sig(1,i) = Xsig_pred_(1,i);
   }


   //
   VectorXd z_pred_mean = VectorXd(dim_z);
   z_pred_mean.fill(0);
   for (size_t i{0}; i < 2*n_aug_+1; ++i) {
     z_pred_mean = z_pred_mean + weights_(i)*z_sig.col(i);
   }

   MatrixXd S = MatrixXd(dim_z, dim_z);
   S.fill(0);
   for(size_t i{0}; i<2*n_aug_+1; i++) {
     VectorXd z_diff = z_sig.col(i) - z_pred_mean;
     S = S + weights_(i)*z_diff*z_diff.transpose();
   }
   S = S + cov_lidar;

   MatrixXd Tc = MatrixXd(n_x_, dim_z);
   Tc.fill(0.0);

   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

     VectorXd z_diff = z_sig.col(i) - z_pred_mean;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }
   MatrixXd K = Tc * S.inverse();
   VectorXd z_diff = z_ - z_pred_mean;
   x_ = x_ + K*z_diff;
   P_ = P_ - K*S*K.transpose();


}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   VectorXd z_ = meas_package.raw_measurements_;

   uint dim_z = 3;
   MatrixXd z_sig = MatrixXd(dim_z, 2*n_aug_+1);



   for(int i = 0; i < 2 * n_aug_ + 1; i++){
     double p_x = Xsig_pred_(0, i);
     double p_y = Xsig_pred_(1, i);
     double v = Xsig_pred_(2, i);
     double yaw = Xsig_pred_(3, i);
     double yawd = Xsig_pred_(4, i);

     double vx = cos(yaw)*v;
     double vy = sin(yaw)*v;

     z_sig(0, i) = sqrt(p_x*p_x + p_y*p_y);                      // r
     z_sig(1, i) = atan2(p_y, p_x);                              // phi
     z_sig(2, i) = (p_x*vx + p_y*vy)/(sqrt(p_x*p_x + p_y*p_y));  // r_dot
   }


   // calculate mean predicted measurement
   VectorXd z_pred_mean = VectorXd(dim_z);
   z_pred_mean.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     z_pred_mean = z_pred_mean + weights_(i)*z_sig.col(i);
   }
   // calculate covariance of predicted measurement
   MatrixXd S = MatrixXd(dim_z, dim_z);
   S.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd z_diff = z_sig.col(i) - z_pred_mean;

     while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
     while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   // add measurement noise covariance matrix
   S = S + cov_radar;

   // UKF update
   // Cross correlation matrix between sigma points in state space
   // and measurement space
   MatrixXd Tc = MatrixXd(n_x_, dim_z);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

     VectorXd z_diff = z_sig.col(i) - z_pred_mean;
     while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
     while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }


   // calculate Kalman gain K
   MatrixXd K = Tc * S.inverse();

   // update state mean and covariance
   // residual
   VectorXd z_diff = z_ - z_pred_mean;
   while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
   while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

   x_ = x_ + K*z_diff;

   P_ = P_ - K*S*K.transpose();

}
