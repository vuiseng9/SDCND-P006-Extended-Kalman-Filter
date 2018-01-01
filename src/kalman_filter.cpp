#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
    TODO:
     * predict the state
    */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
    TODO:
     * update the state by using Kalman Filter equations
    */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    TODO:
        * update the state by using Extended Kalman Filter equations
    */
  	float px_pred = x_[0];
  	float py_pred = x_[1];
  	float vx_pred = x_[2];
  	float vy_pred = x_[3];
	
	float rho_pred = sqrt(pow(px_pred,2) + pow(py_pred,2));
	float phi_pred = atan2(py_pred,px_pred);

  	// avoid division by zero in computing rho_dot_pred
 	if(rho_pred < 0.000001)
        rho_pred = 0.000001;

  	float rho_dot_pred = (px_pred * vx_pred + py_pred * vy_pred) / rho_pred;

  	VectorXd z_pred = VectorXd(3);
  	z_pred << rho_pred, phi_pred, rho_dot_pred;

	VectorXd y = z - z_pred;

	// normalize the angle between -pi to pi
  	if (y(1) >  M_PI)
        y(1) -= 2 * M_PI;

  	if (y(1) < -M_PI)
    	y(1) += 2 * M_PI;
 
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
