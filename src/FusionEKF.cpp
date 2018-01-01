#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ <<   0.0225  , 0     ,
                0       , 0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<   0.09    , 0     , 0     ,
                0       , 0.0009, 0     ,
                0       , 0     , 0.09  ;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Hj_ - dependant on state vector, 
  // initialize by filling non-zero elements with one.
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

	VectorXd init_x = VectorXd(4);
    init_x << 1, 1, 1, 1;

	MatrixXd init_P = MatrixXd(4, 4);
    init_P << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

    MatrixXd init_F = MatrixXd(4, 4);
	init_F << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

    MatrixXd init_Q = MatrixXd(4, 4);
    init_Q << 1, 0, 1, 0,
              0, 1, 0, 1,
              1, 0, 1, 0,
              0, 1, 0, 1;

    /*
    std::cout << "\n####\ninit_x= \n" << init_x << std::endl;
    std::cout << "\n####\ninit_p= \n" << init_P << std::endl;
    std::cout << "\n####\ninit_F= \n" << init_F << std::endl;
    std::cout << "\n####\ninit_Q= \n" << init_Q << std::endl;
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float rho     = measurement_pack.raw_measurements_[0];
        float phi     = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];

        init_x << rho*cos(phi),     // px
                  rho*sin(phi),     // py
                  rho_dot*cos(phi), // vx
                  rho_dot*sin(phi); // vy

        ekf_.Init(init_x, init_P, init_F, Hj_ , R_radar_, init_Q); 

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // Laser - Initialize state.
	    init_x << measurement_pack.raw_measurements_[0],    // px
                  measurement_pack.raw_measurements_[1],    // py
                  0,                                        // vx
                  0;                                        // vy

        ekf_.Init(init_x, init_P, init_F, H_laser_, R_laser_, init_Q);
    }

    // Store/initialize the first measurement timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // * Update the state transition matrix F w.r.t new elapsed time
  // compute time difference between consecutive measurement - expressed in secs
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; 

  previous_timestamp_ = measurement_pack.timestamp_; // store current timestamp for next iteration

  // update state transition matrix F
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // update process noise covariance matrix Q, Q = G*Qv*Gt
  float noise_ax = 9;
  float noise_ay = 9;

  MatrixXd G = MatrixXd(4,2);
  G <<	dt*dt/2	, 0			,
		0		, dt*dt/2	,
		dt		, 0			,
		0		, dt		;

  MatrixXd Qv = MatrixXd(2, 2);
  Qv << noise_ax	, 0			,
		0			, noise_ay	;

  ekf_.Q_ = G*Qv*G.transpose();

  // Prediction call
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);    
  }

  // print the output
  // cout << "\n####Estimated \nx_ = \n" << ekf_.x_ << endl;
  // cout << "\n####Estimated \nP_ = \n" << ekf_.P_ << endl;
}
