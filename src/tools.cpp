#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   	/**
  	TODO:
      * Calculate the RMSE here.
  	*/

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  	/**
  	TODO:
     * Calculate a Jacobian here.
  	*/

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE
	float sumofsq = px*px + py*py;
    float denom1 = sqrt(sumofsq);
    float denom2 = pow(denom1,3);
    
	//check division by zero
	if (fabs(sumofsq) < 0.0001 )
	{
	    cout << "Hey! divide by zero" << endl;
	    return Hj;
	}
	else
    {
        //compute the Jacobian matrix
        Hj <<   px/denom1, py/denom1, 0, 0,
                -1*py/sumofsq, 	px/sumofsq, 0,0,
                py*(vx*py-vy*px)/denom2, px*(vy*px-vx*py)/denom2, px/denom1, py/denom1;
    }

	return Hj;
}
