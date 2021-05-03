#include "CBLAB_Taegun.h"


int main(int argc, char** argv) {
	
	std::ifstream input_x("X.txt");
	std::ifstream input_y("Y.txt");
	int thread_num = 20;
	std::vector<double> vc_1;
	std::vector<double> vc_2;

	double median_vc1 = 0.0;
	double median_vc2 = 0.0;

	std::cout << "Read X, Y" << std::endl;

	Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
	Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));

	std::cout << "EstimateKinship..." << std::endl;

	Eigen::MatrixXd K = estimateKinship(X);

	std::cout << "EstimateVarComp..." << std::endl;

	estimateVarComp(K, X, Y, vc_1, vc_2); //calculate vc_1, vc_2

	std::cout << "Calculate Median" << std::endl;
	median_vc1 = cal_median(vc_1); median_vc2 = cal_median(vc_2);

	std::cout << "Calculate SigmaMatrix" << std::endl;

	cal_SigmaM(K, median_vc1, median_vc2); // after this function, K is sigmaMatrix
	
	X = X.transpose().eval();
	Y = Y.transpose().eval();

	Eigen::MatrixXd Kx = K.transpose().eval() * X;
	Eigen::MatrixXd Ky = K.transpose().eval() * Y;

	//Gamma_cpp(Kx, Ky);
	
	double min_Ky = Ky.minCoeff();
	
	for (int i = 0; i < Ky.rows(); i++) {
	    for (int j = 0; j < Ky.cols(); j++) {
	        Ky(i, j) = Ky(i, j) - min_Ky;
	    }
	}

	std::future<std::vector<std::string>>* thread_Gamma = new std::future<std::vector<std::string>>[thread_num];

	for (int i = 0; i < thread_num; i++) {
		thread_Gamma[i] = std::async(std::launch::async, Gamma_cpp, std::ref(Kx), std::ref(Ky), i, thread_num);
	}
	for (int i = 0; i < thread_num; i++) {
		thread_Gamma[i].wait();
	}
	for (int i = 0; i < thread_num; i++) {
		std::vector<std::string> temp_arr = thread_Gamma[i].get();
		int tmp_size = temp_arr.size();
		for (int k = 0; k < tmp_size; k++) {
			std::cout << temp_arr[k] << std::endl;
		}
	}

	input_x.close();
	input_y.close();
	//input_qrhs.close();
	//input_p.close();

	/*input_UX.close();
	input_HS.close();*/

	return 0;
}