#include "CBLAB_method.h"


int main(int argc, char** argv) {

	cxxopts::Options options("GAMMA-RAY", "A GAMMA analysis tool written in C++");

	options.add_options()
		("x,input_x", "Path to input Genotpye file", cxxopts::value<std::string>())
		("y,input_y", "Path to input Phenotype file", cxxopts::value<std::string>())
		("o,output", "Path to output file", cxxopts::value<std::string>())
		("t,threads", "Number of threads (default: 1)", cxxopts::value<int>()->default_value("1"))
		("p,permutations", "Number of permutations (default: 4)", cxxopts::value<int>()->default_value("4"))
		("h,help", "Print help");

	auto result = options.parse(argc, argv);

	if (result.count("help") || !result.count("input_x") || !result.count("input_y") || !result.count("output")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	std::string input_x_path = result["input_x"].as<std::string>();
	std::string input_y_path = result["input_y"].as<std::string>();
	std::string output_path = result["output"].as<std::string>();
	int thread_num = result["threads"].as<int>();
	int permutation_num = result["permutations"].as<int>();

	std::ifstream input_x(input_x_path);
	std::ifstream input_y(input_y_path);
	std::ofstream output(output_path);

	if (!input_x || !input_y || !output) {
		std::cerr << "Error opening one or more files." << std::endl;
		return 1;
	}

	std::cout << "Input X: " << input_x_path << "\n";
	std::cout << "Input Y: " << input_y_path << "\n";
	std::cout << "Output: " << output_path << "\n";
	std::cout << "Threads: " << thread_num << "\n";
	std::cout << "Permutations: " << permutation_num << "\n";

	// º» ÇÁ·Î±×·¥ ·ÎÁ÷...
	std::vector<double> vc_1;
	std::vector<double> vc_2;

	double median_vc1 = 0.0;
	double median_vc2 = 0.0;

	std::cout << "Read X, Y" << std::endl;

	Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
	Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));

	std::cout << "EstimateKinship..." << std::endl;
//	Eigen::MatrixXd K = read_mat(kinship, count_matrix_row(kinship), count_matrix_col(kinship)); if data sample size too big
	Eigen::MatrixXd K = calculate_kinship(X); 

	std::cout << "EstimateVarComp..." << std::endl;

	estimateVarComp(K, X, Y, vc_1, vc_2); //calculate vc_1, vc_2 /////////////////////////

	std::cout << "Calculate Median" << std::endl;
	median_vc1 = cal_median(vc_1); median_vc2 = cal_median(vc_2);

	std::cout << "Calculate SigmaMatrix" << std::endl;
	
	std::cout << K.rows() << " " << K.cols() << std::endl;

	std::cout << median_vc1 << " " << median_vc2 << std::endl;

	//median_vc1 = 0;
	cal_SigmaM(K, median_vc1, median_vc2); // after this function, K is sigmaMatrix

	std::cout << K.rows() << " " << K.cols() << std::endl;
	
	X = X.transpose().eval();
	Y = Y.transpose().eval();

	Eigen::MatrixXd Kx = K.transpose().eval() * X;
	Eigen::MatrixXd Ky = K.transpose().eval() * Y;

	//Gamma_cpp(Kx, Ky);
	
	std::future<std::vector<std::string>>* thread_Gamma = new std::future<std::vector<std::string>>[thread_num];

	for (int i = 0; i < thread_num; i++) {
		thread_Gamma[i] = std::async(std::launch::async, Gamma_cpp, std::ref(Kx), std::ref(Ky), i, thread_num, permutation_num);
	}
	for (int i = 0; i < thread_num; i++) {
		thread_Gamma[i].wait();
	}
	for (int i = 0; i < thread_num; i++) {
		std::vector<std::string> temp_arr = thread_Gamma[i].get();
		int tmp_size = temp_arr.size();
		for (int k = 0; k < tmp_size; k++) {
			output << temp_arr[k] << std::endl;
		}
	}

	input_x.close();
	input_y.close();
	output.close();
	return 0;

}
