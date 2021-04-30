#include "CBLAB_Taegun.h"


int main() {

	std::ifstream input_x("X.txt");
	std::ifstream input_y("Y.txt");

	/*std::ifstream input_UX("rhs.txt");
	std::ifstream input_HS("HS.txt");
	*/
	//std::ifstream input_qrhs("qrhs.txt");
	//std::ifstream input_p("test_p.txt");

	std::ofstream kinship("kinship_cpp.txt");
	
	std::vector<double> vc_1;
	std::vector<double> vc_2;

	double median_vc1 = 0.0;
	double median_vc2 = 0.0;

	std::cout << "Read X, Y" << std::endl;
	Eigen::MatrixXd X = read_mat(input_x, count_matrix_row(input_x), count_matrix_col(input_x));
	Eigen::MatrixXd Y = read_mat(input_y, count_matrix_row(input_y), count_matrix_col(input_y));
	//Eigen::MatrixXd qrhs = read_mat(input_qrhs, count_matrix_row(input_qrhs), count_matrix_col(input_qrhs));
	//Eigen::MatrixXd p = read_mat(input_p, count_matrix_row(input_p), count_matrix_col(input_p));


	std::cout << "EstimateKinship..." << std::endl;

	Eigen::MatrixXd K = estimateKinship(X);
	kinship << K << std::endl;
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

	Eigen::MatrixXd UX = Eigen::MatrixXd(Ky.rows(), 2);
	Eigen::MatrixXd UX_Q = Eigen::MatrixXd(Ky.rows(), 2);

	for (int z = 0; z < Kx.cols(); z++) {
		std::cout << "z = " << z << std::endl;

		for (int i = 0; i < Ky.rows(); i++) {
			UX(i, 0) = 1; UX(i, 1) = Kx(i, z);
		}

		Eigen::HouseholderQR<Eigen::MatrixXd> qr(UX);

		Eigen::MatrixXd HS_Q = qr.householderQ();

		for (int i = 0; i < Ky.rows(); i++) {
			UX_Q(i, 0) = HS_Q(i, 0); UX_Q(i, 1) = HS_Q(i, 1);
		}

		double min_Ky = Ky.minCoeff();

		Eigen::MatrixXd HS = UX_Q * UX_Q.transpose().eval();


		for (int i = 0; i < Ky.rows(); i++) {
			for (int j = 0; j < Ky.cols(); j++) {
				Ky(i, j) = Ky(i, j) - min_Ky;
			}
		}


		//calculate dissimilarity ^2 dmat
		Eigen::MatrixXd dis_Ky(Ky.rows(), Ky.rows());
		Eigen::MatrixXd dmat(Ky.rows(), Ky.rows());

		for (int i = 0; i < Ky.rows(); i++) {
			for (int k = 0; k < Ky.rows(); k++) {
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (int p = 0; p < Ky.cols(); p++) {
					sum1 = sum1 + std::abs(Ky(i, p) - Ky(k, p));
					sum2 = sum2 + Ky(i, p) + Ky(k, p);
				}
				dis_Ky(i, k) = sum1 / sum2;
				dmat(i, k) = dis_Ky(i, k) * dis_Ky(i, k);
			}
		}


		Eigen::MatrixXd G = -(dmat.colwise() - dmat.rowwise().mean()) / 2;

		double ss_exp_comb = 0.0;

		for (int i = 0; i < Ky.rows(); i++) {
			for (int k = 0; k < Ky.rows(); k++) {
				ss_exp_comb += G(i, k) * HS(i, k);
			}
		}
		// in Gamma.R ss_exp_comb != ss_exp_each, but in our method, ss_exp_comb = ss_exp_each  
	//	std::cout << ss_exp_comb << std::endl;

		Eigen::MatrixXd imat;
		Eigen::MatrixXd TIH_snterm = imat.setIdentity(Ky.rows(), Ky.cols()) - HS;


		double SS_Res = 0.0;
		for (int i = 0; i < Ky.rows(); i++) {
			for (int k = 0; k < Ky.rows(); k++) {
				SS_Res += G(i, k) * TIH_snterm(i, k);
			}
		}

		//std::cout << SS_Res << std::endl;

		//in sample df.Exp = 1; df.Res = 94
		double df_Exp = 1;
		int df_Res = Ky.rows() - 2;
		double F_Mod = (ss_exp_comb / df_Exp) / (SS_Res / df_Res);

		double S1_xx = 0.0; double S1_xy = 0.0;
		double S2_xx = 0.0; double S2_xy = 0.0;

		double S1_x_mean = UX.col(1).mean();
		double S1_y_mean = dis_Ky.col(0).mean();

		for (int i = 0; i < Ky.rows(); i++) {
			S1_xx += (UX(i, 1) - S1_x_mean) * (UX(i, 1) - S1_x_mean);
			S1_xy += (UX(i, 1) - S1_x_mean) * (dis_Ky(i, 0) - S1_y_mean);
		}


		std::cout << "ky_rows " << Ky.rows()<< std::endl;
		std::cout << "ky_cols " << Ky.cols() << std::endl;

		double b1 = (S1_xx * S1_xy) / (S1_xx * S1_xx);
		double b0 = S1_y_mean - b1 * S1_x_mean;

	//	std::cout << "b1 " << b1 << std::endl;
	//	std::cout << "b0 " << b0 << std::endl;

		double p_val = 0.0;
		for (int perm = 2; perm < 5; perm++) {
			int permute_num = std::pow(10, perm);
			std::cout << "permute_num is " << permute_num << std::endl;
			Eigen::MatrixXi permutation = Eigen::MatrixXi(permute_num, Ky.rows());
			auto cur_time = std::chrono::system_clock::now();
			auto duration = cur_time.time_since_epoch();
			auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
			std::vector<int> v(Kx.rows());
			std::iota(v.begin(), v.end(), 1);
			for (int i = 0; i < permute_num; i++) {
				ran.seed(i + millis);
				std::shuffle(v.begin(), v.end(), ran);
				for (int j = 0; j < Ky.rows(); j++) {
					permutation(i, j) = v[j];
				}
			}
			// permutation end;

			//fucntiono f_test
			std::vector<double> f_perms;
			for (int i = 0; i < permute_num; i++) {
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (int j = 0; j < Ky.rows(); j++) {
					for (int k = 0; k < Ky.rows(); k++) {
						sum1 += G((int)permutation(i, j) - 1, (int)permutation(i, k) - 1) * HS(j, k);
						sum2 += G((int)permutation(i, j) - 1, (int)permutation(i, k) - 1) * TIH_snterm(j, k);
					}
				}
				f_perms.push_back(sum1 / (sum2 / df_Res));
				//			std::cout << f_perms[i] << std::endl;
			}

			int count = 0;
			for (int i = 0; i < f_perms.size(); i++) {
				if (f_perms[i] >= F_Mod - 1.490116e-08) count++;
			}
			count++;
			p_val = (double)count / (permute_num + 1);
			if (p_val > 5.0 / permute_num) {
				std::cout << "p_val is " << p_val << " vs " << 5.0 / permute_num <<std::endl;
				break;
			}
		}
		std::cout << F_Mod << " " << p_val << std::endl;

	}
	input_x.close();
	input_y.close();
	kinship.close();
	//input_qrhs.close();
	//input_p.close();


	/*input_UX.close();
	input_HS.close();*/

	return 0;
}