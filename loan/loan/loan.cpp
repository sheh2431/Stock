#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector> 
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <ctime>
#include <direct.h>
#include <time.h>
using namespace std;

vector<int> stock;
vector<double> dprice[50];
vector<int> portfolio[100];
vector<int> share[100];
vector<double> s_remainder[100];
vector<double> fs[100];
vector<double> seperate_fs[100];
vector<int> sol[100];

vector<double> first_day_price; //�O���Ĩ�Ĥ@�ѽ�X����

vector<double> beta;
int s_num;
int max_s = 51; // �̦h��51��bit (0050+��Ȧ�)

int N = 10;				//�ɤl��
int T = 50;			//���禸��
int generation = 10000;	//�@�N��
double r_angle = 0.0004;//���ਤ��
double r_upper = 0.0004;//���ਤ�פW��
double r_lower = 0.0004;//���ਤ�פU��
double p = 0.5;			//��l���v



double fee = 0.001425;	//�������O	0.1425%
double stt = 0.003;		//�ҥ�|		0.3%
double rfr = 0.0087;	//�L���I�Q�v	0.87%
double init_fund = 10000000;
double invest_fund;
int day;
int train_day;
int funds[100];
double funds_remainder[100];
double exp_r[100];
double risk[100];
double fit[100];
double reward[100];
//�f�tder�L����
double sharpe_ratio;
double sha_risk;
double sha_rew;
// ii 
double constant[100];

double G_best = 0;
double G_worst = 0;
double P_best, P_worst;
int best_p;
int worst_p;
vector<int> Gbest_p;
vector<int> Gbest_share;
vector<double> Gbest_s_remainder;
vector<double> Gbest_seperate_fs[100];
vector<double> Gbest_fs;
vector<int> Gbest_portfolio;
double Gbest_reward;
double Gbest_risk;
double Gbest_exp_r;
double Gbest_fund;
int Gbest_gen = 0;

vector<int> Gworst_p;
vector<double> Gworst_fs;
vector<int> Gworst_portfolio;
int Gworst_gen = 0;


double all_best = -9999;
double all_worst = 0;
vector<int> all_Gbest_p;
vector<int> all_Gbest_share;
vector<double> all_Gbest_s_remainder;
vector<double> all_Gbest_seperate_fs[100];
vector<double> all_Gbest_fs;
vector<int> all_Gbest_portfolio;
double all_Gbest_reward;
double all_Gbest_risk;
double all_Gbest_exp_r;
double all_Gbest_fund;
int all_Gbest_times = 0;
int all_Gbest_gen = 0;

int find_times;								//T�����示���ۦPall_best������
double avg_Gbest_gen;				//�p�⥭���̨θѪ��N��
double all_the_Gbest[100];			//�O���C�����窺Gbest
double avg_Gbest;						//�C������Gbest������
double std_Gbest;						//�C������Gbest���зǮt

vector<int> all_Gworst_p;
vector<double> all_Gworst_fs;
vector<int> all_Gworst_portfolio;
int all_Gworst_gen = 0;
//int flag = 0;

vector<int> guided_p;


int period = 1;
vector<double> total_test_fs;
double total_test_exp_r;
double total_test_risk;
double total_test_reward;
double total_test_fit;
int total_test_day;
bool train_period; /* --- �V�m��: true; ���մ�: false --- */
bool train_afford = true;			//�O���V�m�����Ѳ��R���R���_�F���R���_�Otrue; ���@�ӶR���_�Ofalse
bool test_afford = true;			//�O�����մ����Ѳ��R���R���_�F���R���_�Otrue; ���@�ӶR���_�Ofalse
string window[13] = { "Y2Y", "Y2H", "Y2Q", "Y2M", "H2H", "H2Q", "H2M", "H#", "Q2Q", "Q2M", "Q#", "M2M", "M#" };

string sliding_window;

/*
string testing_period  = "�ѻ�/" + sliding_window + '/';

string target = "�V�m+����_��Ȧ�/Gbest_10000_10_50_0.0004(���L����)/" + sliding_window + '/';
string train_output = target + "�V�m��/";
string train_simple_output = target + "train_Gbest_10000_10_50_0.0004_" + sliding_window + "_all.csv";
string test_output = target + "���մ�/";
string test_simple_output = target + "test_Gbest_10000_10_50_0.0004_" + sliding_window + "_all.csv";
string whole_period_result = target + "total_test_result_" + sliding_window + ".csv";
string testing = "2010/01-2016/12";
*/
string stock_price = "2010-2017";

string convert_str(double n) {
	stringstream ss;
	ss << n;
	return ss.str();
}
//string target = "(srand114)NQTS_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + "_������";
//string target = "ENQTS_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_upper) + '-' + convert_str(r_lower);
//string target = "(srand114)�F�o_NQTS_1000_10_1_0.0004";
string target = "(�[�w�s)(srand114)"+stock_price+"_GNQTS_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_upper) + '-' + convert_str(r_lower);

//string target = "(������)0050_2010-2018";
string dir;
string testing = "";
string test_input;
bool t_str = false;
bool first;
int mode = 1; // mode: 0, 0050������; 1, 0050������+�w�s(9999)

bool output_test = false;
bool run_out = false;
bool maximum = true;
bool fee_and_tax = false; //�O�_�Ҽ{����O�ε| (true: ���Ffalse: �_)
bool odd_lots = false;

string method;

void make_file() {
	dir = ".\\�Ĩ�\\" + target;
	_mkdir(dir.c_str());
	if (maximum == true && fee_and_tax == true)
		method = "\\������";
	else if (maximum == true && fee_and_tax == false)
		method = "\\�L�|�L����O";
	else if (maximum == false && fee_and_tax == true)
		method = "\\�t����";
	dir += method;
	_mkdir(dir.c_str());
	string slid_dir;
	for (int i = 0; i < 13; i++) {
		slid_dir = dir + "\\" + window[i];
		_mkdir(slid_dir.c_str());
		string mode;
		mode = slid_dir + "\\�V�m��";
		_mkdir(mode.c_str());
		mode = slid_dir + "\\���մ�";
		_mkdir(mode.c_str());
		mode = slid_dir + "\\Gbest";
		_mkdir(mode.c_str());
	}
}
int read_file(string in) {
	fstream file;
	//string input = testing_period + in;
	//string input = "�ѻ�/0050/t_0050_2010-201706.csv";
	string input = "�ѻ�/"+stock_price+"/" + sliding_window + '/' + in;
	file.open(input);
	string line;
	int c = 0;
	while (getline(file, line, '\n')) {
		istringstream templine(line);
		string data;
		int s = 0;
		while (getline(templine, data, ',')) {
			if (c == 0) stock.push_back(atoi(data.c_str()));
			else {
				dprice[s].push_back(atof(data.c_str()));
				s++;
			}
		}
		s = 0;
		c++;
	}
	file.close();
	//��Ȧ�
	//stock.push_back(9999);
	//stock.push_back(10000);
	return c - 1;
}
void force_measure(int n) {
	for (int i = 0; i < s_num; i++) {
		if (n == 0) {
			if (i == 6|| i == 13 || i == 42) {
				sol[n].push_back(1);
				portfolio[n].push_back(i);
			}
			else sol[n].push_back(0);
		}
		if (n == 1) {
			if (i == 7 || i == 42 || i == 45) {
				sol[n].push_back(1);
				portfolio[n].push_back(i);
			}
			else sol[n].push_back(0);
		}
	}
	/*M2M - 2014_07 �t����
		i == 3 || i == 5 || i == 14 || i == 15 || i == 19 || i == 22 || i == 24 || i == 26 || i == 45 || i == 49
		M2M - 2014_07 ������
		i == 14 || i == 15 || i == 22 || i == 45 || i == 49
	/*------------------------------*/
}
/* --------- ����� ---------- */
void single(int p) {
	for (int i = 0; i < s_num; i++) {
		if (i == p) {
			sol[p].push_back(1);
			portfolio[p].push_back(i);
		}
		else sol[p].push_back(0);
	}
}
/* --------------------------- */

void measure_9999(int n, int g) {
	/* ------------ ���W������̦n��C���Ĥ@�N�Ĥ@�Ӳɤl ----------- */
	if (n == 0 && g == 0) {
		force_measure(n);
	}
	else {
		/* ---------------------------------------------------------- */
		for (int i = 0; i < s_num; i++) {
			double r_num = rand() / double(RAND_MAX);

			if (r_num >= beta[i]) {

				sol[n].push_back(1);
				portfolio[n].push_back(i);
			}
			else sol[n].push_back(0);
		}
	}
	/*
	cout << "particle " << n+1 << ": ";
	for(int i=0; i<s_num; i++) cout << sol[n][i] << " ";
	cout << endl;
	*/
}
bool measure(int n, int g) {
	bool nonzero = false;

	for (int i = 0; i < s_num; i++) {
		/* ------------ ���en������̦n��Ĥ@�N�Ĥ@�Ӳɤl -----------
		if(n == 0 && g == 0 && first == false){
			sol[n].assign(guided_p.begin(),guided_p.end());
			guided_p.clear();
		}
		else{
		/* ---------------------------------------------------------- */
		double r_num = rand() / double(RAND_MAX);

		if (beta[i] > r_num) {
			sol[n].push_back(1);
			portfolio[n].push_back(i);
			nonzero = true;
		}
		else sol[n].push_back(0);
	}

	//}
	/*
	cout << "particle " << n+1 << ": ";
	for(int i=0; i<s_num; i++) cout << sol[n][i] << " ";
	cout << endl;
	*/
	return nonzero;
}

void adaptive_angle(int g) {
	/* ----------- ��ΡA�C1000�N�Τ@�Ө��� -----------
	double t = (r_upper - r_lower)/10;
	r_angle = (g/1000) * t + r_lower;
	/* ------------------------------------------------ */

	/* ------- exp�U���A�A��0-1�ϥ��W�ƹ������� ------- */
	//double x = (double)g/1000;
	double norm = exp(0 - g);
	r_angle = norm * (r_upper - r_lower) + r_lower;
	//if(g == 0 || g == 1 || g == 1000 || g == 2000 || g == 3000 || g == 4000 || g == 5000 || g == 6000 || g == 7000 || g == 8000 || g == 9000 || g == 9999)
	//	cout << "g: " << g << " /n exp(-x): " << norm << " /n r_angle: " << r_angle << endl;
	/* ------------------------------------------------ */
}



int cal_share(double funds, double dprice) {
	int share;
	
	if (fee_and_tax == true) {
		/*Origin*/
		share = floor(funds / (dprice * 1000 + dprice * 1000 * fee));
	}
	else if (fee_and_tax == false){
	/*Check-fee-remainder*/
		share = floor(funds / (dprice * 1000));
	}
	
	return share;
}
/* ---------- �R�s�� ---------- */
int cal_odds(double funds, double dprice) {
	int odds;
	if (fee_and_tax == true) {
		odds = floor(funds / (dprice * (1 + fee)));
	}
	else if(fee_and_tax == false){
		odds = floor(funds / (dprice));
	}
	return odds;
}
/* ---------------------------- */
double cal_remainder(double funds, double dprice, int share) {
	double tmp_remainder;
	if (fee_and_tax == true) {
		/*Origin*/
		if (odd_lots == false)
			tmp_remainder = funds - share * dprice * 1000 - share * dprice * 1000 * fee;
		else if (odd_lots == true)
			tmp_remainder = funds - share * dprice - share * dprice * fee;
	}
	else if (fee_and_tax == false) {
		/*Check-fee-remainder*/
		if (odd_lots == false)
			tmp_remainder = funds - share * dprice * 1000;
		else if (odd_lots == true)
			tmp_remainder = funds - share * dprice;
	}
	
	return tmp_remainder;
	
}
/* ---------- �|�^�ǬO�_�R���_ -----------
bool cal_p_info(int n, double total_fund) { //�@�}�l�粒���զX��p��ӧ��զX��T
	bool flag = true;						//�ݳo�ӧ��զX��쪺�O�_���R����
	int p_num = portfolio[n].size();
	if(p_num == 0) funds[n] = 0;
	else funds[n] = total_fund / p_num; //�����U�ɤ��쪺��
	funds_remainder[n] = total_fund - funds[n] * p_num; //������զX���l��
	for (int i = 0; i < p_num; i++) { //p_num����n�ӧ��զX���X��
		int tmp_stock = portfolio[n][i]; //��n�ӧ��զX��쪺�U�ɪѲ�
		int tmp_share = cal_share(funds[n], dprice[tmp_stock][0]); //����զX�̨C�ɪѲ��i�R�i��
		if(tmp_share == 0){ //�p�G������@�i�R����N���諸�B�z
			flag = false;
			share[n].clear();
			break;
		}
		double tmp_remainder = cal_remainder(funds[n], dprice[tmp_stock][0], tmp_share);
		share[n].push_back(tmp_share);
		s_remainder[n].push_back(tmp_remainder);
	}
	return flag;
}
/* ---------------------------------------- */

/* --------- �p���Ȧ�Q�v�A�����P��������P%�� -------- */
double cal_bank_interest(double principal) {
	/* - ���j�Ȧ業���w�s�Q�v�A�x�ȡB�X�@�B�Ĥ@�B�ثn�B�g�a -
	double interest[5] = {0.0005, 0.0016, 0.0039, 0.0104};
	/* ---�l���w�s�Q�v----- 1�Ӥ�,    3�Ӥ�,   6�Ӥ�,  1�~ --- */
	double interest[4] = { 0.00067, 0.002175, 0.00495, 0.0128 };
	/* ------------------------------------------------------- */
	if (sliding_window == "Y2Y") {
		principal *= (1 + interest[3]);
	}
	else if (sliding_window == "Y2H") {
		if (train_period == true)
			principal *= (1 + interest[3]);
		else
			principal *= (1 + interest[2]);
	}
	else if (sliding_window == "Y2Q") {
		if (train_period == true)
			principal *= (1 + interest[3]);
		else
			principal *= (1 + interest[1]);
	}
	else if (sliding_window == "Y2M") {
		if (train_period == true)
			principal *= (1 + interest[3]);
		else
			principal *= (1 + interest[0]);
	}
	else if (sliding_window == "H2H") {
		principal *= (1 + interest[2]);
	}
	else if (sliding_window == "H2Q") {
		if (train_period == true)
			principal *= (1 + interest[2]);
		else
			principal *= (1 + interest[1]);
	}
	else if (sliding_window == "H2M") {
		if (train_period == true)
			principal *= (1 + interest[2]);
		else
			principal *= (1 + interest[0]);
	}
	else if (sliding_window == "H#") {
		principal *= (1 + interest[2]);
	}
	else if (sliding_window == "Q2Q") {
		principal *= (1 + interest[1]);
	}
	else if (sliding_window == "Q2M") {
		if (train_period == true)
			principal *= (1 + interest[1]);
		else
			principal *= (1 + interest[0]);
	}
	else if (sliding_window == "Q#") {
		principal *= (1 + interest[1]);
	}
	else if (sliding_window == "M2M") {
		principal *= (1 + interest[0]);
	}
	else if (sliding_window == "M#") {
		principal *= (1 + interest[0]);
	}

	return principal;
}
/* ------------------------------------------------------ */
void cal_p_info(int n, double total_fund) { //�@�}�l�粒���զX��p��ӧ��զX��T
	int p_num = portfolio[n].size();
	if (p_num == 0) funds[n] = 0;
	else funds[n] = total_fund / p_num; //�����U�ɤ��쪺��
	funds_remainder[n] = total_fund - funds[n] * p_num; //������զX���l��
	for (int i = 0; i < p_num; i++) { //p_num����n�ӧ��զX���X��
		int tmp_stock = portfolio[n][i]; //��n�ӧ��զX��쪺�U�ɪѲ�
		//cout << "Tmp stock: " << tmp_stock << " ";
		/* ------------ ����Ȧ� ----------- */
		if (mode == 1 && tmp_stock == s_num - 1) {
			share[n].push_back(-1);
			s_remainder[n].push_back(funds[n]);
		}
		/* ----------------------------------- */
		else {
			int tmp_share;
			if (odd_lots == false) {
				tmp_share = cal_share(funds[n], dprice[tmp_stock][0]); //����զX�̨C�ɪѲ��i�R�i��
			}
			else if (odd_lots == true) {
				tmp_share = cal_odds(funds[n], dprice[tmp_stock][0]);//����զX�̨C�ɪѲ��i�R�Ѽ�
			}
			double tmp_remainder = cal_remainder(funds[n], dprice[tmp_stock][0], tmp_share);
			share[n].push_back(tmp_share);
			s_remainder[n].push_back(tmp_remainder);
		}
	}
}

double funds_standardization(int n, int s, int d) {
	double tmp_fs = 0;
	int tmp_stock = portfolio[n][s];

	/* ------------ ����Ȧ� ----------- */
	if (mode == 1 && tmp_stock == s_num - 1) {
		/* ----- �̫�@�Ѻ�Q�� ----- */
		if (d == day - 1) {
			tmp_fs = cal_bank_interest(funds[n]);
		}
		else tmp_fs = funds[n];
		/* -------------------------- */
		/* ----- �ĤG�Ѻ�Q�� -------
		if(d == 0) tmp_fs = funds[n];
		else tmp_fs = cal_bank_interest(funds[n]);
		/* -------------------------- */
		/* ----- �Ĥ@�Ѻ�Q�� -------
		tmp_fs = cal_bank_interest(funds[n]);
		/* -------------------------- */
	}
	/* ----------------------------------- */
	else {
		double tmp_price = dprice[tmp_stock][d];
		int tmp_share = share[n][s];
		int p_num = portfolio[n].size();
		/* ------------------------ �̤p��쬰�i�� -------------------------- */
		if (odd_lots == false) {
			if(maximum == false){
			/* ---- ��Ӥ@��������������k  ---- */
				if (d == 0) { 
					if (fee_and_tax == true) {
						/*Origin*/
						tmp_fs += funds[n] - tmp_share * tmp_price * fee * 1000;
					}
					else if(fee_and_tax == false){
						/*Check-fee*/
						tmp_fs += funds[n];
					}
					
					/*Check-fee-remainder
					tmp_fs += funds[n] - s_remainder[n][s];
					/*Check-remainder
					tmp_fs += funds[n] - tmp_share * tmp_price * fee * 1000 - s_remainder[n][s];
					*/
					
					
				}
				else {
					if (fee_and_tax == true) {
						/*Origin*/
						tmp_fs += tmp_share * tmp_price * 1000 - tmp_share * tmp_price * 1000 * fee
							- tmp_share * tmp_price * 1000 * stt + s_remainder[n][s];
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs += tmp_share * tmp_price * 1000 + s_remainder[n][s];
					}
					/*Check-fee-stt-remainder
					tmp_fs += tmp_share * tmp_price * 1000;
					

					/*Check-remainder
					tmp_fs += tmp_share * tmp_price * 1000 - tmp_share * tmp_price * 1000 * fee
						- tmp_share * tmp_price * 1000 * stt;
						*/
				}
			}
			/* -------------------- END ----------------------*/

			/* ---- �Ĩ��������k�A�B�ର������ ---- */
			else if(maximum == true){
				if (d == 0) { //�Ĩ�Ĥ@�Ѭ���X�Ѳ����������
					first_day_price.push_back(tmp_price);

					if (fee_and_tax == true) {
						/*Origin*/
						double handling = tmp_share * tmp_price * (fee + stt) * 1000;
						if (funds[n] >= handling) {
							tmp_fs = funds[n] - handling;
						}
						else tmp_fs = funds[n];
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = funds[n];
					}
					/*Check-fee-remainder
					tmp_fs = funds[n] - s_remainder[n][s];
					

					/*Check-remainder
					double handling = tmp_share * tmp_price * (fee + stt) * 1000;
					if (funds[n] >= handling) {
						tmp_fs = funds[n] - handling - s_remainder[n][s];
					}
					else tmp_fs = funds[n] - s_remainder[n][s];
					*/
					
				}
				else {//�ĤG�ѫ�Y�R�^�A�Ȫ��N�O��Ĥ@�Ѫѻ������t�A����N�O�Ĥ@�Ѹ������+�Ȫ����t-����O
					if (fee_and_tax == true) {
						/*Origin*/
						tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price) * 1000 - tmp_share * tmp_price * 1000 * fee;
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price) * 1000;
					}
					/*Check-fee-remainder
					tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price) * 1000;
					
					/*Check-remainder
					tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price) * 1000 - tmp_share * tmp_price * 1000 * fee; 
					*/
					

				}
			}
			/* -------------------- END ----------------------*/
		}
		/* ------------------------------------------------------------------ */
		/* ------------------------ �̤p��쬰�Ѽ� -------------------------- */
		else if (odd_lots == true) {
			/*------------��̤p�A��@������k�ۦP�A���R���--------*/
			if(maximum == false){
				if (d == 0) {
					if (fee_and_tax == true) {
						/*Origin*/
						double handling = tmp_share * tmp_price * fee;
						first_day_price.push_back(tmp_price);
						if (funds[n] >= handling) tmp_fs = funds[n] - handling;
						else tmp_fs = funds[n];
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = funds[n];
					}
				}
				else {
					if (fee_and_tax == true) {
						/*Origin*/
						tmp_fs = tmp_share * tmp_price - tmp_share * tmp_price * (fee+stt) + s_remainder[n][s];
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = tmp_share * tmp_price + s_remainder[n][s];

					}
				}
			}
			/*------------------------ END -------------------------*/
			/* ---- �Ĩ��������k�A�B�ର������ ---- */
			else if(maximum == true){
				if (d == 0) { //�Ĩ�Ĥ@�Ѭ���X�Ѳ����������
					if (fee_and_tax == true) {
						/*Origin*/
						double handling = tmp_share * tmp_price * (fee + stt);
						first_day_price.push_back(tmp_price);
						if (funds[n] >= handling) tmp_fs = funds[n] - handling;
						else tmp_fs = funds[n];
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = funds[n];
					}
				}
				else {//�ĤG�ѫ�Y�R�^�A�Ȫ��N�O��Ĥ@�Ѫѻ������t�A����N�O�Ĥ@�Ѹ������+�Ȫ����t-����O
					if (fee_and_tax == true) {
						/*Origin*/
						tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price) - tmp_share * tmp_price * fee;
					}
					else if (fee_and_tax == false) {
						/*Check-fee-stt*/
						tmp_fs = seperate_fs[n][s] + tmp_share * (first_day_price[s] - tmp_price);
					}
				}
			}
			/* -------------------- END ----------------------*/
		}
		/* ------------------------------------------------------------------ */
	}
	return tmp_fs;
}

void daily_fs(int n, int day) {
	int p_num = portfolio[n].size();
	/*Origin*/
	double portfolio_fs = funds_remainder[n];
	/*Check-fee-remainder
	double portfolio_fs = 0/*funds_remainder[n]*/;
	/*Check-fee
	double portfolio_fs = funds_remainder[n];
	/*Check-remainder
	double portfolio_fs = 0;
	*/
	for (int i = 0; i < p_num; i++) {
		double tmp_fs = funds_standardization(n, i, day);
		seperate_fs[n].push_back(tmp_fs);
		portfolio_fs += tmp_fs;
	}
	fs[n].push_back(portfolio_fs);
}

double expected_return(int n, double init_money) {
	double tmp_exp;
	double num = 0, denom = 0;
	for (int i = 1; i <= day; i++) {
		double tmp_num = i * fs[n][i - 1] - i * init_money;
		double tmp_denom = i * i;
		num += tmp_num;
		denom += tmp_denom;
	}
	tmp_exp = num / denom;
	return tmp_exp;
}

double ii_slope(int n) {
	double tmp_m;
	double num = 0, num2 = 0, num3 = 0, denom = 0;
	for (int i = 1; i <= day; i++) {
		double tmp_num = i * fs[n][i - 1];
		double tmp_denom = i * i;
		num += tmp_num;
		num2 += i;
		num3 += fs[n][i - 1];
		denom += tmp_denom;
	}
	num = day * num - (num2 * num3);
	denom = day * denom - num2 * num2;
	tmp_m = num / denom;
	return tmp_m;
}

double ii_constant(int n) {
	double tmp_c;
	double f_sum = 0;
	double d_sum = 0;
	for (int i = 1; i <= day; i++) {
		f_sum += fs[n][i - 1];
		d_sum += i;
	}
	tmp_c = f_sum / day - exp_r[n] * (d_sum / day);
	return tmp_c;
}

double cal_exp_fs(int n, int d, int sw, double init_money) {
	double d_fs;
	if (sw == 0)
		d_fs = exp_r[n] * d + init_money;
	else if (sw == 1)
		d_fs = exp_r[n] * d + constant[n];
	return d_fs;
}

double cal_risk(int n, int sw, double init_money) {
	double differ_sum = 0;
	for (int j = 0; j < day; j++) {
		double d_fs = cal_exp_fs(n, j + 1, sw, init_money);
		double differ = pow((d_fs - fs[n][j]), 2);
		differ_sum += differ;
	}
	double tmp_risk;
	tmp_risk = sqrt(differ_sum / day);
	return tmp_risk;
}

double real_reward(int n, double init_money) {
	double r = fs[n][day - 1] - init_money;
	return r;
}

void sharpe(vector<double> fs, double init_money) {
	double tmp = 0;
	int d = fs.size();
	if (d == 0) {
		sha_rew = 0;
		sha_risk = 0;
		sharpe_ratio = 0;
	}
	else {
		for (int i = 0; i < d; i++)
			tmp += fs[i];
		double avg = tmp / d;
		sha_rew = (fs[d - 1] - init_money) / init_money - rfr;
		tmp = 0;
		for (int i = 0; i < d; i++)
			tmp += pow(fs[i] - avg, 2);
		sha_risk = sqrt(tmp / d) / avg;
		sharpe_ratio = sha_rew / sha_risk;
	}
}

double fitness(double exp_r, double risk) {
	double f;
	if(maximum == false) {
		if (exp_r < 0)
			f = exp_r / risk;
		else
			f = exp_r * risk;
	}
	else if (maximum == true){
		if (exp_r > 0)
			f = exp_r / risk;
		else
			f = exp_r * risk;
	}
	return f;
}

void record(int best) {
	Gbest_share.clear();
	Gbest_s_remainder.clear();
	Gbest_fs.clear();
	Gbest_portfolio.clear();
	for (int i = 0; i < max_s; i++) Gbest_seperate_fs[i].clear();

	Gbest_reward = reward[best];
	Gbest_fund = funds[best];
	Gbest_risk = risk[best];
	Gbest_exp_r = exp_r[best];

	int n = portfolio[best].size();
	for (int i = 0; i < n; i++) {
		int t_share = share[best][i];
		double t_s_remainder = s_remainder[best][i];
		int t_portfolio = portfolio[best][i];
		Gbest_share.push_back(t_share);
		Gbest_s_remainder.push_back(t_s_remainder);
		Gbest_portfolio.push_back(t_portfolio);
		for (int d = 0; d < day; d++) {
			double t_s_fs = seperate_fs[best][d*n + i];
			Gbest_seperate_fs[i].push_back(t_s_fs);
		}
	}
	for (int d = 0; d < day; d++) {
		double t_fs = fs[best][d];
		Gbest_fs.push_back(t_fs);
	}
}

void record_worst(int worst) {
	/*
	Gworst_fs.clear();
	Gworst_portfolio.clear();

	int n = portfolio[worst].size();
	for(int i=0; i<n; i++){
		int t_portfolio = portfolio[worst][i];
		Gworst_portfolio.push_back(t_portfolio);
	}
	for(int d=0; d<day; d++){
		double t_fs = fs[worst][d];
		Gworst_fs.push_back(t_fs);
	}
	*/
}

void update_all_Gbest(int t) {
	all_best = G_best;

	all_Gbest_p.clear();
	all_Gbest_share.clear();
	all_Gbest_s_remainder.clear();
	all_Gbest_fs.clear();
	all_Gbest_portfolio.clear();
	for (int i = 0; i < max_s; i++) all_Gbest_seperate_fs[i].clear();

	all_Gbest_p.assign(Gbest_p.begin(), Gbest_p.end());
	all_Gbest_share.assign(Gbest_share.begin(), Gbest_share.end());

	all_Gbest_s_remainder.assign(Gbest_s_remainder.begin(), Gbest_s_remainder.end());
	all_Gbest_fs.assign(Gbest_fs.begin(), Gbest_fs.end());
	all_Gbest_portfolio.assign(Gbest_portfolio.begin(), Gbest_portfolio.end());
	for (int i = 0; i < max_s; i++) all_Gbest_seperate_fs[i].assign(Gbest_seperate_fs[i].begin(), Gbest_seperate_fs[i].end());

	all_Gbest_times = t;
	all_Gbest_gen = Gbest_gen;
	all_Gbest_reward = Gbest_reward;
	all_Gbest_fund = Gbest_fund;
	all_Gbest_risk = Gbest_risk;
	all_Gbest_exp_r = Gbest_exp_r;
	sharpe(all_Gbest_fs, init_fund);
}
/*
void update_all_Gworst(){
		all_worst = G_worst;

		all_Gworst_p.clear();
		all_Gworst_fs.clear();
		all_Gworst_portfolio.clear();

		all_Gworst_p.assign(Gworst_p.begin(),Gworst_p.end());

		all_Gworst_fs.assign(Gworst_fs.begin(),Gworst_fs.end());
		all_Gworst_portfolio.assign(Gworst_portfolio.begin(),Gworst_portfolio.end());

		all_Gworst_gen = Gworst_gen;
}
*/
void ALL_G(int t) {

	if (G_best == all_best) find_times++;
	else if ((maximum == true && G_best > all_best) || (maximum == false && G_best < all_best)) find_times = 1;

	if ((maximum == true  && G_best > all_best) || (maximum == false  && G_best < all_best)) {
		update_all_Gbest(t);
	}
	/*
	if(all_worst > G_worst){
		update_all_Gworst();
	}
	*/
}

void entanglement(vector<int> best) {
	int len = best.size();
	sol[N + 1].assign(best.begin(), best.end());
	/*
	cout << "best_p: ";
	for(int z=0; z<len; z++) cout << sol[best_p][z] << " ";
	cout << endl;
	cout << "beta: ";
	for(int z=0; z<len; z++) cout << beta[z] << " ";
	cout << endl;
	cout << "Gbest: " << G_best << endl;
	*/
	for (int i = 0; i < len; i++) {
		if (best[i] == 1) {
			//cout << "i: " << i+1 << endl;
			for (int j = 0; j < len; j++) {
				if (best[j] == 0 && beta[j] < 0.5) {
					sol[N + 1][i] = 0;
					sol[N + 1][j] = 1;
					/*
					cout << "test_p: ";
					for(int z=0; z<len; z++) cout << sol[N+1][z] << " ";
					cout << endl;
					*/
					for (int k = 0; k < len; k++) {
						if (sol[N + 1][k] == 1)
							portfolio[N + 1].push_back(k);
					}
					cal_p_info(N + 1, init_fund);
					for (int d = 0; d < day; d++) daily_fs(N + 1, d);
					exp_r[N + 1] = expected_return(N + 1, init_fund);
					risk[N + 1] = cal_risk(N + 1, 0, init_fund);	// 0:�_�I�� , 1:�Ͷխ�
					fit[N + 1] = fitness(exp_r[N + 1], risk[N + 1]);
					reward[N + 1] = real_reward(N + 1, init_fund);
					if (fit[N + 1] > G_best) {
						G_best = fit[N + 1];
						Gbest_p.clear();
						Gbest_p.assign(sol[N + 1].begin(), sol[N + 1].end());
						record(N + 1);
					}
					/*
					cout << "fit[N+1]: " << fit[N+1] << endl;
					cout << "Gbest: " << G_best << endl;
					system("pause");
					*/
					portfolio[N + 1].clear();
					share[N + 1].clear();
					fs[N + 1].clear();
					s_remainder[N + 1].clear();
					seperate_fs[N + 1].clear();
					sol[N + 1][i] = 1;
					sol[N + 1][j] = 0;
				}
			}
		}
	}
}


void update(int gen) {
	P_best = fit[0];
	P_worst = fit[0];
	best_p = 0;
	worst_p = 0;

	for (int i = 1; i < N; i++) {
		if(maximum == true){
			if (fit[i] > P_best) {
				P_best = fit[i];
				best_p = i;
			}
			if (fit[i] < P_worst) {
				P_worst = fit[i];
				worst_p = i;
			}
		}
		else if (maximum == false){
			if (fit[i] < P_best) {
				P_best = fit[i];
				best_p = i;
			}
			if (fit[i] > P_worst) {
				P_worst = fit[i];
				worst_p = i;
			}
		}
	}

	/* ------ ��s�Ĥ@�N��쪺Pworst��Gworst ----- */
	if (gen == 1) {
		G_worst = P_worst;
		Gworst_gen = gen;
		Gworst_p.clear();
		Gworst_p.assign(sol[worst_p].begin(), sol[worst_p].end());
		record_worst(worst_p);
	}
	/* -------------------------------------------- */
	/*
	if(gen == 1){
		if(P_best > 0){
			G_best = P_best;
			Gbest_gen = gen;
			Gbest_p.clear();
			for(int i=0; i<s_num; i++){
				Gbest_p.push_back(sol[best_p][i]);
			}
			record(best_p);
		}
		else if(G_best == 0 && P_best <= 0){
			Gbest_p.clear();
			for(int i=0; i<s_num; i++){
				Gbest_p.push_back(0);
			}
		}
		if(P_worst < 0){
			G_worst = P_worst;
			Gworst_gen = gen;
			Gworst_p.clear();
			for(int i=0; i<s_num; i++){
				Gworst_p.push_back(sol[worst_p][i]);
			}
			record_worst(worst_p);
		}
		else if(G_worst == 0 && P_worst >= 0){
			Gworst_p.clear();
			for(int i=0; i<s_num; i++){
				Gworst_p.push_back(0);
			}
		}
	}
	*/
	//��flag = 1, beta��s�覡���P
	if ((maximum == true && P_best > G_best) || (maximum == false && P_best < G_best)) {
		//beta�Ȥ��઺����
		//flag = 1;
		/* ��sbeta���P�_
		for(int i=0; i<s_num; i++){
			if(sol[best_p][i] > Gbest_p[i]){
				if(beta[i] > 0.5) beta[i] = 1-beta[i];
			}
			else if(sol[best_p][i] < Gbest_p[i]){
				if(beta[i] < 0.5) beta[i] = 1-beta[i];
			}
		}
		*/
		//���`��sGbest
		G_best = P_best;
		Gbest_gen = gen;
		Gbest_p.clear();
		Gbest_p.assign(sol[best_p].begin(), sol[best_p].end());
		record(best_p);
	}

	if ((maximum == true && P_worst < G_worst) || (maximum == false && P_worst > G_worst)) {
		G_worst = P_worst;
		Gworst_gen = gen;
		Gworst_p.clear();
		for (int i = 0; i < s_num; i++) {
			Gworst_p.push_back(sol[worst_p][i]);
		}
		record_worst(worst_p);
	}

	// Gbest�BGworst�ɤ�
	/*
	for (int i = 0; i<s_num; i++) {
		if (Gbest_p[i] > Gworst_p[i])
			beta[i] -= r_angle;
		else if (Gbest_p[i] < Gworst_p[i])
			beta[i] += r_angle;
	}
	*/

	// Gbest�BPworst�ɤ�

	for (int i = 0; i < s_num; i++) {
		if (Gbest_p[i] > sol[worst_p][i]) {
			//if(flag == 1 && beta[i] > 0.5){
			/* ----- NQTS ------ */
			if (beta[i] < 0.5) {
				beta[i] = 1 - beta[i];
				//	flag = 0;
			}
			/* ----------------- */
			beta[i] += r_angle;
		}
		else if (Gbest_p[i] < sol[worst_p][i]) {
			//if(flag == 1 && beta[i] < 0.5){
			/* ----- NQTS ------ */
			if (beta[i] > 0.5) {
				beta[i] = 1 - beta[i];
				//	flag = 0;
			}
			/* ----------------- */
			beta[i] -= r_angle;
		}
	}
	/*
	if(output_test == true){
		cout << "UPDATE_TRUE" << endl;
		fstream file;
		string out = dir + "/" + sliding_window + "/Update_solution_" + file_name;
		cout << file_name << endl;
		cout << out << endl;
		file.open(out, ios::app);
		file << "Gen" << "," << gen << endl;
		file << "Gbest" << "," << G_best << endl;
		file << "Gbest_sol" << endl;
		for(int i=0; i<s_num; i++)
			file << Gbest_p[i] << ",";
		file << endl;
		file << "Worst" << "," << setprecision(16) << P_worst << endl;
		file << "Worst_sol" << endl;
		for(int i=0; i<s_num; i++)
			file << sol[worst_p][i] << ",";
		file << endl;
		file.close();
	}
	*/
	

	//Pbest�BPworst�ɤ�
	/*
	for (int i = 0; i<s_num; i++) {
		if (sol[best_p][i] > sol[worst_p][i])
			beta[i] -= r_angle;
		else if (sol[best_p][i] < sol[worst_p][i])
			beta[i] += r_angle;
	}
	*/
}

/* ------- ��l�̨θѬ���Ȧ� -------- */
void init() {
	s_num = stock.size();  //�s�W1bit -> ��Ȧ�
	if (mode == 0) {
		/* --------- ��lGbest��0 -------------- */
		for (int i = 0; i < s_num; i++) {
			Gbest_p.push_back(0);
			beta.push_back(p);
		}
		exp_r[N] = 0;
		risk[N] = 0;
		fit[N] = 0;
		reward[N] = 0;
		G_best = 0;
		for (int d = 0; d < day; d++) {
			seperate_fs[N].push_back(init_fund);
			fs[N].push_back(init_fund);
		}
		record(N);
		/* ------------------------------------- */
	}
	else if (mode == 1) {
		/* --------- ��lGbest����Ȧ� --------- */
		for (int i = 0; i < s_num; i++) {
			if (i == s_num - 1) {
				sol[N].push_back(1);
				portfolio[N].push_back(i);
			}
			else sol[N].push_back(0);
			beta.push_back(p);
		}
		cal_p_info(N, init_fund);
		for (int d = 0; d < day; d++) daily_fs(N, d);
		exp_r[N] = expected_return(N, init_fund);
		risk[N] = cal_risk(N, 0, init_fund);	// 0:�_�I�� , 1:�Ͷխ�
		fit[N] = fitness(exp_r[N], risk[N]);
		reward[N] = real_reward(N, init_fund);
		G_best = fit[N];
		Gbest_gen = 0;
		Gbest_p.assign(sol[N].begin(), sol[N].end());
		record(N);
		/* -------------------------------------- */
	}
}

void test_portfolio(vector<int> train_stock) {
	for (int i = 0; i < train_stock.size(); i++) {
		portfolio[0].push_back(train_stock[i]);
	}
	/* ----0716���շQ�k----
	cal_p_info(0,invest_fund);
	exp_r[0] = all_Gbest_exp_r;
	for(int d=0; d<day; d++) daily_fs(0,d);
	risk[0] = cal_risk(0,0,init_fund);
	fit[0] = fitness(exp_r[0], risk[0]);
	reward[0] = real_reward(0,invest_fund);
	 -------------------- */
}

double total_expected_return(double init_money) {
	double tmp_exp;
	double num = 0, denom = 0;
	for (int i = 1; i <= total_test_day; i++) {
		double tmp_num = i * total_test_fs[i - 1] - i * init_money;
		double tmp_denom = i * i;
		num += tmp_num;
		denom += tmp_denom;
	}
	tmp_exp = num / denom;
	return tmp_exp;
}

double cal_total_exp_fs(int d, double init_money) {
	double d_fs;
	d_fs = total_test_exp_r * d + init_money;
	return d_fs;
}

double cal_total_risk(double init_money) {
	double differ_sum = 0;
	for (int j = 0; j < total_test_day; j++) {
		double d_fs = cal_total_exp_fs(j + 1, init_money);
		double differ = pow((d_fs - total_test_fs[j]), 2);
		differ_sum += differ;
	}
	double tmp_risk;
	tmp_risk = sqrt(differ_sum / total_test_day);
	return tmp_risk;
}

double total_reward(double init_money) {
	double r = total_test_fs[total_test_day - 1] - init_money;
	return r;
}

void write_All_G(string file_name) {

	fstream file;
	//string out = dir + "/" + file_name;
	string out = dir + "/" + sliding_window + "/�V�m��/" + file_name;
	file.open(out, ios::out);
	file << "�N��" << "," << generation << "\n";
	file << "�ɤl��" << "," << N << "\n";
	file << "���ਤ�פW��" << "," << r_upper << "\n";
	file << "���ਤ�פU��" << "," << r_lower << "\n";
	file << "���ਤ��" << "," << r_angle << "\n";
	file << "���禸��" << "," << T << "\n";
	file << "\n";
	file << "��l���" << "," << fixed << setprecision(12) << init_fund << "\n";
	file << "�̫���" << "," << all_Gbest_fs[day - 1] << "\n";
	file << "�u����S" << "," << all_Gbest_reward << "\n";
	file << "\n";
	file << "�L���ȳ��S" << "," << sha_rew << endl;
	file << "�L���ȭ��I" << "," << sha_risk << endl;
	file << "�L����" << "," << sharpe_ratio << endl;
	file << "\n";
	file << "�w�����S" << "," << all_Gbest_exp_r << "\n";
	file << "���I" << "," << all_Gbest_risk << "\n";
	file << "�_�I��" << "," << all_best << "\n";
	file << "���̨θѥ@�N" << "," << all_Gbest_gen << "\n";
	file << "���̨θѹ���#" << "," << all_Gbest_times << "\n";
	file << "���̨θѦ���" << "," << find_times << "\n";
	file << "\n";
	file << "���զX�ɼ�" << "," << all_Gbest_portfolio.size() << "\n";
	int f = 0;
	for (int i = 0; i < all_Gbest_portfolio.size(); i++) {
		int n = all_Gbest_portfolio[i];
		file << stock[n] << "(" << n << ")" << ",";
		if (i == all_Gbest_portfolio.size() - 1) {
			file << "\n";
			if (f == 0) {
				i = -1;
				file << "Stock#" << ",";
				f = 1;
			}
		}
	}
	f = 0;

	for (int r = 0; r < day + 3; r++) {
		for (int i = 0; i <= all_Gbest_portfolio.size(); i++) {
			if (r == 0 && i == 0) file << "�i��" << ",";
			else if (r == 1 && i == 0) file << "���t���" << ",";
			else if (r == 2 && i == 0) file << "�Ѿl���" << ",";
			else if (i == 0) file << "FS( " << r - 2 << " ),";
			else if (r == 0) file << all_Gbest_share[i - 1] << ",";
			else if (r == 1) file << all_Gbest_fund << ",";
			else if (r == 2) file << all_Gbest_s_remainder[i - 1] << ",";
			else {
				if (i == all_Gbest_portfolio.size())
					file << all_Gbest_seperate_fs[i - 1][r - 3] << "," << all_Gbest_fs[r - 3];
				else
					file << all_Gbest_seperate_fs[i - 1][r - 3] << ",";
			}
		}
		file << "\n";

	}

	file.close();
}


void write_simple() {
	fstream file;
	//string out = train_simple_output;
	string out = dir + "/" + sliding_window + "/train_Gbest_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + '_' + sliding_window + "_all.csv";
	file.open(out, ios::app);

	file << period << "," << all_Gbest_portfolio.size() << ",";
	for (int i = 0; i < all_Gbest_portfolio.size(); i++) {
		int n = all_Gbest_portfolio[i];
		file << stock[n] << "(" << n << ") ";
	}
	file << "," << fixed << setprecision(10) << "�_�I��" << "," << all_best << "," << "�w�����S" << "," << all_Gbest_exp_r << "," << "���I" << "," << all_Gbest_risk;
	file << "," << "�̨θѹ���#" << "," << all_Gbest_times;
	file << "," << "�̨θѥ@�N��" << "," << all_Gbest_gen;
	file << "," << "�̨θѥX�{����" << "," << find_times << "\n";
	file.close();
}

void write_test_result(string file_name) {

	fstream file;
	//string out = test_output + file_name;
	string out = dir + "/" + sliding_window + "/���մ�/" + file_name;
	file.open(out, ios::out);
	file << "�N��" << "," << generation << "\n";
	file << "�ɤl��" << "," << N << "\n";
	file << "���ਤ��" << "," << r_angle << "\n";
	file << "���禸��" << "," << T << "\n";
	file << "\n";
	file << "��l���" << "," << fixed << setprecision(10) << invest_fund << "\n";
	file << "�̫���" << "," << fs[0][day - 1] << "\n";
	file << "�u����S" << "," << reward[0] << "\n";
	file << "\n";
	file << "�L���ȳ��S" << "," << sha_rew << endl;
	file << "�L���ȭ��I" << "," << sha_risk << endl;
	file << "�L����" << "," << sharpe_ratio << endl;
	file << "\n";
	file << "�w�����S" << "," << exp_r[0] << "\n";
	file << "���I" << "," << risk[0] << "\n";
	file << "�_�I��" << "," << fit[0] << "\n";
	file << "\n";
	file << "���զX�ɼ�" << "," << portfolio[0].size() << "\n";
	for (int i = 0; i <= portfolio[0].size(); i++) {
		if (i == 0)file << "Stock#" << ",";
		else {
			int n = portfolio[0][i - 1];
			file << stock[n] << "(" << n << ")";
			if (i == portfolio[0].size()) file << "\n";
			else file << ",";
		}
	}
	int a = portfolio[0].size();
	for (int r = 0; r < day + 3; r++) {
		for (int i = 0; i <= a; i++) {
			if (r == 0 && i == 0) file << "�i��" << ",";
			else if (r == 1 && i == 0) file << "���t���" << ",";
			else if (r == 2 && i == 0) file << "�Ѿl���" << ",";
			else if (i == 0) file << "FS( " << r - 2 << " ),";
			else if (r == 0) file << share[0][i - 1] << ",";
			else if (r == 1) file << funds[0] << ",";
			else if (r == 2) file << s_remainder[0][i - 1] << ",";
			else {
				if (i == a)
					file << seperate_fs[0][(r - 3)*a + (i - 1)] << "," << fs[0][r - 3];
				else
					file << seperate_fs[0][(r - 3)*a + (i - 1)] << ",";
			}
		}
		file << "\n";

	}

	file.close();
}

void write_test_simple() {
	fstream file;
	//string out = test_simple_output;
	string out = dir + "/" + sliding_window + "/test_Gbest_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + '_' + sliding_window + "_all.csv";
	file.open(out, ios::app);

	file << period << "," << portfolio[0].size() << ",";
	for (int i = 0; i < portfolio[0].size(); i++) {
		int n = portfolio[0][i];
		file << stock[n] << "(" << n << ") ";
	}
	file << "," << fixed << setprecision(10) << "�_�I��" << "," << fit[0] << "," << "�w�����S" << "," << exp_r[0] << "," << "���I" << "," << risk[0] << "\n";
	file.close();

}

void write_total_test() {
	fstream file;
	//string out = whole_period_result;
	string out = dir + "/" + sliding_window + "/total_test_result_" + sliding_window + ".csv";
	file.open(out, ios::out);

	file << "���մ��϶�" << "," << testing << "\n";
	file << "�@�N��" << "," << generation << "\n";
	file << "���ਤ��" << "," << r_angle << "\n";
	file << "�ɤl��" << "," << N << "\n";
	file << "���禸��" << "," << T << "\n";
	file << "��l���" << "," << fixed << setprecision(10) << init_fund << "\n";
	file << "�u����S" << "," << total_test_reward << "\n";
	file << "\n";
	file << "�L���ȳ��S" << "," << sha_rew << endl;
	file << "�L���ȭ��I" << "," << sha_risk << endl;
	file << "�L����" << "," << sharpe_ratio << endl;
	file << "\n";
	file << "�w�����S" << "," << total_test_exp_r << "\n";
	file << "���I" << "," << total_test_risk << "\n";
	file << "�_�I��" << "," << total_test_fit << "\n";
	file << "*�o�̪��_�I�ȬO�ιw�����S���W���I" << "\n";
	file << "\n";
	file << "�`�Ѽ�" << "," << total_test_day << "\n";
	file << "�Ѽ�" << "," << "�`�������" << "\n";
	for (int i = 0; i < total_test_day; i++) {
		file << "FS(" << i + 1 << ")" << "," << total_test_fs[i] << "\n";
	}

	file.close();
}

string year[10] = { "2009","2010","2011","2012","2013","2014","2015","2016","2017", "2018" };
string mon[12] = { "01","02","03","04","05","06","07","08","09","10","11","12" };
string quarter[4] = { "Q1","Q2","Q3","Q4" };
string half[2] = { "Q1-Q2","Q3-Q4" };

/* ----------------------------------- */
string file_name(int a, int b) {
	string f_name;
	if (sliding_window == "Y2Y") {
		if (train_period == true)
			f_name = "train_" + year[b] + '(' + year[b] + " Q1).csv";
		else {
			f_name = "test_";
			int tty = b + 1;
			f_name += year[tty] + '(' + year[b] + " Q1).csv";
		}
	}
	else if (sliding_window == "Y2H") {
		if (train_period == true) {
			if (b == 0)
				f_name = "train_" + year[a] + '(' + year[a] + " Q1).csv";
			else if (b == 1)
				f_name = "train_" + year[a] + "_Q3~" + year[a + 1] + "_Q2(" + year[a] + " Q1).csv";
		}
		else {
			if (b == 0)
				f_name = "test_" + year[a + 1] + "_Q1-Q2(" + year[a] + " Q1).csv";
			else if (b == 1)
				f_name = "test_" + year[a + 1] + "_Q3-Q4(" + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "Y2Q") {
		if (train_period == true) {
			if (b == 0)
				f_name = "train_" + year[a] + '(' + year[a] + " Q1).csv";
			else {
				int eqq = b - 1;
				int eyy = a + 1;
				f_name = "train_" + year[a] + '_' + quarter[b] + '~' + year[eyy] + '_' + quarter[eqq] + '(' + year[a] + " Q1).csv";
			}
		}
		else {
			f_name = "test_" + year[a + 1] + '_' + quarter[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "Y2M") {
		if (train_period == true) {
			f_name = "train_";
			if (b == 0)
				f_name += year[a] + '(' + year[a] + " Q1).csv";
			else {
				int tyy = a + 1;
				int emm = b - 1;
				f_name += year[a] + '_' + mon[b] + '~' + year[tyy] + '_' + mon[emm] + '(' + year[a] + " Q1).csv";
			}
		}
		else {
			f_name = "test_";
			int tyy = a + 1;
			f_name += year[tyy] + '_' + mon[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "H2H") {
		if (train_period == true)
			f_name = "train_" + year[a] + '_' + half[b] + '(' + year[a] + " Q1).csv";
		else {
			f_name = "test_";
			if (b == 0)
				f_name += year[a] + '_' + half[1 - b] + '(' + year[a] + " Q1).csv";
			else if (b == 1)
				f_name += year[a + 1] + '_' + half[1 - b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "H2Q") {
		if (train_period == true) {
			f_name = "train_";
			if (b == 3)
				f_name += year[a] + '_' + quarter[b] + '~' + year[a + 1] + '_' + quarter[0] + '(' + year[a] + " Q1).csv";
			else
				f_name += year[a] + '_' + quarter[b] + '-' + quarter[b + 1] + '(' + year[a] + " Q1).csv";
		}
		else {
			f_name = "test_";
			b += 2;
			if (b >= 4) { b -= 4; a += 1; }
			if (b == 0 || b == 1)
				f_name += year[a] + '_' + quarter[b] + '(' + year[a - 1] + " Q1).csv";
			else
				f_name += year[a] + '_' + quarter[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "H2M") {
		if (train_period == true) {
			f_name = "train_";
			int emm = b + 5;
			if (emm >= 12) {
				int eyy = a;
				emm -= 12;
				eyy += 1;
				f_name += year[a] + '_' + mon[b] + '~' + year[eyy] + '_' + mon[emm] + '(' + year[a] + " Q1).csv";
			}
			else
				f_name += year[a] + '_' + mon[b] + '-' + mon[emm] + '(' + year[a] + " Q1).csv";
		}
		else {
			f_name = "test_";
			b += 6;
			if (b >= 12) { b -= 12; a += 1; }
			if (b >= 0 && b <= 5)
				f_name += year[a] + '_' + mon[b] + '(' + year[a - 1] + " Q1).csv";
			else
				f_name += year[a] + '_' + mon[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "H#") {
		if (train_period == true)
			f_name = "train_" + year[a] + '_' + half[b] + '(' + year[a] + " Q1).csv";
		else {
			a += 1;
			f_name = "test_";
			f_name += year[a] + '_' + half[b] + '(' + year[a - 1] + " Q1).csv";
		}
	}
	else if (sliding_window == "Q2Q") {
		if (train_period == true)
			f_name = "train_" + year[a] + '_' + quarter[b] + '(' + year[a] + " Q1).csv";
		else {
			f_name = "test_";
			b += 1;
			if (b == 4) { b = 0; a += 1; }
			int tty;
			if (b == 0) tty = a - 1;
			else tty = a;
			f_name += year[a] + '_' + quarter[b] + '(' + year[tty] + " Q1).csv";
		}
	}
	else if (sliding_window == "Q2M") {
		if (train_period == true) {
			f_name = "train_";
			int tmm = b + 2;
			int tyy = a;
			if (tmm >= 12) {
				tmm -= 12;
				tyy += 1;
				f_name += year[a] + '_' + mon[b] + '~' + year[tyy] + '_' + mon[tmm] + '(' + year[a] + " Q1).csv";
			}
			else
				f_name += year[a] + '_' + mon[b] + '-' + mon[tmm] + '(' + year[a] + " Q1).csv";
		}
		else {
			f_name = "test_";
			b += 3;
			if (b >= 12) { b -= 12; a += 1; }
			if (b <= 2) {
				f_name += year[a] + '_' + mon[b] + '(' + year[a - 1] + " Q1).csv";
			}
			else f_name += year[a] + '_' + mon[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "Q#") {
		if (train_period == true)
			f_name = "train_" + year[a] + '_' + quarter[b] + '(' + year[a] + " Q1).csv";
		else {
			f_name = "test_";
			f_name += year[a + 1] + '_' + quarter[b] + '(' + year[a] + " Q1).csv";
		}
	}
	else if (sliding_window == "M2M") {
		int tyy = a;
		if (train_period == true)
			f_name = "train_";
		else {
			b += 1;
			if (b == 12) { b = 0; a += 1; }
			f_name = "test_";
			if (b == 0) tyy = a - 1;
		}
		f_name += year[a] + '_' + mon[b] + '(' + year[tyy] + " Q1).csv";
	}
	else if (sliding_window == "M#") {
		int tyy = a;
		if (train_period == true) {
			f_name = "train_";
		}
		else {
			f_name = "test_";
			a += 1;
			tyy = a - 1;
		}
		f_name += year[a] + '_' + mon[b] + '(' + year[tyy] + " Q1).csv";
	}
	return f_name;
}

string output_name(string in) {
	string f_name = "result_" + in;
	return f_name;
}

void test_measure(int n) {
	sol[n].assign(sol[n - 10].begin(), sol[n - 10].end());
	portfolio[n].assign(portfolio[n - 10].begin(), portfolio[n - 10].end());
	sol[n].push_back(1);
	portfolio[n].push_back(s_num - 1);

	cal_p_info(n, init_fund);
	for (int d = 0; d < day; d++) daily_fs(n, d);
	exp_r[n] = expected_return(n, init_fund);
	risk[n] = cal_risk(n, 0, init_fund);	// 0:�_�I�� , 1:�Ͷխ�
	fit[n] = fitness(exp_r[n], risk[n]);
	reward[n] = real_reward(n, init_fund);
	cout << "Portfolio: " << endl;
	for (int q = 0; q < portfolio[n].size(); q++)
		cout << portfolio[n][q] << " ";
	cout << endl;
	cout << "init_fund: " << init_fund << endl;
	cout << "Exp: " << fixed << setprecision(13) << exp_r[n] << endl;
	cout << "Risk: " << risk[n] << endl;
	cout << "Reward: " << reward[n] << endl;
	cout << "Train Fitness: " << fit[n] << endl;

}

/* ------------------ �L�Y�N��i�Ӳɤl���զX��T ------------------ */
void test_write(int i) {

	fstream file;
	string out = "test/Q2M_2015_06-08_���P�ѼơB�i�Ƭ��Ͷ�.csv";
	file.open(out, ios::app);
	file << "init_fund" << "," << fixed << setprecision(13) << init_fund << "\n";
	file << "Expected Return" << "," << exp_r[i] << "\n";
	file << "Risk" << "," << risk[i] << "\n";
	file << "Real Reward" << "," << reward[i] << "\n";
	file << "Fitness" << "," << fit[i] << "\n";

	file << "���զX�ɼ�" << "," << portfolio[0].size() << "\n";
	file << "Portfolio " << i + 1 << "\n";
	for (int q = 0; q < portfolio[i].size(); q++)
		file << portfolio[i][q] << ",";
	file << "\n";


	for (int s = 0; s <= portfolio[0].size(); s++) {
		if (s == 0)file << "Stock#" << ",";
		else {
			int n = portfolio[i][s - 1];
			file << stock[n] << "(" << n << ")";
			if (i == portfolio[i].size()) file << "\n";
			else file << ",";
		}
	}
	file << "\n";
	int a = portfolio[i].size();
	for (int r = 0; r < day + 3; r++) {
		for (int s = 0; s <= a; s++) {
			if (r == 0 && s == 0) file << "�i��" << ",";
			else if (r == 1 && s == 0) file << "���t���" << ",";
			else if (r == 2 && s == 0) file << "�Ѿl���" << ",";
			else if (s == 0) file << "FS( " << r - 2 << " ),";
			else if (r == 0) file << share[i][s - 1] << ",";
			else if (r == 1) file << funds[i] << ",";
			else if (r == 2) file << s_remainder[i][s - 1] << ",";
			else {
				if (s == a)
					file << seperate_fs[i][(r - 3)*a + (s - 1)] << "," << fs[i][r - 3];
				else
					file << seperate_fs[i][(r - 3)*a + (s - 1)] << ",";
			}
		}
		file << "\n";

	}


	file.close();
}
/* ------------------------------------------------------------------ */

int flag = true;
int gen = 0;
int times = 0;

void write() {
	fstream file;
	string out = "(srand 114)Gbest_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + ".csv";
	file.open(out, ios::app);

	gen++;
	if (flag == true) {
		file << "Generation" << "," << gen << "\n";
		file << "Best(" << best_p + 1 << ')' << "," << fit[best_p] << "," << "�w�����S" << "," << exp_r[best_p] << "," << "���I" << "," << risk[best_p] << "\n";
		file << "Worst(" << worst_p + 1 << ')' << "," << fit[worst_p] << "," << "�w�����S" << "," << exp_r[worst_p] << "," << "���I" << "," << risk[worst_p] << "\n";
		file << "Gbest(" << Gbest_gen << ')' << "," << G_best << "," << "�w�����S" << "," << Gbest_exp_r << "," << "���I" << "," << Gbest_risk << "\n";
		file << "Best_p" << "\n";
		for (int i = 0; i < s_num - 1; i++)
			file << sol[best_p][i] << ",";
		file << "\n";
		file << "Worst_p" << "\n";
		for (int i = 0; i < s_num - 1; i++)
			file << sol[worst_p][i] << ",";
		file << "\n";
		file << "Gbest_p" << "\n";
		for (int i = 0; i < s_num - 1; i++)
			file << Gbest_p[i] << ",";
		file << "\n";
		file << "Beta" << "\n";
		for (int i = 0; i < s_num - 1; i++)
			file << beta[i] << ",";
		file << "\n";
		file << "\n";
	}
	if (flag == false) {
		times++;
		file << times << "��������̦n" << "\n";
		file << "�̨θѦb��" << all_Gbest_times << "������" << "\n";
		file << "�w�����S" << "," << all_Gbest_exp_r << "\n";
		file << "���I" << "," << all_Gbest_risk << "\n";
		file << "�Ͷխ�" << "," << all_best << "\n";
		file << "\n";
		flag = true;
		gen = 0;
	}
	file.close();
}

void write_fixed_detail(int p, string file_name) {

	fstream file;
	//string out = train_output + file_name;
	//string out = dir + "/" + sliding_window + '/' + convert_str(stock[p]) + '_' + convert_str(p) + '_' + file_name;
	string out = dir + "/" + sliding_window + "/particle" + convert_str(p) + '_' + file_name;
	file.open(out, ios::out);
	file << "�N��" << "," << generation << "\n";
	file << "�ɤl��" << "," << N << "\n";
	file << "���ਤ��" << "," << r_angle << "\n";
	file << "���禸��" << "," << T << "\n";
	file << "\n";
	file << "��l���" << "," << fixed << setprecision(12) << init_fund << "\n";
	file << "�̫���" << "," << fs[p][day - 1] << "\n";
	file << "�u����S" << "," << reward[p] << "\n";
	file << "\n";
	file << "�L���ȳ��S" << "," << sha_rew << endl;
	file << "�L���ȭ��I" << "," << sha_risk << endl;
	file << "�L����" << "," << sharpe_ratio << endl;
	file << "\n";
	file << "�w�����S" << "," << exp_r[p] << "\n";
	file << "���I" << "," << risk[p] << "\n";
	file << "�_�I��" << "," << fit[p] << "\n";
	file << "���զX�ɼ�" << "," << portfolio[p].size() << "\n";
	int f = 0;
	for (int i = 0; i < portfolio[p].size(); i++) {
		int n = portfolio[p][i];
		file << stock[n] << "(" << n << ")" << ",";
		if (i == portfolio[p].size() - 1) {
			file << "\n";
			if (f == 0) {
				i = -1;
				file << "Stock#" << ",";
				f = 1;
			}
		}
	}
	f = 0;

	int a = portfolio[p].size();
	for (int r = 0; r < day + 3; r++) {
		for (int i = 0; i <= a; i++) {
			if (r == 0 && i == 0) file << "�i��" << ",";
			else if (r == 1 && i == 0) file << "���t���" << ",";
			else if (r == 2 && i == 0) file << "�Ѿl���" << ",";
			else if (i == 0) file << "FS( " << r - 2 << " ),";
			else if (r == 0) file << share[p][i - 1] << ",";
			else if (r == 1) file << funds[p] << ",";
			else if (r == 2) file << s_remainder[p][i - 1] << ",";
			else {
				if (i == a)
					file << seperate_fs[p][(r - 3)*a + (i - 1)] << "," << fs[p][r - 3];
				else
					file << seperate_fs[p][(r - 3)*a + (i - 1)] << ",";
			}
		}
		file << "\n";

	}
	file.close();
}

void write_beta(int t, int g, string input) {
	int s = input.length();
	string f_name = "";
	for (int i = 0; i < s - 4; i++) f_name += input[i];
	fstream file;
	string out = dir + "/" + sliding_window + "/" + sliding_window + '_' + f_name + "_beta.csv";
	file.open(out, ios::app);
	file << "����#" << "," << t << ",";
	file << "�@�N#" << "," << g << "\n";
	for (int i = 0; i < s_num; i++) {
		file << beta[i] << ",";
	}
	file << "\n";
	file.close();
}
/* ----- ��X�C������̨θ� ----- */
void write_Gbest(int t, string input) {
	int s = input.length();
	string f_name = "";
	for (int i = 0; i < s - 4; i++) f_name += input[i];
	fstream file;
	//string out = "test/" +sliding_window+'_'+f_name+"_Gbest.csv";
	string out = dir + "/" + sliding_window + "/Gbest/" + sliding_window + '_' + f_name + "_Gbest.csv";
	file.open(out, ios::app);
	file << "����#" << "," << t << ",";
	file << "Gbest_GEN" << "," << Gbest_gen << ",";
	file << "Gbest" << "," << G_best << ",";
	file << "Portfolio" << ",";
	for (int i = 0; i < Gbest_portfolio.size(); i++) {
		file << stock[Gbest_portfolio[i]] << "(" << Gbest_portfolio[i] << ")" << " ";
	}
	file << "\n";
	file.close();
}

/* ----- ��X�Y������C�N�̨θ� ----- */
void write_every_gen_Gbest(int t, int gen, string input) {
	int s = input.length();
	string f_name = "";
	for (int i = 0; i < s - 4; i++) f_name += input[i];
	fstream file;
	//string out = "test/" +sliding_window+'_'+f_name+"_Gbest.csv";
	string out = dir + "/" + sliding_window + "/Gbest/" + sliding_window + '_' + f_name + "exp_" + convert_str(t) + "_every_gen_Gbest.csv";
	file.open(out, ios::app);
	file << "Gen" << "," << gen << ",";
	file << "Gbest_GEN" << "," << Gbest_gen << ",";
	file << "Gbest" << "," << G_best << ",";
	file << "Portfolio" << ",";
	for (int i = 0; i < Gbest_portfolio.size(); i++) {
		file << stock[Gbest_portfolio[i]] << "(" << Gbest_portfolio[i] << ")" << " ";
	}
	file << "\n";
	file.close();
}

/* ------ ��X�Y�϶��C������̨θ� ----- */
void write_every_exp_Gbest(int t, string input) {
	int s = input.length();
	string f_name = "";
	for (int i = 0; i < s - 4; i++) f_name += input[i];
	fstream file;
	//string out = "test/" +sliding_window+'_'+f_name+"_Gbest.csv";
	string out = dir + "/" + sliding_window + "/Gbest/" + sliding_window + '_' + f_name + "_every_exp_Gbest.csv";
	file.open(out, ios::app);
	file << "Experiment" << "," << t << ",";
	file << "Gbest_GEN" << "," << Gbest_gen << ",";
	file << "Gbest" << "," << G_best << ",";
	file << "Portfolio" << ",";
	for (int i = 0; i < Gbest_portfolio.size(); i++) {
		file << stock[Gbest_portfolio[i]] << "(" << Gbest_portfolio[i] << ")" << " ";
	}
	file << "\n";
	file.close();
}

void write_all_train_simple() {
	fstream file;
	//string out = train_simple_output;
	string out = dir + "/train_Gbest_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + "_all.csv";
	file.open(out, ios::app);
	if (period == 1) file << sliding_window << "\n";
	file << period << "," << all_Gbest_portfolio.size() << ",";
	for (int i = 0; i < all_Gbest_portfolio.size(); i++) {
		int n = all_Gbest_portfolio[i];
		file << stock[n] << "(" << n << ") ";
	}
	file << "," << fixed << setprecision(10) << "�_�I��" << "," << all_best << "," << "�w�����S" << "," << all_Gbest_exp_r << "," << "���I" << "," << all_Gbest_risk;
	file << "," << "�̨θѹ���#" << "," << all_Gbest_times;
	file << "," << "�̨θѥ@�N��" << "," << all_Gbest_gen;
	file << "," << "�̨θѥX�{����" << "," << find_times;
	file << "," << "�̨θѥ���" << "," << avg_Gbest;
	file << "," << "�̨θѼзǮt" << "," << std_Gbest;
	file << "," << "�C��������̨θѥ����@�N��" << "," << avg_Gbest_gen;
	if (train_afford == false) {
		file << "," << period << "," << "**���R���_���Ѳ�";
	}
	file << "\n";
	file.close();
}

/* ---------- ��X�Ҧ��ưʵ������մ����զX²�n��T�� ------------ */
void write_all_test_simple() {
	fstream file;
	//string out = test_simple_output;
	string out = dir + "/test_Gbest_" + convert_str(generation) + '_' + convert_str(N) + '_' + convert_str(T) + '_' + convert_str(r_angle) + "_all.csv";
	file.open(out, ios::app);
	if (period == 1) file << sliding_window << "\n";
	file << period << "," << portfolio[0].size() << ",";
	for (int i = 0; i < portfolio[0].size(); i++) {
		int n = portfolio[0][i];
		file << stock[n] << "(" << n << ") ";
	}
	file << "," << fixed << setprecision(10) << "�_�I��" << "," << fit[0] << "," << "�w�����S" << "," << exp_r[0] << "," << "���I" << "," << risk[0];
	if (test_afford == false) file << "," << period << "," << "**���R���_���Ѳ�";
	file << "\n";
	file.close();

}


int a, b;
int break_a, break_b;

void sliding() {
		int y_a;
		if(stock_price.length() == 11)
			y_a = stock_price[stock_price.length() - 3] - '0' +1;
		else if (stock_price.length() == 9)
			y_a = stock_price[stock_price.length() - 1] - '0' + 1;

		period = 1;
		if (sliding_window == "Y2Y") {
			a = 1; b = y_a;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "Y2H") {
			a = y_a; b = 2;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "Y2Q") {
			a = y_a; b = 4;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "Y2M") {
			a = y_a; b = 12;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "H2H") {
			a = y_a+1; b = 2;
			break_a = y_a; break_b = 1;
		}
		else if (sliding_window == "H2Q") {
			a = y_a+1; b = 4;
			break_a = y_a; break_b = 2;
		}
		else if (sliding_window == "H2M") {
			a = y_a+1; b = 12;
			break_a = y_a; break_b = 6;
		}
		else if (sliding_window == "H#") {
			a = y_a; b = 2;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "Q2Q") {
			a = y_a+1; b = 4;
			break_a = y_a; break_b = 3;
		}
		else if (sliding_window == "Q2M") {
			a = y_a+1; b = 12;
			break_a = y_a; break_b = 9;
		}
		else if (sliding_window == "Q#") {
			a = y_a; b = 4;
			break_a = a; break_b = b;
		}
		else if (sliding_window == "M2M") {
			a = y_a+1; b = 12;
			break_a = y_a; break_b = 11;
		}
		else if (sliding_window == "M#") {
			a = y_a; b = 12;
			break_a = a; break_b = b;
		}
}

/* ---------- �p��Gbest�����B�зǮt ----------- */
void cal_Gbest_info() {
	double sum = 0;
	for (int i = 0; i < T; i++)
		sum += all_the_Gbest[i];
	avg_Gbest = sum / T;
	double tmp = 0;
	for (int i = 0; i < T; i++)
		tmp += pow((all_the_Gbest[i] - avg_Gbest), 2);
	std_Gbest = sqrt(tmp / T);
}

int main() {
	double start_time = clock();
	//srand(time(NULL));
	make_file();
	/*	�@pe[13] = { Y2Y(0), Y2H(1), Y2Q(2), Y2M(3), H2H(4), H2Q(5), H2M(6), H*(7), Q2Q(8), Q2M(9), Q*(10), M2M(11), M*(12)};�@*/
	for (int pe = 0; pe < 13; pe++) {
		srand(114);
		run_out = false;
		//first = true;
		sliding_window = window[pe];
		cout << sliding_window << endl;
		sliding();
		invest_fund = init_fund;
		for (int yy =0; yy < a; yy++) {
			for (int pp = 0; pp < b; pp++) {
				if (yy == break_a && pp == break_b) break;
				if (sliding_window == "H2H") {
					if (yy == 0 && pp == 0) { pp = 1; }
				}
				else if (sliding_window == "H2Q") {
					if (yy == 0 && pp == 0) { pp = 2; }
				}
				else if (sliding_window == "H2M") {
					if (yy == 0 && pp == 0) { pp = 6; }
				}
				else if (sliding_window == "Q2Q") {
					if (yy == 0 && pp == 0) { pp = 3; }
				}
				else if (sliding_window == "Q2M") {
					if (yy == 0 && pp == 0) { pp = 9; }
				}
				else if (sliding_window == "M2M") {
					if (yy == 0 && pp == 0) { pp = 11; }
				}
				train_period = true;
				string input_file = file_name(yy, pp);
				//string input_file = "0050.csv";
				cout << input_file << endl;

				/* ------------------�V�m��----------------- */

				day = read_file(input_file);
				find_times = 0;
				//stock.push_back(9999);
		/* ------------------------------------------ */
				for (int t = 0; t < T; t++) {
					/*
					if(yy == 5 && pp == 6 && t == 0) output_test = true;
					else output_test = false;
					*/
					init();
					
					for (int g = 0; g < generation; g++) {
						//if(run_out == true) break;
						//adaptive_angle(g);
						for (int i = 0; i < N; i++) {
							first_day_price.clear();
							bool nonzero = true;
							nonzero = measure(i, g);
							//if(output_test == false) continue;
							//force_measure(i);
							//single(i);
							cal_p_info(i, init_fund);
							for (int d = 0; d < day; d++) daily_fs(i, d);

							if (nonzero == false) {
								exp_r[i] = 0;
								risk[i] = 0;
								fit[i] = 0;
								reward[i] = 0;
							}
							else {
								exp_r[i] = expected_return(i, init_fund);
								risk[i] = cal_risk(i, 0, init_fund);	// 0:�_�I�� , 1:�Ͷխ�
								fit[i] = fitness(exp_r[i], risk[i]);
								reward[i] = real_reward(i, init_fund);
							}
							/* ------- ��X�C�Ӳɤl��T ------- 
							for(int k=0; k<portfolio[i].size(); k++)
								cout << stock[portfolio[i][k]] << ",  (" << portfolio[i][k] << ") , ";
							cout << endl;
							cout << "fitness: "<< fit[i] << endl;
							cout << "exp_r: "<< exp_r[i] << endl;
							cout << "risk:"<< risk[i] << endl;
							system("PAUSE");
							/*
							write_fixed_detail(i, input_file);
							fstream file;
							string out = dir + "/" + sliding_window + "/single_stock_" + input_file;
							file.open(out, ios::app);
							for(int k=0; k<portfolio[i].size(); k++)
								file << stock[portfolio[i][k]] << "," << portfolio[i][k];
							file << endl;
							file << "fitness" << "," << fit[i] << endl;
							file << "exp_r" << "," << exp_r[i] << endl;
							file << "risk" << "," << risk[i] << endl;
							file.close();
							/* ---------------------------------- */
							
							
							/*
							if(g == 9999){
							for(int k=0; k<portfolio[i].size(); k++)
								cout << portfolio[i][k] << " ";
							cout << endl;
							cout << "fitness: " << fit[i] << endl;
							cout << "exp_r: " << exp_r[i] << endl;
							cout << "risk: " << risk[i] << endl;
							}
							*/
						}
						/* --------- ��s�̨θѤ�beta�}�C ---------- */
						update(g + 1);
						/*
						if(output_test == true){
							write_every_gen_Gbest(t+1, g+1, input_file);
							write_beta(t+1,g+1,input_file);
						}
						*/
						/* ----------------------------------------- */
						/* -------------�M����N�V�m���------------ */
						for (int i = 0; i <= N; i++) {
							portfolio[i].clear();
							share[i].clear();
							s_remainder[i].clear();
							fs[i].clear();
							sol[i].clear();
							seperate_fs[i].clear();
						}
						/* ----------------------------------------- */
					}
					/*
					if(output_test == true)
						write_Gbest(t+1,input_file);
					*/
					ALL_G(t + 1); //��s���v�H�ӳ̨θ�
					/*
					cout << "����#" << t+1 << endl;
					cout << "Gbest: " << G_best << endl;
					cout << "Gbest_Gen: " << Gbest_gen << endl;
					cout << endl;
					cout << "�̨θ�" << endl;
					cout << "BEST_EXP#" << all_Gbest_times << endl;
					cout << "BEST_GEN: " << all_Gbest_gen << endl;
					cout << "BEST: " << all_best << endl;
					cout << endl;
					*/
					/* -------------�M���Ӧ��V�m���------------ */
					G_best = 0;
					Gbest_gen = 0;
					beta.clear();
					Gbest_p.clear();
					Gbest_fs.clear();
					Gbest_share.clear();
					Gbest_portfolio.clear();
					Gbest_s_remainder.clear();
					for (int i = 0; i < max_s; i++) Gbest_seperate_fs[i].clear();

					G_worst = 0;
					Gworst_p.clear();
					Gworst_fs.clear();
					Gworst_portfolio.clear();

					sharpe_ratio = 0;
					sha_rew = 0;
					sha_risk = 0;
					/* ----------------------------------------- */

				}
				/* ------------------------------------------------- */
				avg_Gbest_gen /= T;							//�p����Gbest�������@�N��
				cal_Gbest_info();								//�p��Gbest�Ѫ������B�зǮt
				/* --------- �ݦ��S���R���_���Ѳ� ---------- */
				for (int q = 0; q < all_Gbest_share.size(); q++) {
					if (all_Gbest_share[q] == 0)
						train_afford = false;
				}
				/* ------------��X�V�m���ɮ�--------------- */
				string out = output_name(input_file);
				if (all_Gbest_portfolio.size() != 0) {
					write_All_G(out);
				}
				write_simple();
				write_all_train_simple();
				/* ----------------------------------------- */


				/* -----------�M���V�m���ѻ����------------ */
				stock.clear();
				for (int i = 0; i < max_s; i++) dprice[i].clear();
				/* ----------------------------------------- */

				/* ------------------���մ�----------------- */
				train_period = false;
				test_input = file_name(yy, pp);
				//test_input  = "0050_2010-201706.csv";
				cout << test_input << endl;
				if (t_str == false) {
					testing += test_input;
					t_str = true;
				}
				train_day = day;
				first_day_price.clear();
				day = read_file(test_input);
				if (all_Gbest_portfolio.size() == 0 || run_out == true) {
					for (int d = 0; d < day; d++)
						total_test_fs.push_back(invest_fund);
				}
				else {
					test_portfolio(all_Gbest_portfolio);
					cal_p_info(0, invest_fund);
					for (int d = 0; d < day; d++) {
						daily_fs(0, d);
						total_test_fs.push_back(fs[0][d]);
					}

					exp_r[0] = expected_return(0, invest_fund);
					risk[0] = cal_risk(0, 0, invest_fund);	// 0:�_�I�� , 1:�Ͷխ�
					fit[0] = fitness(exp_r[0], risk[0]);
					reward[0] = real_reward(0, invest_fund);
					sharpe(fs[0], invest_fund);

					//cout << "Test Fitness: " << fit[0] << endl;

				/* ----------------------------------------- */

					for (int q = 0; q < share[0].size(); q++) {
						if (share[0][q] == 0) test_afford = false;
					}
				/* ------------��X���մ��ɮ�--------------- */

					string test_out = output_name(test_input);
					write_test_result(test_out);
					write_test_simple();
					write_all_test_simple();

					/* ----------------------------------------- */
					invest_fund = fs[0][day - 1];	// �C�q���մ�������n���b�@�_�A�U�@�q���մ��������B�O�o�����մ��̫᪺�������
					//init_fund = fs[0][day-1];	// �V�m������l����δ��մ��W�@�q�̫�@�Ѫ�����}�l�V�m
					if(invest_fund <= 0) run_out = true;
					else run_out= false;

				}
				period++;
				/* -----------�M�����մ��Ѳ����------------ */
				seperate_fs[0].clear();
				portfolio[0].clear();
				share[0].clear();
				s_remainder[0].clear();
				fs[0].clear();

				sharpe_ratio = 0;
				sha_rew = 0;
				sha_risk = 0;
				/* ----------------------------------------- */

				/* ------------�M����e�g�����------------- */
				stock.clear();
				for (int i = 0; i < max_s; i++) dprice[i].clear();
				all_best = 0;
				all_Gbest_p.clear();
				all_Gbest_fs.clear();
				all_Gbest_share.clear();
				all_Gbest_portfolio.clear();
				all_Gbest_s_remainder.clear();
				for (int i = 0; i < max_s; i++) all_Gbest_seperate_fs[i].clear();
				/*
				all_worst = 0;
				all_Gworst_p.clear();
				all_Gworst_fs.clear();
				all_Gworst_portfolio.clear();
				*/
				/* ----------------------------------------- */
			}
		}
		/* ----------- �p���q���մ��϶� ------------ */

		total_test_day = total_test_fs.size();
		total_test_exp_r = total_expected_return(init_fund);
		total_test_risk = cal_total_risk(init_fund);
		total_test_reward = total_reward(init_fund);
		total_test_fit = fitness(total_test_exp_r, total_test_risk);
		sharpe(total_test_fs, init_fund);
		testing += " - " + test_input;
		write_total_test();
		/* ------------------------------------------- */
		total_test_fs.clear();
		testing = "";
		t_str = false;
	}
	double end_time = clock();
	double time = (end_time - start_time) / CLOCKS_PER_SEC;
	if(time < 60)
		cout << time << "��" << endl;
	else{
		int hour = time / 3600;
		int min  = (time - hour * 3600) / 60;
		int sec = time - hour * 3600 - min * 60;
		cout << hour << " �p�� " << min << " �� " << sec << " ��" << endl;
	}

	system("pause");
	return 0;
}