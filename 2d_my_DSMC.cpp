// myDSMC.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//
//質量だけ無次元化してない
#define _USE_MATH_DEFINES // for C 
//#define M_PI 3.1415926535

#include "stdafx.h"
#include <stdio.h>
#include <cmath>
#include <random>
#include <iostream> //入出力ライブラリ
#include <fstream> //iostreamのファイル入出力をサポート
#include <time.h>     // for clock()
#include <vector>       // ヘッダファイルインクルード

#include "mydef.h"		//define を集めてある場所

/*-------------------------------------------------関数の宣言--------------*/
//int body_collision(int mol, double dtfly, double part_x, double part_y, double part_cx, double part_cy, double part_hit[7]);
int particle_injection(double xbuf1, double xbuf2, double ybuf1, double ybuf2, double part_xy[][2], double part_c[][3]);
int re_indexing(double part_xy[][2], double part_c[][3], int lc[][my], int lcr[nmax], int  lc0[][my], int i_j_cel[][2]);
int collision(double part_xy[][2], double part_c[][3], double col[mx][my], int lc[][my], int lcr[nmax], int lc0[][my]);
double uniform_random();
double min_thit(double thit[4]);
double injection_update(double xzone1, double xzone2, double yzone1, double yzone2, double part_xy[][2], double part_c[][3]);
double body_force_and_heat(double part_hit[], double c_x, double c_y, double xbody_force_x[][B_my], double xbody_force_y[][B_my], double xbody_energy[][B_my], double ybody_force_x[][B_mx], double ybody_force_y[][B_mx], double ybody_energy[][B_mx]);
double line_body_collsion(int mol, double body_point_pair[][2][2], double dtfly, double part_x, double part_y, double part_cx, double part_cy, double part_hit[]);
double rotation(double theta, double xy[], double after_xy[]);
double angle_vectors(double vect_x, double vect_y);
double body_rotate(double theta);
double delete_particle_in_body();
double calc_torque(double hit_x, double hit_y, double u, double v);
double intersec_lines(double a_bef[2], double b_bef[2], double a_aft[2], double b_aft[2], double c[2]);
double near_point_on_line(double a_bef[2], double b_bef[2], double d[2]);
double rotate_collsion_line(double dtheta, double sum_torque);
double c_new_generate(double a_dash[2], double b_dash[2], double theta_sign, double v_part_new[2]);
double check_particle_in_body();
/*-------------------------------------------------------------------------*/
//物体の回転はあり？なし？
const int flag_body_rotate = 1;

//Flow Condition
const double rgas = 287.0;
const double v0 = 8000.0;
const double t0 = 200.0;
const double d = 0.1;
const double twall = 200.0;
const double akn = 0.4;
const double vref = sqrt(2.0*rgas*t0);
const double tref = d / vref;
const double r_uni = 8.314;						//一般気体定数
const double mol_mass = 29.0;
//高度300km
const double p0 = 8.7704e-6;						//3.2011*10^-2  Pa
const double rho = 1.916e-11;

////計算領域
//const double xmin = -5.0;
//const double xmax = 5.0;
//const double ymin = -5.0;
//const double ymax = 5.0;
//const double zmin = -5.0;
//const double zmax = 5.0;

//計算領域
const double xmin = -0.5;
const double xmax = 0.5;
const double ymin = -0.5;
const double ymax = 0.5;
const double zmin = -0.5;
const double zmax = 0.5;

////物体の大きさ
//const double xbody1 = -0.5;
//const double xbody2 = 0.5;
//const double ybody1 = -0.5;
//const double ybody2 = 0.5;

//物体の大きさ
const double xbody1 = -0.15 ;
const double xbody2 = 0.15 ;
const double ybody1 = -0.05;
const double ybody2 = 0.05;

const int sides = 8;				//物体の辺の数

//double body_point[sides][2][2] = { { { xbody1,ybody1 },{ xbody1,ybody2 } },
//{ { xbody1,ybody2 },{ xbody2,ybody2 } },
//{ { xbody2,ybody2 },{ xbody2,ybody1 } },
//{ { xbody2,ybody1 },{ xbody1,ybody1 } } };

//double body_point[sides][2][2] = { { { -1.5,0.5 },{ 1.45,0.5 } },
//{ { 1.45,0.5 },{ 1.45,3.5 } },
//{ { 1.45,3.5 },{ 1.5,3.5 } },
//{ {1.5,3.5},{1.5,-3.5} },
//{ { 1.5,-3.5 },{ 1.45,-3.5 } },
//{ { 1.45,-3.5 },{ 1.45,-0.5 } },
//{ {1.45,-0.5},{-1.5,-0.5} },
//{ {-1.5,-0.5},{-1.5,0.5} } };
double body_point[sides][2][2] = { { { -0.15,0.05 },{ 0.145,0.05 } },
{ { 0.145,0.05 },{ 0.145,0.35 } },
{ { 0.145,0.35 },{ 0.15,0.35 } },
{ { 0.15,0.35 },{ 0.15,-0.35 } },
{ { 0.15,-0.35 },{ 0.145,-0.35 } },
{ { 0.145,-0.35 },{ 0.145,-0.05 } },
{ { 0.145,-0.05 },{ -0.15,-0.05 } },
{ { -0.15,-0.05 },{ -0.15,0.05 } } };



//物体の慣性モーメント
const double I = 1.0e-11;							//質量の無次元化はする？

													//計算領域の詳細確定
const double dx = (xmax - xmin) / double(mx);
const double dy = (ymax - ymin) / double(my);
const double body_dx = (xbody2 - xbody1) / double(B_mx);
const double body_dy = (ybody2 - ybody1) / double(B_my);

//繰り返し回数
const int nlast = 5000;

//サンプリングするステップ数
const int nlp = 100;

//set n*
const int ns = 20000;

//dtを決定する
const double dtv = (d*fmin(dx, dy) / (v0 + sqrt(2.0*rgas*t0))) / tref;
const double dtc = 0.2*akn;
//二つのうち小さい方を採用
const double dt = fmin(dtv, dtc);

//Gmaxを定義する
const double fgmax = 10.0;
const double gmax = fgmax*sqrt(rgas*t0) / vref;

//分子の質量
const double mass_particle = p0*d*d*1.0 / (r_uni*t0)*mol_mass*0.001 / ns;				//単位はkg
const double mass_part_2 = rho*d*d / ns;												//無次元化はする？？

																						//領域内にある粒子の個数
int nmol = 0;												//<----危ないかも

															//////////////////大きな配列の置き場所/////////////////////////

double(*particle_x_y)[2] = new double[nmax][2];				//粒子の位置 [x座標][y座標]		ただし、[0]の値は使わないとする。
double(*particle_x_y_new)[2] = new double[nmax][2];
double(*particle_c)[3] = new double[nmax][3];				//粒子の速度 [速度u][速度v]		上と同じ
double(*particle_c_new)[3] = new double[nmax][3];
double(*grid_density)[my] = new double[mx][my];					//グリッド内の粒子密度
double(*velocity_data) = new double[nmax];

//格子の座標を定義
double xgrid[mx + 1][my + 1] = {};
double ygrid[mx + 1][my + 1] = {};
double xcg[mx][my] = {};
double ycg[mx][my] = {};
double area[mx][my] = {};
double col[mx][my] = {};

////格子内での平均値
double rou[mdx][mdy] = {};
double u[mdx][mdy] = {};
double v[mdx][mdy] = {};
double smol[mdx][mdy] = {};
double u1[mdx][mdy] = {};
double v1[mdx][mdy] = {};
double w1[mdx][mdy] = {};
double su1[mdx][mdy] = {};
double sv1[mdx][mdy] = {};
double sw1[mdx][mdy] = {};

////Body 表面にグリッドを張る
double body_xgrid[B_mx] = {};
double body_ygrid[B_my] = {};

double xbody_force_x[2][B_my] = {};
double xbody_force_y[2][B_my] = {};
double xbody_energy[2][B_my] = {};
double ybody_force_x[2][B_mx] = {};
double ybody_force_y[2][B_mx] = {};
double ybody_energy[2][B_mx] = {};


int main()
{

	//時間計測
	clock_t start = clock();    // スタート時間

								/*(1)	Computationl Conditions*/
								//粒子の諸元を配列に

								//上の二次元配列を初期化
	for (int i = 0; i < nmax; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			particle_x_y[i][j] = 0.0;
			particle_x_y_new[i][j] = 0.0;
		}
	}

	for (int i = 0; i < nmax; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			particle_c[i][j] = 0.0;
			particle_c_new[i][j] = 0.0;
		}
	}

	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			grid_density[i][j] = 0;
		}
	}

	///////////////////大量の配列のもともと置いてあった場所///////////////////////////

	/*(2)	格子の座標値を確定 & dt*/

	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			xgrid[i][j] = xmin + dx*double(i);
			ygrid[i][j] = ymin + dy*double(j);
		}
	}

	//セル中心や面積etc.を定義
	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			xcg[i][j] = (xgrid[i][j] + xgrid[i + 1][j] + xgrid[i][j + 1] + xgrid[i + 1][j + 1]) / 4.0;
			ycg[i][j] = (ygrid[i][j] + ygrid[i + 1][j] + ygrid[i][j + 1] + ygrid[i + 1][j + 1]) / 4.0;
			area[i][j] = dx*dy;
			col[i][j] = fgmax*dt / (4.0*area[i][j] * double(ns)*akn);
		}
	}

	for (int i = 0; i < B_mx; i++)				//0.5ずれてる
	{
		body_xgrid[i] = xbody1 + body_dx*double(i);
	}
	for (int j = 0; j < B_my; j++)
	{
		body_ygrid[j] = ybody1 + body_dy*double(j);
	}


	/*(3)	計算条件を示す*/

	printf(">> Two-Dimensional Rarefied Gas Flow Simulation....\n");
	printf(">> by using DSMC with Max - Collision - Number Method..\n");
	printf(">> for 2D Flow of Monoatomic Gas\n");
	printf(">> around Rectangular Body\n ");
	printf(" \n");
	printf(">> Computational Cells.........\n");
	printf("    No. of Cell       = %d X %d\n", mx, my);
	printf("    Comp. Domain    x = %f, %f\n", xmin, xmax);
	printf("                    y = %f, %f\n", ymin, ymax);
	printf("    Body            x = %f, %f\n", xbody1, xbody2);
	printf("                    y = %f, %f\n", ybody1, ybody2);
	printf("    Cell length     dx = %f, dy = %f\n", dx, dy);
	printf(">> Flow Conditions.............\n");
	printf("    V0   (m/s)        = %f\n", v0);
	printf("    T0   (K)          = %f\n", t0);
	printf("    Kn                = %f\n", akn);
	printf("    Twall(K)          = %f\n", twall);
	printf("    particle_mass     = %e\n", mass_particle);
	printf("    particle_mass_2        = %e\n", mass_part_2);
	printf(">> Computational Conditions....\n");
	printf("    Factor for Gmax   = %f\n", fgmax);
	printf("    N*                = %d\n", ns);
	printf("    Time Step Size    = %f\n", dt);
	printf("    No. of Time Steps = %d\n", nlast);
	printf("    Vref              = %f\n", vref);
	//	printf("    Steps to reach steady state (guess) \n");
	//	printf("                      = %f\n",nstdy);
	//	printf("    Sampling Interval = %f\n",);
	//	printf("    No. of Time Steps = %f\n");
	printf("\n");

	//物体の初期位置の回転
	double init_theta = M_PI / 3.0;
	double body_omega = 0;
	double d_theta = 0;
	double body_theta = init_theta;					//物体の現在の迎角を格納

													//bodyの回転角を格納しておく場所
	std::vector<std::vector<double>> body_theta_list;								// ローカル変数として、v を生成
	std::vector<std::vector<double>> body_torque_list;
	body_theta_list = std::vector<std::vector<double>>(nlast + 1, std::vector<double>(2, 0));		//二次元配列[time][theta]
	body_torque_list = std::vector<std::vector<double>>(nlast + 1 , std::vector<double>(2, 0));			//[time][torque] 有次元化した値を格納
																									//初期値を格納
	body_theta_list[0][1] = body_theta * 180 / M_PI;

	body_rotate(init_theta);

	/*(4)	Time Marching*/

	int nstep = 1;
	for (nstep = 1; nstep <= nlast; nstep++)
	{
		double time = dt*double(nstep);

		/*物体の回転*/
		double sum_torque = 0;					//全ての粒子が物体に及ぼす力積，右周り正になっているので注意
		//printf("body_omega: %f", body_omega);

		//printf("theta; %e\n", d_theta);
		//printf("dtheta : %f\n", d_theta);
		if (flag_body_rotate == 1)
		{
			d_theta = body_omega*dt;
			body_theta += d_theta;
			sum_torque = rotate_collsion_line(d_theta, sum_torque);
			body_rotate(d_theta);
		}

//		delete_particle_in_body();
		//check_particle_in_body();


		//bodyの回転角を格納
		body_theta_list[nstep][0] = time*tref;
		body_theta_list[nstep][1] = body_theta * 180 / M_PI;

		//粒子番号の定義
		int mol = 1;
		while (mol <= nmol)
		{
			/*物体と粒子の衝突*/
			double dtfly = dt;

			particle_x_y_new[mol][0] = particle_x_y[mol][0] + particle_c[mol][0] * dtfly;
			particle_x_y_new[mol][1] = particle_x_y[mol][1] + particle_c[mol][1] * dtfly;

			double part_hit[7] = {};				//更新するときに新しい値を格納する場所

			line_body_collsion(mol, body_point, dtfly, particle_x_y[mol][0], particle_x_y[mol][1], particle_c[mol][0], particle_c[mol][1], part_hit);

			//iwcol != 0 なら残り時間分粒子を飛ばす
			while (int(part_hit[0]) != 0)
			{
				double part1_torque = 0;
				double part2_torque = 0;
				//物体にかかる力の計算
				body_force_and_heat(part_hit, particle_c[mol][0], particle_c[mol][1], xbody_force_x, xbody_force_y, xbody_energy, ybody_force_x, ybody_force_y, ybody_energy);
				part1_torque = calc_torque(part_hit[2], part_hit[3], particle_c[mol][0], particle_c[mol][1]);		//入射する粒子によるトルク
				part2_torque = -calc_torque(part_hit[2], part_hit[3], part_hit[4], part_hit[5]);						//反射する粒子によるトルク
				
				sum_torque += part1_torque + part2_torque;
				//if (std::isnan(sum_torque)) {
				//	printf("sum_nan\n");
				//	printf("mol:%d\n", mol);
				//	printf("xy: %f,%f\n", part_hit[2], part_hit[3]);
				//	printf("bef_vel: %f,%f\n", particle_c[mol][0], particle_c[mol][1]);
				//	printf("aft_vel: %f,%f\n", part_hit[4], part_hit[5]);
				//	printf("sum_torque: %e\n", sum_torque);
				//}
				//
																													//壁衝突判定のためのパラメタを元に戻す
				part_hit[0] = 0.0;

				//壁衝突の時までデータを更新
				//dtflyはtwcol分減って、残りは上と同じことをする
				//printf("%d\n", mol);

				particle_x_y[mol][0] = part_hit[2];
				particle_x_y[mol][1] = part_hit[3];
				particle_c[mol][0] = part_hit[4];
				particle_c[mol][1] = part_hit[5];
				particle_c[mol][2] = part_hit[6];
				dtfly = dtfly - part_hit[1];

				////ここから上と同じ
				particle_x_y_new[mol][0] = particle_x_y[mol][0] + particle_c[mol][0] * dtfly;
				particle_x_y_new[mol][1] = particle_x_y[mol][1] + particle_c[mol][1] * dtfly;

				line_body_collsion(mol, body_point, dtfly, particle_x_y[mol][0], particle_x_y[mol][1], particle_c[mol][0], particle_c[mol][1], part_hit);
				//ここまで
			}

			/* Delete Particles*/
			if (((particle_x_y_new[mol][0] <= xmin) || (particle_x_y_new[mol][0] >= xmax)) ||
				((particle_x_y_new[mol][1] <= ymin) || (particle_x_y_new[mol][1] >= ymax)))
			{
				particle_x_y[mol][0] = particle_x_y[nmol][0];
				particle_x_y[mol][1] = particle_x_y[nmol][1];
				particle_c[mol][0] = particle_c[nmol][0];
				particle_c[mol][1] = particle_c[nmol][1];
				particle_c[mol][2] = particle_c[nmol][2];
				nmol = nmol - 1;
				mol = mol - 1;
			}

			//粒子番号の更新は最後に
			mol += 1;
		}

		//物体回転の角速度を計算と，トルクの格納
		body_torque_list[nstep][0] = time*tref;
		body_torque_list[nstep][1] = sum_torque * vref / (dt * tref);
		
		if (flag_body_rotate == 1)
		{
			body_omega = sum_torque / (dt*I);
		}

		////物体内に粒子が入ったかチェック
		//for (int m = 1; m <= nmol; m++)
		//{
		//	if (((particle_x_y_new[m][0] > xbody1) && (particle_x_y_new[m][0] < xbody2)) &&
		//		((particle_x_y_new[m][1] > ybody1) && (particle_x_y_new[m][1] < ybody2)))// &&
		//		//((particle_x_y[m][0] < xbody1) || (particle_x_y[m][0] > xbody2) ||
		//		//(particle_x_y[m][1] < ybody1) || (particle_x_y[m][1] > ybody2)))
		//	{
		//		printf("x = %f, y = %f\n", particle_x_y[m][0], particle_x_y[m][1]);
		//		printf("new_x = %f, new_y = %f\n", particle_x_y_new[m][0], particle_x_y_new[m][1]);
		//	}
		//}

		//update position
		for (int m = 1; m <= nmol; m++)
		{
			particle_x_y[m][0] = particle_x_y_new[m][0];
			particle_x_y[m][1] = particle_x_y_new[m][1];
		}

		/*Particle Injection*/

		//バッファゾーンの指定
		double xbuf1 = xmin - dt*(4.0*sqrt(2.0*rgas*t0) + v0) / vref;
		double xbuf2 = xmax + dt*fmax(0.0, (4.0*sqrt(2.0*rgas*t0) - v0) / vref);
		double ybuf1 = xmin - dt*4.0*sqrt(2.0*rgas*t0) / vref;
		double ybuf2 = xmax + dt*4.0*sqrt(2.0*rgas*t0) / vref;

		particle_injection(xbuf1, xbuf2, ybuf1, ybuf2, particle_x_y, particle_c);

		//for (int m = 1; m <= nmol; m++)
		//{
		//	printf("%f, %f\n", particle_x_y[m][0], particle_x_y[m][1]);
		//}

		//				                                          5.4 Re - indexing particles
		/*......................................................................*/
		//     lc(i, j) : number of particles in cell(i, j)
		//     lc0(i, j) : total number of particles in cells(1, 1) to(i, j - 1)
		//     lcr(k) : relation between indices k and mol			kが新しいインデックス、mはもともとのインデックス
		//     icel(mol) : cell index i for particle mol
		//     jcel(mol) : cell index j for particle mol
		//     grid_density : 最後にグリッド内の粒子個数を保存する
		/*......................................................................*/

		int(*lc)[my] = new int[mx][my];
		int(*lc0)[my] = new int[mx][my];
		int *lcr;
		lcr = new int[nmax];
		int(*i_j_cel)[2] = new int[nmax][2];

		//初期化
		for (int i = 0; i < mx; i++)
		{
			for (int j = 0; j < my; j++)
			{
				lc[i][j] = 0;
				lc0[i][j] = 0;
			}
		}
		for (int i = 0; i < nmax; i++)
		{
			lcr[i] = 0;
		}
		for (int i = 0; i < nmax; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				i_j_cel[i][j] = 0;
			}
		}

		re_indexing(particle_x_y, particle_c, lc, lcr, lc0, i_j_cel);
		//															5.5 Collision
		collision(particle_x_y, particle_c, col, lc, lcr, lc0);

		//最後のステップで、各グリッド内の粒子の個数を保存する

		//密度を長時間にわたってはかる（時間平均を求める感じ）
		//if (nstep >= nlast - 1000)

		//密度を長時間にわたってはかる（時間平均を求める感じ）
		if ((nstep >= nlast - 500) && (nstep % 10 == 0))
		{
			for (int i = 0; i < mx; i++)
			{
				for (int j = 0; j < my; j++)
				{
					grid_density[i][j] += double(lc[i][j]) / 50.0;
				}
			}
		}

		delete[] lc;
		delete[] lc0;
		delete[] lcr;
		delete[] i_j_cel;

		//速度のデータを集める

		for (int m = 0; m < nmol; m++)
		{
			velocity_data[m] = sqrt(pow(particle_c[m][0], 2.0) + pow(particle_c[m][1], 2.0) + pow(particle_c[m][2], 2.0));
		}

		//															5.8 Monitoring
		if ((nstep%nlp) == 0)
		{
			printf(".....Step %d No. of Particle = %d\n", nstep, nmol);
			printf("                        time = %f\n", time*tref);
			printf("                     d_theta = %e\n", d_theta * 180 / M_PI);
			printf("                  body_theta = %f\n", body_theta * 180 / M_PI);
		}

	}

	//															5.9 calculate Cd
	double sum_body_momentum = 0;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < B_my; j++)
		{
			sum_body_momentum += xbody_force_x[i][j];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < B_mx; j++)
		{
			sum_body_momentum += ybody_force_x[i][j];
		}
	}
	double force = sum_body_momentum / (dt*nlast);
	double cd = force / (0.5*rho*v0 / vref*v0 / vref*d*d);			//速度を無次元化する，無次元化した値からcdを計算する．


	//トルクを時間平均をとる
	double torque_ave = 0;
	for (int n = 1; n <= nlast; n++)
	{
		torque_ave += body_torque_list[n][1] / double(nlast);
	}

																	//for (int m = 1; m <= 10; m++)
																	//{
																	//	printf("%f, %f\n", particle_x_y[m][0], particle_x_y[m][1]);
																	//}

																	//																(6) Final Result
	printf("F = %e\nCd = %f\n", force, cd);
	printf("Torque_ave = %e\n", torque_ave);

	std::ofstream ofs("Result.data", std::ios::out | std::ios::trunc);
	for (int i = 1; i <= nmol; i++)
	{
		ofs << particle_x_y[i][0] << " " << particle_x_y[i][1] << " " << std::endl;
	}

	std::ofstream ofs_density("density.data", std::ios::out | std::ios::trunc);
	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			ofs_density << xgrid[i][j] << " " << ygrid[i][j] << " " << grid_density[i][j] / (double(ns) * dx * dy) << std::endl;
		}
		ofs_density << std::endl;
	}


	std::ofstream ofs_velocity("velocity.data", std::ios::out | std::ios::trunc);
	for (int i = 1; i <= nmol; i++)
	{
		ofs_velocity << velocity_data[i] << std::endl;
	}

	std::ofstream ofs_press("body_pressure.data", std::ios::out | std::ios::trunc);
	{
		for (int i = 0; i < B_my; i++)
		{
			ofs_press << body_ygrid[i] << " " << xbody_force_x[0][i] * dt*nlast*d / vref << " " << xbody_force_y[0][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_press << std::endl;
		for (int i = 0; i < B_my; i++)
		{
			ofs_press << body_ygrid[i] << " " << xbody_force_x[1][i] * dt*nlast*d / vref << " " << xbody_force_y[1][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_press << std::endl;
		for (int i = 0; i < B_mx; i++)
		{
			ofs_press << body_xgrid[i] << " " << ybody_force_x[0][i] * dt*nlast*d / vref << " " << ybody_force_y[0][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_press << std::endl;
		for (int i = 0; i < B_mx; i++)
		{
			ofs_press << body_xgrid[i] << " " << ybody_force_x[1][i] * dt*nlast*d / vref << " " << ybody_force_y[1][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_press << std::endl;
	}

	std::ofstream ofs_energy("body_energy.data", std::ios::out | std::ios::trunc);
	{
		for (int i = 0; i < B_my; i++)
		{
			ofs_energy << body_ygrid[i] << " " << xbody_energy[0][i] *dt*nlast*d/vref << std::endl;
		}
		ofs_energy << std::endl;
		for (int i = 0; i < B_my; i++)
		{
			ofs_energy << body_ygrid[i] << " " << xbody_energy[1][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_energy << std::endl;
		for (int i = 0; i < B_mx; i++)
		{
			ofs_energy << body_xgrid[i] << " " << ybody_energy[0][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_energy << std::endl;
		for (int i = 0; i < B_mx; i++)
		{
			ofs_energy << body_xgrid[i] << " " << ybody_energy[1][i] * dt*nlast*d / vref << std::endl;
		}
		ofs_energy << std::endl;
	}

	std::ofstream ofs_theta("body_theta.data", std::ios::out | std::ios::trunc);
	for (int i = 0; i <= nlast; i++)
	{
		ofs_theta << body_theta_list[i][0] << " " << body_theta_list[i][1] << " " << std::endl;
	}

	std::ofstream ofs_etc("body_etc.data", std::ios::out | std::ios::trunc);
		ofs_etc << "force" << " " << force << std::endl;
		ofs_etc << "cd" << " " << cd  << std::endl;
		ofs_etc << "torque_ave" << torque_ave << std::endl;


	//時間計測
	clock_t end = clock();     // 終了時間
	std::cout << "main_runtime = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

	return 0;
}
//
//int body_collision(int mol, double dtfly, double part_x, double part_y, double part_cx, double part_cy, double part_hit[7])		//part_hit は更新用にmain関数内で用意
//{
//	// thit は衝突が起こるまでにかかる時間
//	double thit[4] = {};
//	thit[0] = (xbody1 - part_x) / part_cx;
//	thit[1] = (xbody2 - part_x) / part_cx;
//	thit[2] = (ybody1 - part_y) / part_cy;
//	thit[3] = (ybody2 - part_y) / part_cy;
//
//	double iwcol = 0.0;		//どちらの壁に当たったか判定
//
//	//double thit_min = min_thit(thit);		//最小のthitを求める. 但し、thitは正であるものに限る
//	double thit_min = 10000.0;				//適当に大きい数を使用
//
//											//xbodyで当たるか判定する
//
//	if (part_cx != 0.0)
//	{
//		//xbody1で当たるか判定
//		if ((thit[0] > 0.0) && (thit[0] <= dtfly))
//		{
//			double yhit1 = part_y + part_cy * thit[0];
//			if ((yhit1 >= ybody1) && (yhit1 <= ybody2))
//			{
//				if (thit[0] < thit_min)
//				{
//					thit_min = thit[0];
//					iwcol = 1.0;
//					double twhit = thit[0];
//					double xwhit = xbody1;
//					double ywhit = yhit1;
//					//壁では拡散反射モデルを採用
//					//パラメタ
//					double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double bb = 2.0 * M_PI*uniform_random();
//
//					double cxhit = -sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double cyhit = aa*cos(bb);
//					double czhit = aa*sin(bb);
//
//					//返り値の用意
//					part_hit[0] = iwcol;
//					part_hit[1] = twhit;
//					part_hit[2] = xwhit;
//					part_hit[3] = ywhit;
//					part_hit[4] = cxhit;
//					part_hit[5] = cyhit;
//					part_hit[6] = czhit;
//
//					
//				}
//			}
//		}
//	}
//
//	//xbody2で当たるか判定
//	if (part_cx != 0.0)
//	{
//		if ((thit[1] > 0.0) && (thit[1] <= dtfly))
//		{
//			double yhit2 = part_y + part_cy * thit[1];
//			if ((yhit2 >= ybody1) && (yhit2 <= ybody2))
//			{
//				if (thit[1] < thit_min)
//				{
//					thit_min = thit[1];
//					iwcol = 2.0;
//					double twhit = thit[1];
//					double xwhit = xbody2;
//					double ywhit = yhit2;
//					//壁では拡散反射モデルを採用
//					//パラメタ
//					double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double bb = 2.0 * M_PI*uniform_random();
//
//					double cxhit = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double cyhit = aa*cos(bb);
//					double czhit = aa*sin(bb);
//
//					//返り値の用意
//					part_hit[0] = iwcol;
//					part_hit[1] = twhit;
//					part_hit[2] = xwhit;
//					part_hit[3] = ywhit;
//					part_hit[4] = cxhit;
//					part_hit[5] = cyhit;
//					part_hit[6] = czhit;
//
//					
//				}
//			}
//		}
//	}
//
//											//ybodyで当たるか判定する
//	if (part_cy != 0.0)
//	{
//		//ybody1で当たるか判定
//		if ((thit[2] > 0.0) && (thit[2] <= dtfly))
//		{
//			double xhit1 = part_x + part_cx * thit[2];
//			if ((xhit1 >= xbody1) && (xhit1 <= xbody2))
//			{
//				if (thit[2] < thit_min)
//				{
//					thit_min = thit[2];
//					iwcol = 3.0;
//					double twhit = thit[2];
//					double xwhit = xhit1;
//					double ywhit = ybody1;
//					//壁では拡散反射モデルを採用
//					//パラメタ
//					double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double bb = 2.0 * M_PI*uniform_random();
//
//					double cyhit = -sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double cxhit = aa*cos(bb);
//					double czhit = aa*sin(bb);
//
//					//返り値の用意
//					part_hit[0] = iwcol;
//					part_hit[1] = twhit;
//					part_hit[2] = xwhit;
//					part_hit[3] = ywhit;
//					part_hit[4] = cxhit;
//					part_hit[5] = cyhit;
//					part_hit[6] = czhit;
//
//					
//				}
//			}
//		}
//	}
//	if (part_cy != 0.0)
//	{
//		//ybody2で当たるか判定
//		if ((thit[3] > 0.0) && (thit[3] <= dtfly))
//		{
//			double xhit2 = part_x + part_cx * thit[3];
//			if ((xhit2 >= xbody1) && (xhit2 <= xbody2))
//			{
//				if (thit[3] < thit_min)
//				{
//					thit_min = thit[3];
//					iwcol = 4.0;
//					double twhit = thit[3];
//					double xwhit = xhit2;
//					double ywhit = ybody2;
//					//壁では拡散反射モデルを採用
//					//パラメタ
//					double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double bb = 2.0 * M_PI*uniform_random();
//
//					double cyhit = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
//					double cxhit = aa*cos(bb);
//					double czhit = aa*sin(bb);
//
//					//返り値の用意
//					part_hit[0] = iwcol;
//					part_hit[1] = twhit;
//					part_hit[2] = xwhit;
//					part_hit[3] = ywhit;
//					part_hit[4] = cxhit;
//					part_hit[5] = cyhit;
//					part_hit[6] = czhit;
//					
//				}
//
//			}
//
//		}
//	}
//
//	return 0;
//}

int particle_injection(double xbuf1, double xbuf2, double ybuf1, double ybuf2, double part_xy[][2], double part_c[][3])
{
	//injection from buffer 1
	{
		double xzone1 = xbuf1;
		double xzone2 = xmin;
		double yzone1 = ybuf1;
		double yzone2 = ybuf2;
		injection_update(xzone1, xzone2, yzone1, yzone2, part_xy, part_c);
	}

	//injection from buffer 2
	{
		double xzone1 = xmax;
		double xzone2 = xbuf2;
		double yzone1 = ybuf1;
		double yzone2 = ybuf2;
		injection_update(xzone1, xzone2, yzone1, yzone2, part_xy, part_c);
	}

	//injection from buffer 3
	{
		double xzone1 = xmin;
		double xzone2 = xmax;
		double yzone1 = ybuf1;
		double yzone2 = ymin;
		injection_update(xzone1, xzone2, yzone1, yzone2, part_xy, part_c);
	}

	//injection from buffer 4
	{
		double xzone1 = xmin;
		double xzone2 = xmax;
		double yzone1 = ymax;
		double yzone2 = ybuf2;
		injection_update(xzone1, xzone2, yzone1, yzone2, part_xy, part_c);
	}
	return 0;
}

int re_indexing(double part_xy[][2], double part_c[][3], int lc[][my], int lcr[nmax], int lc0[][my], int i_j_cel[][2])
{
	//粒子mがどこの領域(i,j)の中にあるか判定
	//その個数をlcに格納している

	/*printf("%d", nmol);*/

	for (int m = 1; m <= nmol; m++)
	{
		int i = int((part_xy[m][0] - xmin) / dx);
		int j = int((part_xy[m][1] - ymin) / dy);

		//if (m % 10 == 0)
		//{
		//	printf("%d\n", m);
		//}
		//printf("i = %d, j = %d\n", i, j);

		i_j_cel[m][0] = i;
		i_j_cel[m][1] = j;
		lc[i][j] = lc[i][j] + 1;
	}

	//(0,0)から(i,j-1)までの個数をlc0に格納
	//これは次の処理のための布石
	int nsum = 0;
	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			lc0[i][j] = nsum;
			nsum = nsum + lc[i][j];
			lc[i][j] = 0;
		}
	}

	//(0,0)から(i,j-1)までの箱の中の粒子の個数　+ (i,j)の箱の中ですでに数えた個数　を k　としている
	//それを m と結び付けている
	for (int m = 1; m <= nmol; m++)
	{
		int i = i_j_cel[m][0];
		int j = i_j_cel[m][1];

		//if (m % 10 == 0)
		//{
		//	printf("%d\n", m);
		//	printf("i = %d, j = %d\n", i, j);
		//}

		lc[i][j] = lc[i][j] + 1;
		int k = lc0[i][j] + lc[i][j];
		//printf("k = %d, lc = %d, lc0 = %d\n", k, lc[i][j], lc0[i][j]);
		lcr[k] = m;				//<------おかしなアクセス
	}

	return 0;
}

int collision(double part_xy[][2], double part_c[][3], double col[mx][my], int lc[][my], int lcr[nmax], int lc0[][my])
{
	for (int i = 0; i < mx; i++)
	{
		for (int j = 0; j < my; j++)
		{
			//最大衝突数の設定
			double acol = col[i][j] * double(lc[i][j] * (lc[i][j] - 1));
			int ncol = int(acol);
			if (uniform_random() < (acol - double(ncol)))
			{
				ncol = ncol + 1;
			}

			//ncol<=0　や　lc[i][j]<=1　の場合は除く
			if ((ncol > 0) && (lc[i][j] > 1))
			{
				for (int nc = 0; nc < ncol; nc++)
				{
					//mは粒子番号l,nを選ぶためのパラメタ
					int m = int(uniform_random()*double(lc[i][j])) + lc0[i][j] + 1;
					int l = lcr[m];

					m = int(uniform_random()*double(lc[i][j])) + lc0[i][j] + 1;
					int n = lcr[m];

					//l=nとなったら粒子番号を再設定
					while (l == n)
					{
						m = int(uniform_random()*double(lc[i][j])) + lc0[i][j] + 1;
						n = lcr[m];

						//printf("lc = %d, lc0 = %d", lc[i][j], lc0[i][j]);
						//printf("m = %d, n = %d\n", m, n);

					}
					//							<-- (l,n) is collision pair.

					//gm は相対速度の大きさ
					double gm = sqrt(pow((part_c[l][0] - part_c[n][0]), 2.0) + pow((part_c[l][1] - part_c[n][1]), 2.0) + pow((part_c[l][2] - part_c[n][2]), 2.0));

					//最大衝突法
					if (uniform_random() <= (gm / gmax))
					{
						double cs = 1.0 - 2.0*uniform_random();
						double sn = sqrt(1.0 - pow(cs, 2.0));
						double bb = 2.0*M_PI*uniform_random();
						double rx = sn*cos(bb);
						double ry = sn*sin(bb);
						double rz = cs;

						double cxt = part_c[l][0] + part_c[n][0];
						double cyt = part_c[l][1] + part_c[n][1];
						double czt = part_c[l][2] + part_c[n][2];

						//粒子が衝突した後の速度に更新
						part_c[l][0] = 0.5*(cxt - gm*rx);
						part_c[l][1] = 0.5*(cyt - gm*ry);
						part_c[l][2] = 0.5*(czt - gm*rz);
						part_c[n][0] = 0.5*(cxt + gm*rx);
						part_c[n][1] = 0.5*(cyt + gm*ry);
						part_c[n][2] = 0.5*(czt + gm*rz);

					}
				}
			}
		}
	}
	return 0;
}

//一様分布を返す関数
//ランダム値を毎回出すためには毎回変数に代入すること
double uniform_random()
{
	static std::random_device rd;
	static std::mt19937 mt(rd());
	static std::uniform_real_distribution<double> dist(0.0, 1.0);

	return dist(mt);
}

//thitの最小値を求めるだけのゴミ関数
double min_thit(double thit[4])
{
	double min = 100000.0;

	for (int i = 0; i < 4; i++)
	{
		if ((thit[i] < min) && (thit[i] > 0))
		{
			min = thit[i];
		}
	}

	return min;
}

double injection_update(double xzone1, double xzone2, double yzone1, double yzone2, double part_xy[][2], double part_c[][3])
{
	double ainjec = (xzone2 - xzone1)*(yzone2 - yzone1)*double(ns);
	double ninjec = int(ainjec);
	if (uniform_random() < (ainjec - double(ninjec)))
	{
		ninjec = ninjec + 1;
	}

	for (int nn = 0; nn < ninjec; nn++)
	{
		double xst = xzone1 + (xzone2 - xzone1)*double(uniform_random());
		double yst = yzone1 + (yzone2 - yzone1)*double(uniform_random());

		double aa = sqrt(-2.0*rgas*t0*log(uniform_random()));
		double bb = 2 * M_PI*log(uniform_random());
		double cxt = (aa*cos(bb) + v0) / vref;
		double cyt = aa*sin(bb) / vref;

		aa = sqrt(-2.0*rgas*t0*log(uniform_random()));
		bb = 2 * M_PI*log(uniform_random());
		double czt = aa*cos(bb) / vref;

		double xin = xst + dt*cxt;
		double yin = yst + dt*cyt;

		//計算領域内に入るか判定するパラメタ
		double crtitx = (xin - xmin)*(xin - xmax);
		double crtity = (yin - ymin)*(yin - ymax);
		if ((crtitx < 0.0) && ((crtity < 0.0)))
		{
			//各種値をnmolの後ろに付けたし
			nmol = nmol + 1;
			part_c[nmol][0] = cxt;
			part_c[nmol][1] = cyt;
			part_c[nmol][2] = czt;
			part_xy[nmol][0] = xin;
			part_xy[nmol][1] = yin;
		}
		//printf("cxt = %f, cyt = %f, czt = %f\n", cxt, cyt, czt);
	}
	return 0.0;
}

double body_force_and_heat(double part_hit[], double c_x, double c_y, double xbody_force_x[][B_my], double xbody_force_y[][B_my], double xbody_energy[][B_my], double ybody_force_x[][B_mx], double ybody_force_y[][B_mx], double ybody_energy[][B_mx])
{
	int hit_num = int(part_hit[0]);
	double s = part_hit[7];
	double xhit = part_hit[2];
	double yhit = part_hit[3];

	double momentum_gap_x = mass_part_2*(part_hit[4] - c_x);
	double momentum_gap_y = mass_part_2*(part_hit[5] - c_y);
	double energy_gap = 0.5*mass_part_2*((pow(c_x, 2.0) + pow(c_y, 2.0)) - (pow(part_hit[4], 2.0) + pow(part_hit[5], 2.0)));

	//xbodyに衝突したとき
	if (hit_num == 1)
	{
		int i = int(s*B_my);
		xbody_force_x[hit_num - 1][i] += momentum_gap_x;
		xbody_force_y[hit_num - 1][i] += momentum_gap_y;
		xbody_energy[hit_num - 1][i] += energy_gap;
	}
	if (hit_num == 3)
	{
		int i = int(s*B_my);
		xbody_force_x[hit_num - 2][i] += momentum_gap_x;
		xbody_force_y[hit_num - 2][i] += momentum_gap_y;
		xbody_energy[hit_num - 2][i] += energy_gap;
	}

	//ybodyに衝突したとき
	if (hit_num == 2)
	{
		int i = int(s*B_mx);
		ybody_force_x[hit_num - 1][i] += momentum_gap_x;		// -3としているのは配列の要素数が2しかないのに揃えるため
		ybody_force_y[hit_num - 1][i] += momentum_gap_y;
		ybody_energy[hit_num - 1][i] += energy_gap;
	}
	if (hit_num == 4)
	{
		int i = int(s*B_mx);
		ybody_force_x[hit_num - 4][i] += momentum_gap_x;		// -3としているのは配列の要素数が2しかないのに揃えるため
		ybody_force_y[hit_num - 4][i] += momentum_gap_y;
		ybody_energy[hit_num - 4][i] += energy_gap;
	}

	return 0;
}

double line_body_collsion(int mol, double body_point_pair[][2][2], double dtfly, double part_x, double part_y, double part_cx, double part_cy, double part_hit[])
{
	//sはbody上の位置かを調べるパラメタ，０～１の範囲
	//tは粒子の１ステップ上の軌跡を表すパラメタ，０～１の範囲
	double s, t;
	double t_min = 1000.0;			//適当に大きな数．意味はない．
	for (int n = 0; n < sides; n++)
	{
		double a = body_point_pair[n][0][0];
		double b = body_point_pair[n][0][1];
		double c = body_point_pair[n][1][0];
		double d = body_point_pair[n][1][1];

		s = ((part_x - a)*part_cy - (part_y - b)*part_cx) / ((c - a)*part_cy - (d - b)*part_cx);
		t = ((part_x - a)*(d - b)*part_cx - (part_y - b)*(c - a)*part_cx)
			/ (part_cx*dt*((c - a)*part_cy - (d - b)*part_cx));
		if ((0 <= s) && (s <= 1) && (1.0e-9 < t) && (t <= 1) && (t < t_min))
		{
			//最小のtの時にのみ衝突
			t_min = t;
			double iwcol = double(n) + 1.0;
			double twhit = t*dt;
			double xwhit = (1 - s)*a + s*c;
			double ywhit = (1 - s)*b + s*d;
			//壁では拡散反射モデルを採用
			//パラメタ
			double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
			double bb = 2.0 * M_PI*uniform_random();

			//(1,0)方向に向かって放出するのが初期
			double cxhit = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
			double cyhit = aa*cos(bb);
			double czhit = aa*sin(bb);

			double c_xy_before[2] = { cxhit, cyhit };
			double theta_body = angle_vectors(c - a, d - b);
			double c_xy_after[2];

			//拡散反射する方向に傾ける
			rotation(theta_body + M_PI / 2, c_xy_before, c_xy_after);

			//速度更新
			cxhit = c_xy_after[0];
			cyhit = c_xy_after[1];

			//返り値の用意
			part_hit[0] = iwcol;
			part_hit[1] = twhit;
			part_hit[2] = xwhit;
			part_hit[3] = ywhit;
			part_hit[4] = cxhit;
			part_hit[5] = cyhit;
			part_hit[6] = czhit;
			part_hit[7] = s;

			//printf("iwcol = %f\n", iwcol);
			//printf("%14.15f\n", twhit);
			//printf("%f\n", xwhit);
			//printf("%f\n", ywhit);
			//printf("%f\n", cxhit);
			//printf("%f\n", cyhit);
			//printf("%f\n", czhit);

		}
	}
	return 0;
}

double rotation(double theta, double xy[], double after_xy[])
{
	//theta はラジアンで
	after_xy[0] = cos(theta)*xy[0] - sin(theta)*xy[1];
	after_xy[1] = sin(theta)*xy[0] + cos(theta)*xy[1];

	return 0;
}

double angle_vectors(double vect_x, double vect_y)
{
	double abs_vect = pow(vect_x, 2.0) + pow(vect_y, 2.0);
	double cos_theta = vect_x / pow(abs_vect, 0.5);
	if (vect_y > 0)
	{
		return acos(cos_theta);
	}
	if (vect_y < 0) {
		return -acos(cos_theta);
	}
	else
	{
		return acos(cos_theta);
	}
}

double body_rotate(double theta)
{
	double buffer[2];					//座標変換のために一時的に使う変数置き場
										//使う直線の数は指定する必要あり
	for (int i = 0; i < sides; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			rotation(theta, body_point[i][j], buffer);
			body_point[i][j][0] = buffer[0];
			body_point[i][j][1] = buffer[1];
		}
	}
	return 0;
}

double delete_particle_in_body()
{
	for (int m = 1; m <= nmol; m++)
	{
		for (int i = 0; i <= sides; i++)
		{
			double p1x = body_point[i][0][0];
			double p1y = body_point[i][0][1];
			double p2x = body_point[i][1][0];
			double p2y = body_point[i][1][1];

			double x = particle_x_y[m][0];
			double y = particle_x_y[m][1];

			double a_cross = (p1x - p2x)*(y - p1y) - (x - p1x)*(p1y - p2y);
			double b_cross = (-p1x)*y - x*(-p1y);
			double c_cross = p2x*(y - p2y) - (x - p2x)*p2y;

			//三角形の中に含まれる時，粒子を削除
			if (((a_cross >= 0) && (b_cross >= 0) && (c_cross >= 0)) ||
				((a_cross <= 0) && (b_cross <= 0) && (c_cross <= 0)))
			{
				particle_x_y[m][0] = particle_x_y[nmol][0];
				particle_x_y[m][1] = particle_x_y[nmol][1];
				particle_c[m][0] = particle_c[nmol][0];
				particle_c[m][1] = particle_c[nmol][1];
				particle_c[m][2] = particle_c[nmol][2];
				nmol = nmol - 1;
			}
		}
	}
	return 0;
}

double calc_torque(double hit_x, double hit_y, double u, double v)
{
	double part_momentum;			//一つの粒子の物体に及ぼす力積
	double alpha;
	double beta;
	double gamma;

	alpha = angle_vectors(hit_x, hit_y);
	beta = angle_vectors(u, v);
	gamma = beta - alpha;

	part_momentum = mass_particle*pow(u*u + v*v, 0.5)*sin(gamma)*pow(hit_x*hit_x + hit_y*hit_y, 0.5);

	return part_momentum;
}

double intersec_lines(double a_bef[2], double b_bef[2], double a_aft[2], double b_aft[2], double c[2])
{
	//a_bef と b_bef を結ぶ直線と，a_aft と b_bef　を結ぶ直線の交点を求める関数
	double a1 = a_bef[0];
	double a2 = a_bef[1];
	double b1 = b_bef[0];
	double b2 = b_bef[1];
	double aa1 = a_aft[0];
	double aa2 = a_aft[1];
	double bb1 = b_aft[0];
	double bb2 = b_aft[1];

	c[0] = ((a1*b2 - a2*b1)*(bb1 - aa1) - (aa1*bb2 - aa2*bb1)*(b1 - a1)) / ((bb1 - aa1)*(b2 - a2) - (b1 - a1)*(bb2 - aa2));
	c[1] = ((a2*b1 - a1*b2)*(bb2 - aa2) - (aa2*bb1 - aa1*bb2)*(b2 - a2)) / ((bb2 - aa2)*(b1 - a1) - (b2 - a2)*(bb1 - aa1));
	//printf("c_frac, %f,%f\n", ((bb1 - aa1)*(b2 - a2) - (b1 - a1)*(bb2 - aa2)), ((bb2 - aa2)*(b1 - a1) - (b2 - a2)*(bb1 - aa1)));

	return 0;
}

double near_point_on_line(double a_bef[2], double b_bef[2], double d[2])
{
	double a1 = a_bef[0];
	double a2 = a_bef[1];
	double b1 = b_bef[0];
	double b2 = b_bef[1];

	d[0] = ((a1*b2 - a2*b1)*(b2 - a2)) / ((b2 - a2)*(b2 - a2) + (b1 - a1)*(b1 - a1));
	d[1] = -((a1*b2 - a2*b1)*(b1 - a1)) / ((b2 - a2)*(b2 - a2) + (b1 - a1)*(b1 - a1));

	return 0;
}

double rotate_collsion_line(double dtheta, double sum_torque)
{
	double torque = sum_torque;

	if ( (dtheta > -1.0e-8) && (dtheta < 1.0e-8 ) )
	{
		return 0;
	}
	else {
		for (int mol = 1; mol <= nmol; mol++)
		{
			//最初にどの辺に当たるか判定．回転角のcosをとり，大きい方が回転角はちいさい．
			double cos_dphi_min = 1;

			//新しい粒子の座標を格納する場所
			double x_part_new;
			double y_part_new;
			double v_part_new[2] = {};

			//本当に衝突したときの値を格納
			double xy_part_true[2];
			double v_part_true[2];

			//flag はx_part_new などが更新されたときに1となり，更新するルーチンを動かす．
			int hit_flag = 0;
			double partxy[2] = { particle_x_y[mol][0],particle_x_y[mol][1] };
			double r_part_sq = partxy[0] * partxy[0] + partxy[1] * partxy[1];

			for (int i = 0; i < sides; i++)
			{
				//回転前の座標
				//Aの方がBより原点から遠い方になるようにする．
				double dummy_r_a_sq = body_point[i][0][0] * body_point[i][0][0] + body_point[i][0][1] * body_point[i][0][1];
				double dummy_r_b_sq = body_point[i][1][0] * body_point[i][1][0] + body_point[i][1][1] * body_point[i][1][1];
				double a1[2] = { body_point[i][0][0],body_point[i][0][1] };
				double b1[2] = { body_point[i][1][0],body_point[i][1][1] };

				//A点の方がB点の方より遠いように指定する
				if (dummy_r_a_sq < dummy_r_b_sq)
				{
					a1[0] = body_point[i][1][0];
					a1[1] = body_point[i][1][1];
					b1[0] = body_point[i][0][0];
					b1[1] = body_point[i][0][1];
				}
				double r_a_sq = a1[0] * a1[0] + a1[1] * a1[1];
				double r_b_sq = b1[0] * b1[0] + b1[1] * b1[1];

				//回転後の座標
				double a2[2];
				double b2[2];
				rotation(dtheta, a1, a2);
				rotation(dtheta, b1, b2);

				//交点の座標
				double c[2];
				intersec_lines(a1, b1, a2, b2, c);
				//printf("dtheta:%e\n", dtheta);
				//printf("c: %f,%f\n", c[0], c[1]);

				//直線に垂線を下した足
				double d1[2];
				double d2[2];
				near_point_on_line(a1, b1, d1);
				near_point_on_line(a2, b2, d2);
				double r_d_sq = d1[0] * d1[0] + d1[1] * d1[1];
				double d2_theta = angle_vectors(d2[0], d2[1]);

				//直線の間はどこか判定するためのサンプル
				double e1[2];
				double e2[2];
				rotation(dtheta / 2, a1, e1);
				rotation(dtheta / 2, b1, e2);

				//c,dがABの間にあるか判定するパラメタ
				double c_par_1 = (c[0] - a1[0])*(c[0] - b1[0]);
				double c_par_2 = (c[1] - a1[1])*(c[1] - b1[1]);
				double d_par_1 = (d1[0] - a1[0])*(d1[0] - b1[0]);
				double d_par_2 = (d1[1] - a1[1])*(d1[1] - b1[1]);

				//オーバーフローが起きるほど大きな値は入れないと思って，積で判定
				//eとparticleが同じ側にあるか調べる
				//AとA'の間にあるかのパラメタ
				double e1_par_1 = (b1[0] - a1[0])*e1[1] - (b1[1] - a1[1])*e1[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double e1_par_2 = (b2[0] - a2[0])*e1[1] - (b2[1] - a2[1])*e1[0] + a2[0] * b2[1] - a2[1] * b2[0];

				//BとB'の間にあるかのパラメタ
				double e2_par_1 = (b1[0] - a1[0])*e2[1] - (b1[1] - a1[1])*e2[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double e2_par_2 = (b2[0] - a2[0])*e2[1] - (b2[1] - a2[1])*e2[0] + a2[0] * b2[1] - a2[1] * b2[0];

				//原点側にあるかどうかの判定
				double o_par_1 = a1[0] * b1[1] - a1[1] * b1[0];
				double o_par_2 = a2[0] * b2[1] - a2[1] * b2[0];

				//particleが同じ側にあるかのパラメタ
				double p1_par_1 = (b1[0] - a1[0])*partxy[1] - (b1[1] - a1[1])*partxy[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double p1_par_2 = (b2[0] - a2[0])*partxy[1] - (b2[1] - a2[1])*partxy[0] + a2[0] * b2[1] - a2[1] * b2[0];

				//dがABの間にあるとき
				if ((d_par_1 <= 0) && (d_par_2 <= 0))
				{
					if (dtheta < 0)
					{
						//(a)の時
						//if (partxy[0] <= c[0])
						//{
						if ((e1_par_1*p1_par_1 > 0) && (e1_par_2*p1_par_2 > 0) && (r_part_sq < r_a_sq))
						{
							double s2 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_a_sq - r_d_sq, 0.5);
							c_new_generate(a2, b2, d2_theta, v_part_new);
							x_part_new = (1 - s2)*d2[0] + s2*a2[0] + v_part_new[0] * 1.0e-8;
							y_part_new = (1 - s2)*d2[1] + s2*a2[1] + v_part_new[1] * 1.0e-8;
							hit_flag = 1;
							//もしこの方向で進んだ時に，まだ線分sweep内にいるなら逆方向に変換する．
							double sample_par_1 = (b1[0] - a1[0]) * y_part_new - (b1[1] - a1[1]) * x_part_new + a1[0] * b1[1] - a1[1] * b1[0];
							double sample_par_2 = (b2[0] - a2[0]) * y_part_new - (b2[1] - a2[1]) * x_part_new + a2[0] * b2[1] - a2[1] * b2[0];
							double sample_r_sq = x_part_new*x_part_new + y_part_new*y_part_new;
							if ((e1_par_1*sample_par_1 > 0) && (e1_par_2*sample_par_2 > 0) && (sample_r_sq < r_a_sq))
							{
								c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
								x_part_new = (1 - s2)*d2[0] + s2*a2[0] + v_part_new[0] * 1.0e-8;
								y_part_new = (1 - s2)*d2[1] + s2*a2[1] + v_part_new[1] * 1.0e-8;
							}
						}
						//}

						//(b)の時
						//if (partxy[0] > c[0])
						//{
						if ((e2_par_1*p1_par_1 > 0) && (e2_par_2*p1_par_2 > 0) && (r_part_sq < r_b_sq))
						{
							double s1 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_b_sq - r_d_sq, 0.5);
							//c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
							c_new_generate(a2, b2, d2_theta , v_part_new);
							x_part_new = (1 - s1)*d2[0] + s1*b2[0] + v_part_new[0] * 1.0e-8;
							y_part_new = (1 - s1)*d2[1] + s1*b2[1] + v_part_new[1] * 1.0e-8;
							hit_flag = 1;
							//もしこの方向で進んだ時に，まだ線分sweep内にいるなら逆方向に変換する．
							double sample_par_1 = (b1[0] - a1[0]) * y_part_new - (b1[1] - a1[1]) * x_part_new + a1[0] * b1[1] - a1[1] * b1[0];
							double sample_par_2 = (b2[0] - a2[0]) * y_part_new - (b2[1] - a2[1]) * x_part_new + a2[0] * b2[1] - a2[1] * b2[0];
							double sample_r_sq = x_part_new*x_part_new + y_part_new*y_part_new;
							if ((e1_par_1*sample_par_1 > 0) && (e1_par_2*sample_par_2 > 0) && (sample_r_sq < r_b_sq))
							{
								c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
								x_part_new = (1 - s1)*d2[0] + s1*b2[0] + v_part_new[0] * 1.0e-8;
								y_part_new = (1 - s1)*d2[1] + s1*b2[1] + v_part_new[1] * 1.0e-8;
							}
						}
						//}

					}
					if (dtheta > 0)
					{

						//(a')の時
						//if (partxy[0] <= c[0])
						//{
						if ((e1_par_1*p1_par_1 > 0) && (e1_par_2*p1_par_2 > 0) && (r_part_sq < r_a_sq))
						{
							double s2 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_a_sq - r_d_sq, 0.5);
							//c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
							c_new_generate(a2, b2, d2_theta , v_part_new);
							x_part_new = (1 - s2)*d2[0] + s2*a2[0] + v_part_new[0] * 1.0e-8;
							y_part_new = (1 - s2)*d2[1] + s2*a2[1] + v_part_new[1] * 1.0e-8;
							hit_flag = 1;
							//もしこの方向で進んだ時に，まだ線分sweep内にいるなら逆方向に変換する．
							double sample_par_1 = (b1[0] - a1[0]) * y_part_new - (b1[1] - a1[1]) * x_part_new + a1[0] * b1[1] - a1[1] * b1[0];
							double sample_par_2 = (b2[0] - a2[0]) * y_part_new - (b2[1] - a2[1]) * x_part_new + a2[0] * b2[1] - a2[1] * b2[0];
							double sample_r_sq = x_part_new*x_part_new + y_part_new*y_part_new;
							if ((e1_par_1*sample_par_1 > 0) && (e1_par_2*sample_par_2 > 0) && (sample_r_sq < r_a_sq))
							{
								c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
								x_part_new = (1 - s2)*d2[0] + s2*a2[0] + v_part_new[0] * 1.0e-8;
								y_part_new = (1 - s2)*d2[1] + s2*a2[1] + v_part_new[1] * 1.0e-8;
							}
						}
						//}

						//(b')の時
						//if (partxy[0] > c[0])
						//{
						if ((e2_par_1*p1_par_1 > 0) && (e2_par_2*p1_par_2 > 0) && (r_part_sq < r_b_sq))
						{
							double s1 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_b_sq - r_d_sq, 0.5);
							c_new_generate(a2, b2, d2_theta, v_part_new);
							x_part_new = (1 - s1)*d2[0] + s1*b2[0] +v_part_new[0] * 1.0e-8;
							y_part_new = (1 - s1)*d2[1] + s1*b2[1] +v_part_new[1] * 1.0e-8;
							hit_flag = 1;
							//もしこの方向で進んだ時に，まだ線分sweep内にいるなら逆方向に変換する．
							double sample_par_1 = (b1[0] - a1[0]) * y_part_new - (b1[1] - a1[1]) * x_part_new + a1[0] * b1[1] - a1[1] * b1[0];
							double sample_par_2 = (b2[0] - a2[0]) * y_part_new - (b2[1] - a2[1]) * x_part_new + a2[0] * b2[1] - a2[1] * b2[0];
							double sample_r_sq = x_part_new*x_part_new + y_part_new*y_part_new;
							if ((e1_par_1*sample_par_1 > 0) && (e1_par_2*sample_par_2 > 0) && (sample_r_sq < r_b_sq))
							{
								c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
								x_part_new = (1 - s1)*d2[0] + s1*b2[0] + v_part_new[0] * 1.0e-8;
								y_part_new = (1 - s1)*d2[1] + s1*b2[1] + v_part_new[1] * 1.0e-8;
							}
						}
						//}

					}

					//角度がOD,OD'の間に入っているか確かめる
					double angle_part = angle_vectors(partxy[0], partxy[1]);
					double angle_d1 = angle_vectors(d1[0], d1[1]);

					//-180と180の間の判定のためにパラメタを複数用意
					double theta_par_1 = (angle_part - angle_d1)*(angle_part - angle_d1 - dtheta);
					double theta_par_2 = (angle_part - angle_d1 - 2 * M_PI)*(angle_part - angle_d1 - dtheta - 2 * M_PI);
					double theta_par_3 = (angle_part - angle_d1 + 2 * M_PI)*(angle_part - angle_d1 - dtheta + 2 * M_PI);

					if ((theta_par_1 < 0) || (theta_par_2 < 0) || (theta_par_3 < 0))
					{
						if ((o_par_1*p1_par_1 > 0) && (o_par_2*p1_par_2 > 0) && (r_part_sq > r_d_sq))
						{
							if (dtheta < 0)
							{
								//cもABうちなら無条件にBD上につく
								if ((c_par_1 <= 0) && (c_par_2 <= 0))
								{
									double s1 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_b_sq - r_d_sq, 0.5);
									x_part_new = (1 - s1)*d2[0] + s1*b2[0];
									y_part_new = (1 - s1)*d2[1] + s1*b2[1];

									c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
									hit_flag = 1;
								}
								//cがAB上に無いなら，距離によって当たる場所が変化する
								else
								{
									if (r_part_sq >= r_b_sq)
									{
										double s2 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_a_sq - r_d_sq, 0.5);
										x_part_new = (1 - s2)*d2[0] + s2*a2[0];
										y_part_new = (1 - s2)*d2[1] + s2*a2[1];

										c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
										hit_flag = 1;
									}
									else
									{
										double s1 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_b_sq - r_d_sq, 0.5);
										x_part_new = (1 - s1)*d2[0] + s1*b2[0];
										y_part_new = (1 - s1)*d2[1] + s1*b2[1];

										c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
										hit_flag = 1;
									}
								}
							}
							if (dtheta > 0)
							{
								//cがAB上にあれば無条件に
								if ((c_par_1 <= 0) && (c_par_2 <= 0))
								{
									double s2 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_a_sq - r_d_sq, 0.5);
									x_part_new = (1 - s2)*d2[0] + s2*a2[0];
									y_part_new = (1 - s2)*d2[1] + s2*a2[1];

									c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
									hit_flag = 1;
								}
								//cがAB上に無いなら，距離によって当たる場所が変化する
								else
								{
									if (r_part_sq >= r_a_sq)
									{
										double s1 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_b_sq - r_d_sq, 0.5);
										x_part_new = (1 - s1)*d2[0] + s1*b2[0];
										y_part_new = (1 - s1)*d2[1] + s1*b2[1];

										c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
										hit_flag = 1;
									}
									else
									{
										double s2 = pow(r_part_sq - r_d_sq, 0.5) / pow(r_a_sq - r_d_sq, 0.5);
										x_part_new = (1 - s2)*d2[0] + s2*a2[0];
										y_part_new = (1 - s2)*d2[1] + s2*a2[1];

										c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
										hit_flag = 1;
									}
								}
							}
						}
					}
				}

				//c,dがともにABの外にあるとき
				else
				{
					//(a)の時
					if ((e1_par_1*p1_par_1 > 0) && (e1_par_2*p1_par_2 > 0) && (r_part_sq > r_b_sq) && (r_part_sq < r_a_sq))
					{
						double s3 = (pow(r_part_sq - r_d_sq, 0.5) - pow(r_b_sq - r_d_sq, 0.5)) / (pow(r_a_sq - r_d_sq, 0.5) - pow(r_b_sq - r_d_sq, 0.5));
						c_new_generate(a2, b2, d2_theta, v_part_new);
						x_part_new = s3*a2[0] + (1 - s3)*b2[0] + v_part_new[0] * 1.0e-8;
						y_part_new = s3*a2[1] + (1 - s3)*b2[1] + v_part_new[1] * 1.0e-8;
						//もしこの方向で進んだ時に，まだ線分sweep内にいるなら逆方向に変換する．
						double sample_par_1 = (b1[0] - a1[0]) * y_part_new - (b1[1] - a1[1]) * x_part_new + a1[0] * b1[1] - a1[1] * b1[0];
						double sample_par_2 = (b2[0] - a2[0]) * y_part_new - (b2[1] - a2[1]) * x_part_new + a2[0] * b2[1] - a2[1] * b2[0];
						double sample_r_sq = x_part_new*x_part_new + y_part_new*y_part_new;
						if ((e1_par_1*sample_par_1 > 0) && (e1_par_2*sample_par_2 > 0) && (sample_r_sq > r_b_sq) && (sample_r_sq < r_a_sq))
						{
							c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
							x_part_new = (1 - s3)*a2[0] + s3*b2[0] + v_part_new[0] * 1.0e-8;
							y_part_new = (1 - s3)*a2[1] + s3*b2[1] + v_part_new[1] * 1.0e-8;
						}
						//if (dtheta < 0)
						//{
						//	c_new_generate(a2, b2, d2_theta, v_part_new);
						//	hit_flag = 1;
						//}
						//if (dtheta > 0)
						//{
						//	c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
						//	hit_flag = 1;
						//}
					}

					double angle_part = angle_vectors(partxy[0], partxy[1]);
					double angle_b1 = angle_vectors(b1[0], b1[1]);

					//-180と180の間の判定のためにパラメタを複数用意
					double theta_par_1 = (angle_part - angle_b1)*(angle_part - angle_b1 - dtheta);
					double theta_par_2 = (angle_part - angle_b1 - 2 * M_PI)*(angle_part - angle_b1 - dtheta - 2 * M_PI);
					double theta_par_3 = (angle_part - angle_b1 + 2 * M_PI)*(angle_part - angle_b1 - dtheta + 2 * M_PI);

					if ((theta_par_1 < 0) || (theta_par_2 < 0) || (theta_par_3 < 0))
					{
						if ((o_par_1*p1_par_1 > 0) && (o_par_2*p1_par_2 > 0) && (r_part_sq > r_b_sq))
						{
							double s3 = (pow(r_part_sq - r_d_sq, 0.5) - pow(r_b_sq - r_d_sq, 0.5)) / (pow(r_a_sq - r_d_sq, 0.5) - pow(r_b_sq - r_d_sq, 0.5));
							x_part_new = s3*a2[0] + (1 - s3)*a2[0];
							y_part_new = s3*a2[1] + (1 - s3)*a2[1];

							if (dtheta < 0)
							{
								c_new_generate(a2, b2, d2_theta, v_part_new);
								hit_flag = 1;
							}
							if (dtheta > 0)
							{
								c_new_generate(a2, b2, d2_theta - M_PI, v_part_new);
								hit_flag = 1;
							}
						}
					}
				}

				//hit_flagが立っていたら更新
				if (hit_flag == 1)
				{
					//最短時間で衝突したときが正しい
					double cos_temp = (partxy[0] * x_part_new + partxy[1] * y_part_new) / (r_part_sq);
					if (cos_temp < cos_dphi_min)
					{
						//cosの最小値を更新
						cos_dphi_min = cos_temp;

						double torque1 = 0;
						double torque2 = 0;
						torque1 = calc_torque(particle_x_y[mol][0], particle_x_y[mol][1], particle_c[mol][0], particle_c[mol][1]);
						torque2 = -calc_torque(x_part_new, y_part_new, v_part_new[0], v_part_new[1]);

						torque += torque1 + torque2;
						//本当に衝突したときの値を格納
						xy_part_true[0] = x_part_new;
						xy_part_true[1] = y_part_new;
						v_part_true[0] = v_part_new[0];
						v_part_true[1] = v_part_new[1];
					}
				}
			}
			for (int i = 0; i < sides; i++)
			{

				//回転前の座標
				//Aの方がBより原点から遠い方になるようにする．
				double dummy_r_a_sq = body_point[i][0][0] * body_point[i][0][0] + body_point[i][0][1] * body_point[i][0][1];
				double dummy_r_b_sq = body_point[i][1][0] * body_point[i][1][0] + body_point[i][1][1] * body_point[i][1][1];
				double a1[2] = { body_point[i][0][0],body_point[i][0][1] };
				double b1[2] = { body_point[i][1][0],body_point[i][1][1] };
				if (dummy_r_a_sq < dummy_r_b_sq)
				{
					a1[0] = body_point[i][1][0];
					a1[1] = body_point[i][1][1];
					b1[0] = body_point[i][0][0];
					b1[1] = body_point[i][0][1];
				}
				double r_a_sq = a1[0] * a1[0] + a1[1] * a1[1];
				double r_b_sq = b1[0] * b1[0] + b1[1] * b1[1];

				//回転後の座標
				double a2[2];
				double b2[2];
				rotation(dtheta, a1, a2);
				rotation(dtheta, b1, b2);

				//交点の座標
				double c[2];
				intersec_lines(a1, b1, a2, b2, c);
				//printf("dtheta:%e\n", dtheta);
				//printf("c: %f,%f\n", c[0], c[1]);

				//直線に垂線を下した足
				double d1[2];
				double d2[2];
				near_point_on_line(a1, b1, d1);
				near_point_on_line(a2, b2, d2);
				double r_d_sq = d1[0] * d1[0] + d1[1] * d1[1];
				double d2_theta = angle_vectors(d2[0], d2[1]);

				//直線の間はどこか判定するためのサンプル
				double e1[2];
				double e2[2];
				rotation(dtheta / 2, a1, e1);
				rotation(dtheta / 2, b1, e2);

				//c,dがABの間にあるか判定するパラメタ
				double c_par_1 = (c[0] - a1[0])*(c[0] - b1[0]);
				double c_par_2 = (c[1] - a1[1])*(c[1] - b1[1]);
				double d_par_1 = (d1[0] - a1[0])*(d1[0] - b1[0]);
				double d_par_2 = (d1[1] - a1[1])*(d1[1] - b1[1]);

				//オーバーフローが起きるほど大きな値は入れないと思って，積で判定
				//eとparticleが同じ側にあるか調べる
				//AとA'の間にあるかのパラメタ
				double e1_par_1 = (b1[0] - a1[0])*e1[1] - (b1[1] - a1[1])*e1[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double e1_par_2 = (b2[0] - a2[0])*e1[1] - (b2[1] - a2[1])*e1[0] + a2[0] * b2[1] - a2[1] * b2[0];

				//BとB'の間にあるかのパラメタ
				double e2_par_1 = (b1[0] - a1[0])*e2[1] - (b1[1] - a1[1])*e2[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double e2_par_2 = (b2[0] - a2[0])*e2[1] - (b2[1] - a2[1])*e2[0] + a2[0] * b2[1] - a2[1] * b2[0];

				//原点側にあるかどうかの判定
				double o_par_1 = a1[0] * b1[1] - a1[1] * b1[0];
				double o_par_2 = a2[0] * b2[1] - a2[1] * b2[0];

				//particleが同じ側にあるかのパラメタ
				double p1_par_1 = (b1[0] - a1[0])*partxy[1] - (b1[1] - a1[1])*partxy[0] + a1[0] * b1[1] - a1[1] * b1[0];
				double p1_par_2 = (b2[0] - a2[0])*partxy[1] - (b2[1] - a2[1])*partxy[0] + a2[0] * b2[1] - a2[1] * b2[0];

				double p1x = a2[0];
				double p1y = a2[1];
				double p2x = b2[0];
				double p2y = b2[1];

				double x = xy_part_true[0];
				double y = xy_part_true[1];

				double a_cross = (p1x - p2x)*(y - p1y) - (x - p1x)*(p1y - p2y);
				double b_cross = (-p1x)*y - x*(-p1y);
				double c_cross = p2x*(y - p2y) - (x - p2x)*p2y;

				double p2_par_1 = (b1[0] - a1[0])*y - (b1[1] - a1[1])*x + a1[0] * b1[1] - a1[1] * b1[0];
				double p2_par_2 = (b2[0] - a2[0])*y - (b2[1] - a2[1])*x + a2[0] * b2[1] - a2[1] * b2[0];

				////三角形の中に含まれる時，粒子を削除
				//if (((a_cross > 0) && (b_cross > 0) && (c_cross > 0)) ||
				//	((a_cross < 0) && (b_cross < 0) && (c_cross < 0)))
				//{
				//	printf("hitflag : %d\n", hit_flag);
				//	printf("i:%d\n", i);
				//	printf("e1par: %f,%f\n", e1_par_1, e1_par_2);
				//	printf("e2par: %f,%f\n", e2_par_1, e2_par_2);
				//	printf("p_par : %f,%f\n", p1_par_1, p1_par_2);
				//	printf("p2_par : %e,%e\n", p2_par_1, p2_par_2);
				//	printf("parpar : %e\n", (e2_par_1*p1_par_1));
				//	printf("r_a_sq,r_b_sq,r_p : %f, %f, %f\n", r_a_sq, r_b_sq, r_part_sq);
				//	printf("c : %f,%f\n", c[0], c[1]);
				//	printf("bef_body : (%f,%f),(%f,%f)\n", a1[0], a1[1], b1[0], b1[1]);
				//	printf("p1(%f,%f), p2(%f,%f)\n", p1x, p1y, p2x, p2y);
				//	printf("part_xy_bef : %f,%f\n", particle_x_y[mol][0], particle_x_y[mol][1]);
				//	printf("part_xy_aft : %f,%f\n",xy_part_true[0], xy_part_true[1]);
				//	printf("\n");
				//}
			}
			if (hit_flag == 1)
			{
				double torque1 = 0;
				double torque2 = 0;
				torque1 = calc_torque(particle_x_y[mol][0], particle_x_y[mol][1], particle_c[mol][0], particle_c[mol][1]);
				torque2 = -calc_torque(x_part_new, y_part_new, v_part_new[0], v_part_new[1]);

				torque += torque1 + torque2;

				particle_x_y[mol][0] = xy_part_true[0];
				particle_x_y[mol][1] = xy_part_true[1];
				particle_c[mol][0] = v_part_true[0];
				particle_c[mol][1] = v_part_true[1];
			}
		}

		return torque;
	}
}

double c_new_generate(double a_dash[2], double b_dash[2], double theta_c, double v_part_new[2])
{
	//theta_cは線分に対して垂直な方向へのベクトルの角度

	double a = a_dash[0];
	double b = a_dash[1];
	double c = b_dash[0];
	double d = b_dash[1];

	//壁では拡散反射モデルを採用
	//パラメタ
	double aa = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
	double bb = 2.0 * M_PI * uniform_random();

	//(1,0)方向に向かって放出するのが初期
	double cxhit = sqrt(-2.0*rgas*twall*log(uniform_random())) / vref;
	double cyhit = aa*cos(bb);
	double czhit = aa*sin(bb);

	double c_xy_before[2] = { cxhit, cyhit };
	//double theta_body = angle_vectors(c - a, d - b);
	double c_xy_after[2];

	//拡散反射する方向に傾ける
	rotation(theta_c, c_xy_before, c_xy_after);

	v_part_new[0] = c_xy_after[0];
	v_part_new[1] = c_xy_after[1];

	return 0;
}

double check_particle_in_body()
{
	//長方形の中に粒子があるか，チェックする
	for (int m = 1; m <= nmol; m++)
	{
		for (int i = 0; i < sides; i++)
		{
			double p1x = body_point[i][0][0];
			double p1y = body_point[i][0][1];
			double p2x = body_point[i][1][0];
			double p2y = body_point[i][1][1];

			double x = particle_x_y[m][0];
			double y = particle_x_y[m][1];

			double a_cross = (p1x - p2x)*(y - p1y) - (x - p1x)*(p1y - p2y);
			double b_cross = (-p1x)*y - x*(-p1y);
			double c_cross = p2x*(y - p2y) - (x - p2x)*p2y;

			//三角形の中に含まれる時，粒子を削除
			if (((a_cross >= 0) && (b_cross >= 0) && (c_cross >= 0)) ||
				((a_cross <= 0) && (b_cross <= 0) && (c_cross <= 0)))
			{
				printf("i:%d\n", i);
				printf("p1(%f,%f), p2(%f,%f)\n",p1x, p1y, p2x, p2y);
				printf("part_xy : %f,%f\n", particle_x_y[m][0], particle_x_y[m][1]);
				printf("\n");
			}
		}
	}
	return 0;
}

//double simple_rotate_collision(double dtheta)
//{
//	if ((dtheta > -1.0e-8) && (dtheta < 1.0e-8))
//	{
//		return 0;
//	}
//	else {
//		for (int i = 0; i < sides; i++)
//		{
//			double a1[2] = { body_point[i][0][0],body_point[i][0][1] };
//			double b1[2] = { body_point[i][1][0],body_point[i][1][1] };
//		}
//	}
//
//}