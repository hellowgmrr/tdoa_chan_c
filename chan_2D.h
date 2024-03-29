#ifndef chan_2D_h
#define chan_2D_h
#include<matrixmul.h>
#include<math.h>
#include<invert3x3_c.h>
#include<invert2x2_c.h>
#include<transposed.h>
#include<stdio.h>

/*chan算法函数定义
Arguments    :
*                const double anchor_position[8]-->anchor坐标，此处用了4个anchor，输入数据时，依次输入{x1,y1,x2,y2,x3,y3,x4,y4}
*                const double ai_2_tag_minus_a1_2_tag
-->ai_2_tag_minus_a1_2_tag[0]:anchor2与tag距离  与   anchor1与tag距离 之差
*                                   -->ai_2_tag_minus_a1_2_tag[1]:anchor3与tag距离  与   anchor1与tag距离 之差
*                                   -->ai_2_tag_minus_a1_2_tag[2]:anchor4与tag距离  与   anchor1与tag距离 之差
*                                   在实际过程中，由到达时间差乘以光速算得
*                double tag_position[2]  -->输出tag经过chan算法求得的坐标
* Return Type  : void
*/
void chan_2D_algrithm(double* anchor_position, double* ai_2_tag_minus_a1_2_tag, double* tag_position)
{
	int i = 0, j = 0, k = 0;//计数用
	double Q[9] = { 1, 0.5, 0.5 , 0.5, 1, 0.5 , 0.5, 0.5,1 };
	double Q_inver[9] = { 0,0,0,0,0,0,0,0,0 };
	double Ga[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	double Anchor_pos[2][4] = { { 0,0,0,0 },{ 0,0,0,0 } };
	double h[3] = { 0,0,0 };
	double Ga_[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	//将输入的anchor坐标由一维数组转换为二维数组，方便后续使用
	for (i = 0; i < 4; i++)
	{
		Anchor_pos[0][i] = anchor_position[i * 2];
		Anchor_pos[1][i] = anchor_position[i * 2 + 1];
	}

	//m = (anchor_position(1, 1)) ^ 2 + (anchor_position(2, 1)) ^ 2
	double m = pow(anchor_position[0], 2) + pow(anchor_position[1], 2);

	/*for i = 1: num_of_anchor - 1
	K(i) = anchor_position(1, i + 1) ^ 2 + anchor_position(2, i + 1) ^ 2;
	end*/
	double K[3] = { 0,0,0 };
	for (i = 0; i < 3; ++i)
		K[i] = pow(anchor_position[i * 2 + 2], 2) + pow(anchor_position[i * 2 + 3], 2);

	/*for i = 1: num_of_anchor - 1
	Ga(i, 1) = -(anchor_position(1, i + 1) - anchor_position(1, 1));
	Ga(i, 2) = -(anchor_position(2, i + 1) - anchor_position(2, 1));
	Ga(i, 3) = -R(i);
	end*/
	for (i = 0; i < 3; i++)
	{
		Ga[i][0] = -(Anchor_pos[0][i + 1] - Anchor_pos[0][0]);
		Ga[i][1] = -(Anchor_pos[1][i + 1] - Anchor_pos[1][0]);
		Ga[i][2] = -ai_2_tag_minus_a1_2_tag[i];
	}/////////////////////

	/*for i = 1: num_of_anchor - 1
	h(i) = 0.5*(R(i) ^ 2 - K(i) + m);
	end*/
	for (i = 0; i < 3; i++)
		h[i] = 0.5 * (pow(ai_2_tag_minus_a1_2_tag[i], 2) - K[i] + m);

	/*Za0 = inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h';%(14b）*/
	//将上述步骤拆解开来一步步实现：1）求矩阵转置； 2）矩阵求逆； 3）矩阵相乘求解
	//1)求矩阵转置
	transposed(Ga, Ga_);/////////////

	//2)求矩阵求逆inv(Q)
	invert3x3(Q, Q_inver);

	//matrix_inverse(Q, Q_inver);
	// 3）矩阵相乘求解
	Matrix Ga_m;
	Ga_m.columns = 3;
	Ga_m.rows = 3;
	k = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			Ga_m.elements[k++] = Ga_[i][j];//////////////

	Matrix Q_inver_m;
	Q_inver_m.columns = 3;
	Q_inver_m.rows = 3;
	for (i = 0; i < 9; i++)
		Q_inver_m.elements[i] = Q_inver[i];

	Matrix Ga_Qinv_m;
	Ga_Qinv_m.columns = 3;
	Ga_Qinv_m.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_Qinv_m.elements[i] = 0;
	//Ga'*inv(Q)
	multiply(&Ga_m, &Q_inver_m, &Ga_Qinv_m);

	Matrix Gam;
	Gam.columns = 3;
	Gam.rows = 3;
	k = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			Gam.elements[k++] = Ga[i][j];

	Matrix Ga_mQinv_mGam;
	Ga_mQinv_mGam.columns = 3;
	Ga_mQinv_mGam.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_mQinv_mGam.elements[i] = 0;
	//Ga'*inv(Q)*Ga	
	multiply(&Ga_Qinv_m, &Gam, &Ga_mQinv_mGam);

	double Ga_mQinv_mGam_inv[9] = { 0 };
	invert3x3(Ga_mQinv_mGam.elements, Ga_mQinv_mGam_inv);

	Matrix Ga_mQinv_mGam_inv_m;
	Ga_mQinv_mGam_inv_m.columns = 3;
	Ga_mQinv_mGam_inv_m.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_mQinv_mGam_inv_m.elements[i] = Ga_mQinv_mGam_inv[i];

	Matrix Ga_mQinv_mGam_inv_mGa_;
	Ga_mQinv_mGam_inv_mGa_.columns = 3;
	Ga_mQinv_mGam_inv_mGa_.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_mQinv_mGam_inv_mGa_.elements[i] = 0;
	//inv(Ga'*inv(Q)*Ga)*Ga'
	multiply(&Ga_mQinv_mGam_inv_m, &Ga_m, &Ga_mQinv_mGam_inv_mGa_);

	Matrix Ga_mQinv_mGam_inv_mGa_Qinv;
	Ga_mQinv_mGam_inv_mGa_Qinv.columns = 3;
	Ga_mQinv_mGam_inv_mGa_Qinv.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_mQinv_mGam_inv_mGa_Qinv.elements[i] = 0;
	//inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)
	multiply(&Ga_mQinv_mGam_inv_mGa_, &Q_inver_m, &Ga_mQinv_mGam_inv_mGa_Qinv);

	Matrix hm;
	hm.columns = 1;
	hm.rows = 3;
	for (i = 0; i < 3; i++)
		hm.elements[i] = h[i];

	Matrix Za0;
	Za0.columns = 1;
	Za0.rows = 3;
	for (i = 0; i < 3; i++)
		Za0.elements[i] = 0;
	//Za0 = inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h'
	multiply(&Ga_mQinv_mGam_inv_mGa_Qinv, &hm, &Za0);

	/*B = eye(num_of_anchor - 1);
	for i = 1: num_of_anchor - 1
		B(i, i) = sqrt((anchor_position(1, i + 1) - Za0(1)) ^ 2 + (anchor_position(2, i + 1) - Za0(2)) ^ 2);
	end*/
	double C[9] = { 0,0,0,0,0,0,0,0,0 };
	double B[9] = { 0,0,0,0,0,0,0,0,0 };
	for (i = 0; i < 3; i++)
		C[i] = sqrt(pow(anchor_position[2 * i + 2] - Za0.elements[0], 2) + pow(anchor_position[2 * i + 3] - Za0.elements[1], 2));
	k = 0;
	for (i = 0; i < 3; i++)
	{
		B[k] = C[i];
		k = k + 4;
	}

	//FI = B * Q * B;
	Matrix Bm;
	Bm.columns = 3;
	Bm.rows = 3;
	for (i = 0; i < 9; i++)
		Bm.elements[i] = B[i];

	Matrix Qm;
	Qm.columns = 3;
	Qm.rows = 3;
	for (i = 0; i < 9; i++)
		Qm.elements[i] = Q[i];

	Matrix BmQm;
	BmQm.columns = 3;
	BmQm.rows = 3;
	for (i = 0; i < 9; i++)
		BmQm.elements[i] = 0;

	multiply(&Bm, &Qm, &BmQm);

	Matrix FIm;
	FIm.columns = 3;
	FIm.rows = 3;
	for (i = 0; i < 9; i++)
		FIm.elements[i] = 0;

	multiply(&BmQm, &Bm, &FIm);

	//Za1 = inv(Ga'*inv(FI)*Ga)*Ga'*inv(FI)*h';%（14a）
	//先求inv（FI）
	double FI_inv[9] = { 0 };
	invert3x3(FIm.elements, FI_inv);
	//Ga'*inv(FI)
	Matrix FI_invm;
	FI_invm.columns = 3;
	FI_invm.rows = 3;
	for (i = 0; i < 9; i++)
		FI_invm.elements[i] = FI_inv[i];

	Matrix Ga_FI_invm;
	Ga_FI_invm.columns = 3;
	Ga_FI_invm.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_FI_invm.elements[i] = 0;
	multiply(&Ga_m, &FI_invm, &Ga_FI_invm);

	//Ga'*inv(FI)*Ga
	Matrix Ga_FI_invmGam;
	Ga_FI_invmGam.columns = 3;
	Ga_FI_invmGam.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_FI_invmGam.elements[i] = 0;
	multiply(&Ga_FI_invm, &Gam, &Ga_FI_invmGam);

	//inv(Ga'*inv(FI)*Ga)
	double Ga_FI_invmGam_inv[9] = { 0 };
	invert3x3(Ga_FI_invmGam.elements, Ga_FI_invmGam_inv);

	//inv(Ga'*inv(FI)*Ga)*Ga'
	Matrix Ga_FI_invmGam_invm;
	Ga_FI_invmGam_invm.columns = 3;
	Ga_FI_invmGam_invm.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_FI_invmGam_invm.elements[i] = Ga_FI_invmGam_inv[i];

	Matrix Ga_FI_invmGam_invmGa_m;
	Ga_FI_invmGam_invmGa_m.columns = 3;
	Ga_FI_invmGam_invmGa_m.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_FI_invmGam_invmGa_m.elements[i] = 0;

	multiply(&Ga_FI_invmGam_invm, &Ga_m, &Ga_FI_invmGam_invmGa_m);

	//inv(Ga'*inv(FI)*Ga)*Ga'*inv(FI)
	Matrix Ga_FI_invmGam_invmGa_mFI_invm;
	Ga_FI_invmGam_invmGa_mFI_invm.columns = 3;
	Ga_FI_invmGam_invmGa_mFI_invm.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_FI_invmGam_invmGa_mFI_invm.elements[i] = 0;

	multiply(&Ga_FI_invmGam_invmGa_m, &FI_invm, &Ga_FI_invmGam_invmGa_mFI_invm);

	//Za1 = inv(Ga'*inv(FI)*Ga)*Ga'*inv(FI)*h';
	Matrix Za1;
	Za1.columns = 3;
	Za1.rows = 3;
	for (i = 0; i < 9; i++)
		Za1.elements[i] = 0;

	multiply(&Ga_FI_invmGam_invmGa_mFI_invm, &hm, &Za1);

	//if Za1(3) < 0 % 此处需要理解一下这的判别条件
	//	Za1(3) = abs(Za1(3));
	//end
	if (Za1.elements[2] < 0)
		Za1.elements[2] = -Za1.elements[2];

	//   CovZa = inv(Ga'*inv(FI)*Ga);
	Matrix CovZa = Ga_FI_invmGam_invm;

	//sB = eye(3);
	//sB(1, 1) = Za1(1) - anchor_position(1, 1);
	//sB(2, 2) = Za1(2) - anchor_position(2, 1);
	//sB(3, 3) = Za1(3);
	Matrix sB;
	sB.columns = 3;
	sB.rows = 3;
	for (i = 0; i < 9; i++)
		sB.elements[i] = 0;
	sB.elements[0] = Za1.elements[0] - anchor_position[0];
	sB.elements[4] = Za1.elements[1] - anchor_position[1];
	sB.elements[8] = Za1.elements[2];

	Matrix fourtimessB;
	fourtimessB.columns = 3;
	fourtimessB.rows = 3;
	for (i = 0; i < 9; i++)
		fourtimessB.elements[i] = 0;
	//sFI = 4 * sB * CovZa * sB;
	// 4 * sB
	for (i = 0; i < 9; i++)
		fourtimessB.elements[i] = 4 * sB.elements[i];
	//4 * sB * CovZa
	Matrix fourtimessBCovZam;
	fourtimessBCovZam.columns = 3;
	fourtimessBCovZam.rows = 3;
	for (i = 0; i < 9; i++)
		fourtimessBCovZam.elements[i] = 0;

	multiply(&fourtimessB, &CovZa, &fourtimessBCovZam);

	Matrix sFIm;
	sFIm.columns = 3;
	sFIm.rows = 3;
	for (i = 0; i < 9; i++)
		sFIm.elements[i] = 0;

	//sFI =;
	multiply(&fourtimessBCovZam, &sB, &sFIm);

	//sGa = [1, 0; 0, 1; 1, 1];
	Matrix sGam;
	sGam.columns = 2;
	sGam.rows = 3;
	for (i = 0; i < 6; i++)
		sGam.elements[i] = 1;
	sGam.elements[1] = 0;
	sGam.elements[2] = 0;

	//    sh  = [(Za1(1)-anchor_position(1, 1))^2; (Za1(2)-anchor_position(2, 1))^2; Za1(3)^2];%x1,y1,在此处也进行了运算；
	Matrix shm;
	shm.columns = 1;
	shm.rows = 3;
	for (i = 0; i < 3; i++)
		shm.elements[i] = 0;
	shm.elements[0] = pow(Za1.elements[0] - anchor_position[0], 2);
	shm.elements[1] = pow(Za1.elements[1] - anchor_position[1], 2);
	shm.elements[2] = pow(Za1.elements[2], 2);

	// Za2 = inv(sGa'*inv(sFI)*sGa)*sGa'*inv(sFI)*sh;%这个是论文当中的（22a）
	//inv(sFI)  
	Matrix sFIm_invm;
	sFIm_invm.columns = 3;
	sFIm_invm.rows = 3;
	for (i = 0; i < 9; i++)
		sFIm_invm.elements[i] = 0;

	invert3x3(sFIm.elements, sFIm_invm.elements);
	//sGa'
	Matrix sGa_m;
	sGa_m.columns = 3;
	sGa_m.rows = 2;
	for (i = 0; i < 6; i++)
		sGa_m.elements[i] = 1;
	sGa_m.elements[1] = 0;
	sGa_m.elements[3] = 0;
	//sGa'*inv(sFI)
	Matrix sGa_msFIm_invmm;
	sGa_msFIm_invmm.columns = 3;
	sGa_msFIm_invmm.rows = 2;
	for (i = 0; i < 3; i++)
		sGa_msFIm_invmm.elements[i] = 0;
	multiply(&sGa_m, &sFIm_invm, &sGa_msFIm_invmm);

	//sGa'*inv(sFI)*sGa
	Matrix sGa_msFIm_invmmsGam;
	sGa_msFIm_invmmsGam.columns = 3;
	sGa_msFIm_invmmsGam.rows = 2;
	for (i = 0; i < 4; i++)
		sGa_msFIm_invmmsGam.elements[i] = 0;
	multiply(&sGa_msFIm_invmm, &sGam, &sGa_msFIm_invmmsGam);
	//inv(sGa'*inv(sFI)*sGa)
	Matrix sGa_msFIm_invmmsGam_invm;
	sGa_msFIm_invmmsGam_invm.columns = 2;
	sGa_msFIm_invmmsGam_invm.rows = 2;
	for (i = 0; i < 4; i++)
		sGa_msFIm_invmmsGam_invm.elements[i] = 0;
	invert2x2(sGa_msFIm_invmmsGam.elements, sGa_msFIm_invmmsGam_invm.elements);
	//inv(sGa'*inv(sFI)*sGa)*sGa'
	Matrix sGa_msFIm_invmmsGam_invmsGa_m;
	sGa_msFIm_invmmsGam_invmsGa_m.columns = 3;
	sGa_msFIm_invmmsGam_invmsGa_m.rows = 2;
	for (i = 0; i < 4; i++)
		sGa_msFIm_invmmsGam_invmsGa_m.elements[i] = 0;
	multiply(&sGa_msFIm_invmmsGam_invm, &sGa_m, &sGa_msFIm_invmmsGam_invmsGa_m);
	//inv(sGa'*inv(sFI)*sGa)*sGa'*inv(sFI)
	Matrix sGa_msFIm_invmmsGam_invmsGa_minvsFIm;
	sGa_msFIm_invmmsGam_invmsGa_minvsFIm.columns = 3;
	sGa_msFIm_invmmsGam_invmsGa_minvsFIm.rows = 2;
	for (i = 0; i < 4; i++)
		sGa_msFIm_invmmsGam_invmsGa_minvsFIm.elements[i] = 0;
	multiply(&sGa_msFIm_invmmsGam_invmsGa_m, &sFIm_invm, &sGa_msFIm_invmmsGam_invmsGa_minvsFIm);

	//Za2 = inv(sGa'*inv(sFI)*sGa)*sGa'*inv(sFI)*sh
	Matrix Za2;
	Za2.columns = 1;
	Za2.rows = 2;
	for (i = 0; i < 2; i++)
		Za2.elements[i] = 0;
	multiply(&sGa_msFIm_invmmsGam_invmsGa_minvsFIm, &shm, &Za2);
	//Za = sqrt(abs(Za2));
	Matrix Za;
	Za.columns = 1;
	Za.rows = 2;
	for (i = 0; i < 2; i++)
		Za.elements[i] = sqrt(abs(Za2.elements[i]));

	// reout(1)=anchor_position(1, 1);
	//reout(2) = anchor_position(2, 1);
	//out = abs(Za) + reout;
	tag_position[0] = abs(Za.elements[0]) + anchor_position[0];
	tag_position[1] = abs(Za.elements[1]) + anchor_position[1];
	if (Za1.elements[0] < anchor_position[0])
		tag_position[0] = -abs(Za.elements[0]) + anchor_position[0];
	if (Za1.elements[1] < anchor_position[1])
		tag_position[1] = -abs(Za.elements[1]) + anchor_position[1];

	for (i = 0; i < 2; i++)
		printf("%f \n\n", Za.elements[i]);

}

#endif
