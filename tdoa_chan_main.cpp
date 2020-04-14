#include<stdio.h>
#include<string.h>
#include<iostream>
#include<math.h>
#include<matrixmul.h>
#include<invert3x3_c.h>
//chan�㷨��������
#define N 3
void chan_2D_algrithm(int number_of_anchor, double * anchor_position, double *  ai_2_tag_minus_a1_2_tag, double * tag_position);
void transposed(double Ga[3][3], double(*Ga_)[3]);
//void matrix_inverse(double a[N][N], double b[N][N]);
//����һ���Ƿ񱣴���ȷ===

/*chan�㷨��������
Arguments    : double num_of_anchor   -->anchor����
*                const double anchor_position[8]-->anchor���꣬�˴�����4��anchor����������ʱ����������{x1,y1,x2,y2,x3,y3,x4,y4}
*                const double ai_2_tag_minus_a1_2_tag
-->ai_2_tag_minus_a1_2_tag[0]:anchor2��tag����  ��   anchor1��tag���� ֮��
*                                   -->ai_2_tag_minus_a1_2_tag[1]:anchor3��tag����  ��   anchor1��tag���� ֮��
*                                   -->ai_2_tag_minus_a1_2_tag[2]:anchor4��tag����  ��   anchor1��tag���� ֮��
*                                   ��ʵ�ʹ����У��ɵ���ʱ�����Թ������
*                double tag_position[2]  -->���tag����chan�㷨��õ�����
* Return Type  : void
*/
void chan_2D_algrithm(int number_of_anchor, double * anchor_position, double *  ai_2_tag_minus_a1_2_tag, double * tag_position)
{
	/*tag_position[0]= 1;
	tag_position[1] = 15;
	*/
	int i = 0, j = 0,k=0;//������
	/*double Q[3][3] = { { 1, 0.5, 0.5 },{ 0.5, 1, 0.5 },{ 0.5, 0.5, 1 } };
	double Q_inver[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };*/
	double Q[9] = { 1, 0.5, 0.5 , 0.5, 1, 0.5 , 0.5, 0.5,1 };
	double Q_inver[9] = { 0,0,0,0,0,0,0,0,0 };
	double Ga[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	double Anchor_pos[2][4] = { { 0,0,0,0 },{ 0,0,0,0 } };
	double h[3] = { 0,0,0 };
	double Ga_[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	//double Ga_Qinv[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	//�������anchor������һά����ת��Ϊ��ά���飬�������ʹ��
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
	{
		K[i] = pow(anchor_position[i * 2 + 2], 2) + pow(anchor_position[i * 2 + 3], 2);

	}

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
	}

	/*for i = 1: num_of_anchor - 1
	h(i) = 0.5*(R(i) ^ 2 - K(i) + m);
	end*/
	for (i = 0; i < 3; i++)
	{
		h[i] = 0.5*(pow(ai_2_tag_minus_a1_2_tag[i], 2) - K[i] + m);
	}

	/*Za0 = inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h';%(14b��*/
	//�����������⿪��һ����ʵ�֣�1�������ת�ã� 2���������棻 3������������
	//1)�����ת��
	transposed(Ga, Ga_);
	
	//2)���������inv(Q)
	invert3x3(Q, Q_inver);
	
	//matrix_inverse(Q, Q_inver);
	// 3������������
	Matrix Ga_m;
	Ga_m.columns = 3;
	Ga_m.rows = 3;
	k = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			Ga_m.elements[k++] = Ga_[i][j];
		}

	Matrix Q_inver_m;
	Q_inver_m.columns = 3;
	Q_inver_m.rows = 3;
	for (i = 0; i < 9; i++)
			Q_inver_m.elements[i] = Q_inver[i];

	Matrix Ga_Qinv_m;
	Ga_Qinv_m.columns = 3;
	Ga_Qinv_m.rows = 3;
	for (i = 0; i < 9; i++)
			Ga_Qinv_m.elements[i] =0;
	//Ga'*inv(Q)
	multiply(&Ga_m, &Q_inver_m, &Ga_Qinv_m);

	Matrix Gam;
	Gam.columns = 3;
	Gam.rows = 3;
	k = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			Gam.elements[k++] = Ga[i][j];
		}

	Matrix Ga_mQinv_mGam;
	Ga_mQinv_mGam.columns = 3;
	Ga_mQinv_mGam.rows = 3;
	for (i = 0; i < 9; i++)
		Ga_mQinv_mGam.elements[i] = 0;
	//Ga'*inv(Q)*Ga	
	multiply(&Ga_Qinv_m,&Gam, &Ga_mQinv_mGam);
	
	double Ga_mQinv_mGam_inv[9] = {0};
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
	{
		C[i] = sqrt(pow(anchor_position[2 * i + 2] - Za0.elements[0], 2) + pow(anchor_position[2 * i + 3] - Za0.elements[1], 2));
	}
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

	//Za1 = inv(Ga'*inv(FI)*Ga)*Ga'*inv(FI)*h';%��14a��
	//����inv��FI��
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

	//if Za1(3) < 0 % �˴���Ҫ���һ������б�����
	//	Za1(3) = abs(Za1(3));
	//end
	if (Za1.elements[2] < 0)
		Za1.elements[2] = -Za1.elements[2];

	//   CovZa = inv(Ga'*inv(FI)*Ga);
	Matrix CovZa= Ga_FI_invmGam_invm;

	//sB = eye(3);
	//sB(1, 1) = Za1(1) - anchor_position(1, 1);
	//sB(2, 2) = Za1(2) - anchor_position(2, 1);
	//sB(3, 3) = Za1(3);
	Matrix sB;
	sB.columns = 3;
	sB.rows = 3;
	for (i = 0; i < 9; i++)
	{
		sB.elements[i] = 0;
	}
	sB.elements[0] = Za1.elements[0] - anchor_position[0];
	sB.elements[4] = Za1.elements[1] - anchor_position[1];
	sB.elements[8] = Za1.elements[2] ;

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
	{
		fourtimessBCovZam.elements[i] = 0;
	}
	multiply(&fourtimessB, &CovZa, &fourtimessBCovZam);

	Matrix sFIm;
	sFIm.columns = 3;
	sFIm.rows = 3;
	for (i = 0; i < 9; i++)
	{
		sFIm.elements[i] = 0;
	}
	//sFI =;
	multiply(&fourtimessBCovZam, &sB, &sFIm);

	//sGa = [1, 0; 0, 1; 1, 1];
	Matrix sGam;
	sGam.columns = 2;
	sGam.rows = 3;
	for (i = 0; i < 6; i++)
	{
		sGam.elements[i] = 1;
	}
	sGam.elements[1] = 0;
	sGam.elements[2] = 0;

	//    sh  = [(Za1(1)-anchor_position(1, 1))^2; (Za1(2)-anchor_position(2, 1))^2; Za1(3)^2];%x1,y1,�ڴ˴�Ҳ���������㣻
	Matrix shm;
	shm.columns = 1;
	shm.rows = 3;
	for (i = 0; i < 3; i++)
	{
		shm.elements[i] = 0;
	}
	shm.elements[0] = pow(Za1.elements[0] - anchor_position[0], 2);
	shm.elements[1] = pow(Za1.elements[1] - anchor_position[1], 2);
	shm.elements[2] = pow(Za1.elements[2] , 2);

	for (i = 0; i <3; i++)
		printf("%f \n\n", shm.elements[i]);
	
	//for (i = 0; i < 9; i++)
	//	printf("%d %d %f \n\n", FIm.columns, FIm.rows, FIm.elements[i]);
	//

	//for (i = 0; i < 3; i++)
	//	for (j = 0; j < 3; j++)
	//	{
	//		printf("%f \n\n\n\n\n\n\n\n", Ga[i][j]);
	//	}
	//for (i = 0; i < 3; i++)
	//	for (j = 0; j < 3; j++)
	//	{
	//		printf("%f \n\n", Ga_[i][j]);
	//	}
	//printf("%f", m);
}

//3*3������ת�õĺ�������
//�˴�ѧ��һ�㣬c���Եĺ�������ֱ�ӷ��ض�ά����
void transposed(double Ga[3][3], double(*Ga_)[3])
{
	int i = 0, j = 0;//������
	for (i = 0; i<3; i++)
		for (j = 0; j < 3; j++)
		{
			Ga_[j][i] = Ga[i][j];
		}

}



//���������
int main()
{
	//����chan�㷨�����������������������������
	int number_of_anchor = 4;
	double anchor_position_[8] = { 2,10,10,20,25,0,30,20 };
	double ai_2_tag_minus_a1_2_tag[3] = { 2.81024967,19.41311123,20.70992026 };
	double tag_position_[2] = { 0,0 };


	//����chan�㷨
	chan_2D_algrithm(number_of_anchor, anchor_position_, ai_2_tag_minus_a1_2_tag, tag_position_);
	//������к��tag������
	printf("��������X=%f Y=%f\n\n", tag_position_[0], tag_position_[1]);
	system("pause");
	return 0;

}