#include<stdio.h>
#include<string.h>
#include<iostream>
#include<math.h>
//chan算法函数声明
void chan_2D_algrithm(int number_of_anchor, double * anchor_position, double *  ai_2_tag_minus_a1_2_tag, double * tag_position);
void transposed(double Ga[3][3], double(*Ga_)[3]);

/*chan算法函数定义
Arguments    : double num_of_anchor   -->anchor数量
*                const double anchor_position[8]-->anchor坐标，此处用了4个anchor，输入数据时，依次输入{x1,y1,x2,y2,x3,y3,x4,y4}
*                const double ai_2_tag_minus_a1_2_tag
-->ai_2_tag_minus_a1_2_tag[0]:anchor2与tag距离  与   anchor1与tag距离 之差
*                                   -->ai_2_tag_minus_a1_2_tag[1]:anchor3与tag距离  与   anchor1与tag距离 之差
*                                   -->ai_2_tag_minus_a1_2_tag[2]:anchor4与tag距离  与   anchor1与tag距离 之差
*                                   在实际过程中，由到达时间差乘以光速算得
*                double tag_position[2]  -->输出tag经过chan算法求得的坐标
* Return Type  : void
*/
void chan_2D_algrithm(int number_of_anchor, double * anchor_position, double *  ai_2_tag_minus_a1_2_tag, double * tag_position)
{
	/*tag_position[0]= 1;
	tag_position[1] = 15;
*/
	int i = 0,j=0;//计数用
	double Q[9] = { 1, 0.5, 0.5,0.5, 1, 0.5,0.5, 0.5, 1 };
	double Ga[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
	double Anchor_pos[2][4]= { { 0,0,0,0 },{ 0,0,0,0 } };
	double h[3] = {0,0,0};
	double Ga_[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
	//将输入的anchor坐标由一维数组转换为二维数组，方便后续使用
	for (i = 0; i < 4; i++)
	{
		Anchor_pos[0][i] = anchor_position[i * 2];
		Anchor_pos[1][i] = anchor_position[i * 2 + 1];
	}

	//m = (anchor_position(1, 1)) ^ 2 + (anchor_position(2, 1)) ^ 2
	double m = pow(anchor_position[0],2) + pow(anchor_position[1],2);

	/*for i = 1: num_of_anchor - 1
		K(i) = anchor_position(1, i + 1) ^ 2 + anchor_position(2, i + 1) ^ 2;
	end*/
	double K[3] = { 0,0,0 };
	for (i = 0; i < 3;++i)
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

	/*Za0 = inv(Ga'*inv(Q)*Ga)*Ga'*inv(Q)*h';%(14b）*/
	//将上述步骤拆解开来一步步实现：1）求矩阵转置； 2）矩阵求逆； 3）矩阵相乘求解
	//1)求矩阵转置
	transposed(Ga,Ga_);
/*
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			printf("%f \n\n", Ga_[i][j]);
		}*/
	
	//printf("%f", m);

}

//3*3矩阵求转置的函数定义
//此处学到一点，c语言的函数不能直接返回二维数组
void transposed(double Ga[3][3],double (*Ga_)[3])
{
	int i = 0, j = 0;//计数用
	for(i=0;i<3;i++)
		for (j = 0; j < 3; j++)
		{
			Ga_[j][i]=Ga[i][j];
		}
	
}

//主函数入口
int main()
{
	//定义chan算法所需参数，各参数含义如其名所述
	int number_of_anchor = 4;
	double anchor_position_[8] = { 2,10,10,20,25,0,30,20 };
	double ai_2_tag_minus_a1_2_tag[3] = { 2.81024967,19.41311123,20.70992026 };
	double tag_position_[2] = { 0,0 };


	//运行chan算法
	chan_2D_algrithm(number_of_anchor, anchor_position_, ai_2_tag_minus_a1_2_tag, tag_position_);
	//输出运行后的tag坐标结果
	printf("计算结果是X=%f Y=%f\n\n", tag_position_[0], tag_position_[1]);
	system("pause");
	return 0;
	
}