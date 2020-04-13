#include<stdio.h>
#include<string.h>
#include<iostream>
#include<math.h>
//chan算法函数声明
#define N 3
void chan_2D_algrithm(int number_of_anchor, double * anchor_position, double *  ai_2_tag_minus_a1_2_tag, double * tag_position);
void transposed(double Ga[3][3], double(*Ga_)[3]);
void matrix_inverse(double a[N][N], double b[N][N]);


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
	double Q[3][3] = { { 1, 0.5, 0.5 },{ 0.5, 1, 0.5 },{ 0.5, 0.5, 1 } };
	double Q_inver[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };
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
	//2)求矩阵求逆inv(Q)
	//matrix_inverse( Q, Q_inver);
	matrix_inverse(Ga, Ga_);


	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			printf("%f \n\n", Ga[i][j]);
		}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			printf("%f \n\n", Ga_[i][j]);
		}
	/*for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			printf("%f \n\n", Q_inver[i][j]);
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

//求解3*3矩阵的逆
void matrix_inverse(double a[N][N], double b[N][N])
{
	int i=0, j=0, k=0;
	double max=0, temp=0;
	// 定义一个临时矩阵t
	double t[N][N];
	// 将a矩阵临时存放在矩阵t[n][n]中
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			t[i][j] = a[i][j];
		}
	}
	// 初始化B矩阵为单位矩阵
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			b[i][j] = (i == j) ? 1 : 0;//此处删除了1之前的（float）
		}
	}
	// 进行列主消元，找到每一列的主元
	for (i = 0; i<N; i++)
	{
		max = t[i][i];
		// 用于记录每一列中的第几个元素为主元
		k = i;
		// 寻找每一列中的主元元素
		for (j = i + 1; j<N; j++)
		{
			if (fabsf(t[j][i]) > fabsf(max))
			{
				max = t[j][i];
				k = j;
			}
		}
		//cout<<"the max number is "<<max<<endl;
		// 如果主元所在的行不是第i行，则进行行交换
		if (k != i)
		{
			// 进行行交换
			for (j = 0; j<N; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;
				// 伴随矩阵B也要进行行交换
				temp = b[i][j];
				b[i][j] = b[k][j];
				b[k][j] = temp;
			}
		}
		if (t[i][i] == 0)
		{
			printf( "\nthe matrix does not exist inverse matrix\n");
			break;
		}
		// 获取列主元素
		temp = t[i][i];
		// 将主元所在的行进行单位化处理
		//cout<<"\nthe temp is "<<temp<<endl;
		for (j = 0; j<N; j++)
		{
			t[i][j] = t[i][j] / temp;
			b[i][j] = b[i][j] / temp;
		}
		for (j = 0; j<N; j++)
		{
			if (j != i)
			{
				temp = t[j][i];
				//消去该列的其他元素
				for (k = 0; k<N; k++)
				{
					t[j][k] = t[j][k] - temp*t[i][k];
					b[j][k] = b[j][k] - temp*b[i][k];
				}
			}

		}

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
	printf("计算   9999结果是X=%f Y=%f\n\n", tag_position_[0], tag_position_[1]);
	system("pause");
	return 0;
	
}