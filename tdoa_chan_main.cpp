#include<stdio.h>
#include<string.h>
#include<iostream>
#include<chan_2D.h>

//���������
int main()
{
	//����chan�㷨�����������������������������
	int number_of_anchor = 4;
	double anchor_position_[8] = { 2,10,10,20,25,0,30,20 };
	double ai_2_tag_minus_a1_2_tag[3] = { 2.81024967,19.41311123,20.70992026 };
	double tag_position_[2] = { 0,0 };

	//����chan�㷨
	chan_2D_algrithm( anchor_position_, ai_2_tag_minus_a1_2_tag, tag_position_);

	//������к��tag������
	printf("��������X=%f Y=%f\n\n", tag_position_[0], tag_position_[1]);
	system("pause");
	return 0;

}