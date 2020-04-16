#ifndef transposed_h
#define transposed_h

void transposed(double Ga[3][3], double(*Ga_)[3])
{
	int i = 0, j = 0;//¼ÆÊıÓÃ
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			Ga_[j][i] = Ga[i][j];
		}

}

#endif