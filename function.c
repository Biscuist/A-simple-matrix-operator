#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include"matrix.h"
#pragma warning(disable:4996)
extern column,line ;
void simplify_matrix(int t_column, int t_line, double* head, double* origin)//head为每行第一个非0元素
{
	if (t_column == 0 || t_line == 0)
		return;
	int i = 0, j = 0;
	double* pmember = head;
	while (fabs(pmember[i * line] - 0) < (1e-6) && t_line > 0)//定位首元素不为零的行
	{
		i++;
		if (i == t_column)                            //如果一列都为0，横着移动一格
		{
			t_line = t_line - 1;
			if (t_line == 0)
				return;
			head = head + 1;
			pmember = head;
			i = 0;
		}
	}
	/*if (t_column == 0 || t_line == 0)
		return;*/
	pmember = pmember + i * line;
	double* temp_cloumn = (double*)calloc(t_line, sizeof(double));
	if (temp_cloumn == NULL)
		exit(false);
	memmove(temp_cloumn, head, t_line * sizeof(double));
	memmove(head, pmember, t_line * sizeof(double));					//将首元素不为0的行移动到“第一行”
	memmove(pmember, temp_cloumn, t_line * sizeof(double));
	pmember = head;
	double tmp = pmember[0];											//保存“第一行”首非零元素
	for (i = 0; i < t_line; i++)
	{
		pmember[i] = pmember[i] / tmp;
	}

	for (j = 0; j < column; j++)
	{
		if (origin + (j * line + line - t_line) == head)            //origin[j * line + line - t_line]为当前矩阵的首元素
			continue;
		tmp = origin[j * line + line - t_line];					//保存当前行的第一个元素
		for (i = 0; i < t_line; i++)
		{
			origin[j * line + i + line - t_line] = origin[j * line + i + line - t_line] - tmp * head[i];
		}
	}
	simplify_matrix(t_column - 1, t_line - 1, head + line + 1, origin);
	free(temp_cloumn);
}

void show_matrix(int column, int line, double* origin)
{
	int i, j;
	for (j = 0; j < column; j++)
	{
		for (i = 0; i < line; i++)
		{
			printf("%-7.2lf ", origin[j * line + i]);
		}
		putchar('\n');
	}
}

double det_value(int t_line, double* head)
{
	if (t_line == 1)
		return *head;
	int i = 0, j = 0;
	bool flag = 0;												//标识是否交换了两行
	double tmp;
	double* pmember = head;
	while (fabs(pmember[i * line] - 0) < (1e-6) && t_line > 0)//定位首元素不为零的行
	{
		i++;
		if (i == t_line)                                 //如果一列都为0，行列式为0
		{
			return 0.0;
		}
	}
	if (i != 0)
	{
		flag = 1;
		pmember = pmember + i * line;
		double* temp_cloumn = (double*)calloc(t_line, sizeof(double));
		if (temp_cloumn == NULL)
			exit(false);
		memmove(temp_cloumn, head, t_line * sizeof(double));
		memmove(head, pmember, t_line * sizeof(double));     //将首元素不为0的行移动到“第一行”
		memmove(pmember, temp_cloumn, t_line * sizeof(double));
		free(temp_cloumn);
		pmember = head;
	}
	for (j = 1; j < t_line; j++)
	{
		tmp = head[j * line];							//保存当前行的第一个元素
		for (i = 0; i < t_line; i++)
		{
			head[j * line + i] = head[j * line + i] - (tmp / head[0]) * head[i];
		}
	}

	return pow(-1, flag) * (*head) * det_value(t_line - 1, head + line + 1);
}
bool Schmidt_orthogonalization(int column, int line, double* origin)
{
	int i, j, k, l, r_matrix = 0;												//r_matrix为矩阵的秩
	double sum = 0;
	double molecule = 0, denominator = 0;											//这两个是公式中的分子和分母
	double* copy_matrix1 = (double*)calloc(column * line, sizeof(double));		//复制一个矩阵用于初等变换求秩
	if (copy_matrix1 == NULL)
		exit(false);
	memcpy(copy_matrix1, origin, line * column * sizeof(double));
	simplify_matrix(column, line, copy_matrix1, copy_matrix1);
	for (j = 0; j < column; j++)
	{
		if (j > r_matrix)
			break;
		for (i = 0; i < line; i++)
		{
			if (fabs(copy_matrix1[j * line + i] - 0) > (1e-6))
			{
				r_matrix++;
				break;
			}
		}
	}
	free(copy_matrix1);
	if (r_matrix < line)
	{
		printf("这就不是一个线性无关组\n");
		return false;
	}
	//正交化
	double* copy_matrix2 = (double*)calloc(column * line, sizeof(double));							//公式里面会用到这个矩阵
	if (copy_matrix2 == NULL)
		exit(false);
	memcpy(copy_matrix2, origin, line * column * sizeof(double));
	for (j = 1; j < line; j++)																		//变换每一列
	{
		for (i = 0; i < column; i++)																//每一列的每一行做相同的变换
		{
			for (k = 0, sum = 0; k < j; k++)														//不同列有着不同长度的加和
			{
				for (l = 0, molecule = 0; l < column; l++)											//计算点乘
				{
					molecule = molecule + copy_matrix2[j + l * line] * origin[k + l * line];			//____________________________________
				}
				for (l = 0, denominator = 0; l < column; l++)
				{
					denominator = denominator + origin[k + l * line] * origin[k + l * line];
				}
				sum = sum + molecule / denominator * origin[k + i * line];
			}
			origin[j + i * line] = origin[j + i * line] - sum;
		}
	}
	free(copy_matrix2);
	//单位化
	for (j = 0; j < line; j++)
	{
		for (i = 0, denominator = 0; i < column; i++)
		{
			denominator = denominator + origin[j + i * line] * origin[j + i * line];
		}
		denominator = sqrt(denominator);
		for (i = 0; i < column; i++)
		{
			origin[j + i * line] = origin[j + i * line] / denominator;
		}
	}
	return true;
}