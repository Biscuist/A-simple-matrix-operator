#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include"matrix.h"
#pragma warning(disable:4996)
int column, line;
int main()
{
	int i, input;
	double value;
	double* origin;
	while (true)
	{
		system("cls");
		printf("->按1计算一个行列式的值\n");
		printf("->按2将矩阵化为行最简形\n");
		printf("->按3对一个基进行施密特正交化\n");
		printf("按其他数字退出程序\n");
		scanf("%d", &input);
		system("cls");
		switch (input)
		{
		case 1:	printf("请输入行列式的阶数：\n");
			scanf("%d", &line);
			system("cls");
			origin = (double*)calloc(line * line, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("输入行列式的元素：\n");
			for (i = 0; i < line * line; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("你输入的行列式为:\n");
			show_matrix(line, line, origin);
			value = det_value(line, origin);
			printf("化简后的行列式为:\n");
			show_matrix(line, line, origin);
			printf("该行列式的值为:%lf\n", value);
			free(origin);
			system("pause");
			break;
		case 2:	printf("依次输入矩阵的行数和列数：\n");
			scanf("%d %d", &column, &line);
			system("cls");
			origin = (double*)calloc(line * column, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("输入矩阵的元素：\n");
			for (i = 0; i < line * column; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("你输入的矩阵为:\n");
			show_matrix(column, line, origin);
			simplify_matrix(column, line, origin, origin);
			printf("化简后的矩阵为:\n");
			show_matrix(column, line, origin);
			free(origin);
			system("pause");
			break;
		case 3: printf("依次输入基的行数和列数：\n");
			scanf("%d %d", &column, &line);
			system("cls");
			origin = (double*)calloc(line * column, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("输入基的元素：\n");
			for (i = 0; i < line * column; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("你输入的基为:\n");
			show_matrix(column, line, origin);
			value = Schmidt_orthogonalization(column, line, origin);
			if (value == true)
			{
				printf("施密特正交化后的标准正交基为:\n");
				show_matrix(column, line, origin);
			}
			free(origin);
			system("pause");
			break;
		default: return 0;
		}
	}
	return 0;
}