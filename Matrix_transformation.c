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
		printf("->��1����һ������ʽ��ֵ\n");
		printf("->��2������Ϊ�������\n");
		printf("->��3��һ��������ʩ����������\n");
		printf("�����������˳�����\n");
		scanf("%d", &input);
		system("cls");
		switch (input)
		{
		case 1:	printf("����������ʽ�Ľ�����\n");
			scanf("%d", &line);
			system("cls");
			origin = (double*)calloc(line * line, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("��������ʽ��Ԫ�أ�\n");
			for (i = 0; i < line * line; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("�����������ʽΪ:\n");
			show_matrix(line, line, origin);
			value = det_value(line, origin);
			printf("����������ʽΪ:\n");
			show_matrix(line, line, origin);
			printf("������ʽ��ֵΪ:%lf\n", value);
			free(origin);
			system("pause");
			break;
		case 2:	printf("������������������������\n");
			scanf("%d %d", &column, &line);
			system("cls");
			origin = (double*)calloc(line * column, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("��������Ԫ�أ�\n");
			for (i = 0; i < line * column; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("������ľ���Ϊ:\n");
			show_matrix(column, line, origin);
			simplify_matrix(column, line, origin, origin);
			printf("�����ľ���Ϊ:\n");
			show_matrix(column, line, origin);
			free(origin);
			system("pause");
			break;
		case 3: printf("�����������������������\n");
			scanf("%d %d", &column, &line);
			system("cls");
			origin = (double*)calloc(line * column, sizeof(double));
			if (origin == NULL)
				exit(false);
			printf("�������Ԫ�أ�\n");
			for (i = 0; i < line * column; i++)
			{
				scanf("%lf", &origin[i]);
			}
			getchar();
			system("cls");
			printf("������Ļ�Ϊ:\n");
			show_matrix(column, line, origin);
			value = Schmidt_orthogonalization(column, line, origin);
			if (value == true)
			{
				printf("ʩ������������ı�׼������Ϊ:\n");
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