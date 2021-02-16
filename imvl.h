/*
        Библиотека IMVL (Матрично - векторные операции)
        Автор      Куликов Игорь
        Дата       15.12.2020
        Файл       imvl.h
        Версия     1.3
        Основные данные:	
                        - размер матрицы (IMVL_dim)
                        - указатель выделенной памяти (IMVL_memory)
                        - указатель на первый свободный элемент памяти (IMVL_start_ptr)
                        - указатель для сброса матрицы (IMVL_next_start)
                        - значения элементов матрицы (IMVL_val)
                        - инвертированная главная диагональ (IMVL_di)
                        - индексы столбцов элементов (IMVL_col_ind)
                        - указатели на начало каждой строки (IMVL_row_ptr)
                        - вектор решения (IMVL_Solve)
                        - вектор правой части (IMVL_right)
        Вспомогательные данные:
                        - указатель на конец списка (IMVL_end_list)
                        - список (IMVL_list)
        Методы:	
                        - создание матрицы (IMVL_CreateMatrix)
                        - создание вектора (IMVL_CreateVector)
                        - уничтожение матрицы (IMVL_DestroyMatrix)
                        - сброс значений матрицы (IMVL_Reset)
                        - создание элемента матрицы (IMVL_CreateElement)
                        - построение портрета матрицы (IMVL_CreatePortrait)
                        - добавление значения в матрицу (IMVL_AddElement)
                        - добавление значения в вектор (IMVL_AddElementVector)
                        - умножение матрицы на вектор (IMVL_MultMatrixVector)
                        - умножение диагональной матрицы на вектор (IMVL_MultDiagMatrixVector)
                        - сложение векторов (IMVL_SummVectorMultConstant)
                        - копирование векторов (IMVL_CopyVector)
                        - скалярное произведение (IMVL_ScalarMultiply)
                        - установка начального приближения (IMVL_SetStartVector)
                        - решение СЛАУ методом сопряжённых градиентов (IMVL_ConjGradientMethod)
                        - сохранение матрицы в CRS формате (IMVL_Save)
*/

/* Подключение необходимых заголовочных файлов */
#include <memory.h>
#include <math.h>

/* Макроопределения типов и операций над ними */
#define byte            char                    // байт данных
#define integer         long int                // целый тип
#define real            double                  // вещественный тип
#define sizeint         sizeof(integer)         // размер целого
#define sizereal        sizeof(real)            // размер вещественного

/* Определения ошибок и сообщений */
#define IMVL_OK                 0       // всё в порядке
#define IMVL_OUT_OF_SPACE       1000    // недостаточно памяти 
#define IMVL_ELEMENT_CREATED    1001    // элемент существует
#define IMVL_EMPTY_ERROR        1002    // непонятная ошибка
#define IMVL_ELEMENT_NOT_FOUND  1003    // элемент не найден
#define IMVL_NULL_DIAGONAL      1004    // ноль на диагонали
#define IMVL_RUN_TO_MAX_ITER    1005    // выход по итерации
#define IMVL_RUN_TO_EPSILON     1006    // выход по невязке
#define IMVL_FILE_NOT_OPEN      1007    // файл не удалось открыть

/* Определение константы точности*/
real    IMVL_EPS;

/* Определение данных */
integer         IMVL_dim;       // размер матрицы 
byte*           IMVL_memory;    // память
byte*           IMVL_start_ptr; // начало свободной памяти
byte*           IMVL_next_start;// начало памяти для сброса
real*           IMVL_val;       // значения элементов
real*           IMVL_di;        // инвертированная диагональ
integer*        IMVL_col_ind;   // номера столбцов
integer*        IMVL_row_ptr;   // начала строк
real*	        IMVL_Solve;     // решение СЛАУ
real*	        IMVL_right;     // вектор правой части
                                        
integer         IMVL_end_list;  // конец списка
integer*        IMVL_list;      // список

/*********************************************************** 
                   Определение методов 
 ***********************************************************/

/***** Создание матрицы *****/
integer IMVL_CreateMatrix( integer Dimension,   // размер матрицы
                           integer DumpMemory,  // кусок памяти
                           real Epsilon         // точность
                         )
{
        /* Локальные перменные */
        integer i;

        /* Определяем точность */
        IMVL_EPS = Epsilon;

        /* Создаём память для программы */
        #ifdef __cplusplus
             IMVL_memory = new byte[DumpMemory];
        #else
             IMVL_memory = (byte *)malloc(DumpMemory);
        #endif
	
        /* Если памяти мало говорим об этом */
        if(IMVL_memory == NULL) return IMVL_OUT_OF_SPACE;

        /* Устанавливаем размер матрицы */
        IMVL_dim = Dimension;

        /* Подготавливаем список */
        IMVL_list = (integer *)IMVL_memory;
        IMVL_end_list = 2*IMVL_dim;
        for(i=0 ; i<IMVL_end_list ; i++) IMVL_list[i] = -1;

        /* Возвращаем успешное завершение */
        return IMVL_OK;
}

/***** Создание вектора *****/
real* IMVL_CreateVector(integer Dim)
{
    real *tptr = (real *)IMVL_start_ptr;
    IMVL_start_ptr += sizereal * Dim;
    return tptr;
}

/***** Создание элемента матрицы *****/
integer IMVL_CreateElement( integer i,  // номер строки
                            integer k   // номер столбца
                          )
{
        /* Локальные переменные */
        integer NextElem;

        /* Начало нужной строки */
        integer CurElem = 2*i;
	
        /* Пока не найдём будем искать */
        while(1)
        {
                /* Если последний элемент в списке и
                   он сам не последний*/
                if( IMVL_list[CurElem+1] == -1 &&
                    IMVL_list[CurElem] != k)
                {
                    /* Создаём элемент в списке */
                    IMVL_list[IMVL_end_list] = k;	
                    IMVL_list[IMVL_end_list+1] = -1;

                    /* Настраиваем указатель на созданный элемент */
                    IMVL_list[CurElem+1] = IMVL_end_list;

                    /* Меняем размер действительного списка */
                    IMVL_end_list += 2;

                    /* Элемент добавлен */
                    return IMVL_OK;
                }
                else
                {
                    /* Берём следующий элемент */
                    NextElem = IMVL_list[CurElem+1];

                    /* Если элемент между двумя соседями */
                    if( IMVL_list[CurElem] < k && 
                        IMVL_list[NextElem] > k)
                    {
                          /* Создаём элемент */
                          IMVL_list[IMVL_end_list] = k;
                          IMVL_list[IMVL_end_list+1] = NextElem;

                          /* Настраиваем указатель на новый элемент */
                          IMVL_list[CurElem+1] = IMVL_end_list;

                          /* Меняем размер действительного списка */
                          IMVL_end_list += 2;

                          /* Элемент добавлен */
                          return IMVL_OK;				
                    }
                   else
                   {
                          /* Если элемент уже существует */
                          if(IMVL_list[CurElem] == k)
                               return IMVL_ELEMENT_CREATED;

                          /* Иначе переходим к следующему элменту */
                          else CurElem = NextElem;
                   }
                }
        }
        return IMVL_EMPTY_ERROR;
}

/***** Создание портрета матрицы *****/
integer IMVL_CreatePortrait()
{
        /* Локальные данные */
        integer CurElem, i;

        /* Номер текущего элемента */
        integer NumElem = 0;

        /* Установим начало размещения портрета */
        IMVL_start_ptr = (byte *)(IMVL_memory+IMVL_end_list*sizeint);

        /* Создадим место под массив указателей на начало строк 
           и индексов элементов*/
        IMVL_row_ptr = (integer *)IMVL_start_ptr;
        IMVL_col_ind = (integer *)(IMVL_start_ptr+(IMVL_dim+1)*sizeint);

        /* Первая строка идёт с начала */
        IMVL_row_ptr[0] = NumElem;

        /* Формирование из списка портрета */
        for(i=0 ; i<IMVL_dim ; i++)
        {
             /* Указатель на список строки */
             CurElem = IMVL_list[2*i+1];

             /* Пока не конец списка */
             while(CurElem!=-1)
             {
                 /* Вносим элементы строки */
                 IMVL_col_ind[NumElem] = IMVL_list[CurElem];

                 /* Увеличиваем количество элементов */
                 NumElem++;

                 /* Переходим на следующий элемент */
                 CurElem = IMVL_list[CurElem+1];
             }

             /* Записываем указатель на следующую строку */
             IMVL_row_ptr[i+1] = NumElem;
        }

        /* Переписываем в начало памяти */
        memcpy(IMVL_list,IMVL_row_ptr,(IMVL_dim+NumElem+1)*sizeint);

        /* Настраиваем указатели */
        IMVL_row_ptr = IMVL_list;
        IMVL_col_ind = IMVL_list + IMVL_dim +1;

        /* Определяем начало массива элементов матрицы и обнуляем его */
        IMVL_val = (real *)(IMVL_row_ptr + IMVL_dim + NumElem + 1);
        memset(IMVL_val,0,NumElem*sizereal);

        /* Определяем начало свободной области памяти и
           записываем туда обнулённую диагональ */
        IMVL_start_ptr = (byte *)(IMVL_val + NumElem);
        IMVL_di = (real *)(IMVL_start_ptr);
        memset(IMVL_di,0,IMVL_dim*sizereal);
        IMVL_start_ptr = (byte *)(IMVL_di + IMVL_dim);

        /* Создаём нулевой вектор правой части */
        IMVL_right = (real *)(IMVL_start_ptr);
        IMVL_start_ptr += IMVL_dim * sizereal;
        memset(IMVL_right,0,IMVL_dim * sizereal);

        /* Запомним откуда начинаются рабочие вектора */
        IMVL_next_start = IMVL_start_ptr;

        return IMVL_OK;
}

/***** Добавление элемента в матрицу *****/
integer IMVL_AddElement( integer i,     // строка
                         integer k,     // столбец
                         real value     // значение
                       )
{
        /* Локальные данные */
        integer ind, end;

        /* Проверяем на принадлежность диагонали */
        if(i==k) IMVL_di[i] += value;

        /* Цикл по строке */
        end = IMVL_row_ptr[i+1];
        for(ind=IMVL_row_ptr[i] ; ind<end ; ind++)
        {
            /* Если элемент найден */
            if(IMVL_col_ind[ind] == k)
            {
               /* Добаваляем значение */
               IMVL_val[ind] += value;
               return IMVL_OK;
            }
        }
        
        return IMVL_ELEMENT_NOT_FOUND;
}

/***** Добавление элемента в вектор правой части *****/
void IMVL_AddElementVector( integer i, real value)
{
     IMVL_right[i] += value;
}

/***** Предобуславливание матрицы *****/
integer IMVL_PreConditioner()
{
        /* Локальные переменные */
        integer i;
	
        /* Преобразование диагонали */
        for(i=0 ; i<IMVL_dim ; i++)
        {
            if(fabs(IMVL_di[i])<IMVL_EPS) 
                return IMVL_NULL_DIAGONAL;
            else
                IMVL_di[i] = 1.0/IMVL_di[i];
        }

        return IMVL_OK;
}

/***** Уничтожение матрицы и векторов *****/
void IMVL_Destroy()
{
        #ifdef __cplusplus
           delete IMVL_memory;
        #else
           free(IMVL_memory);
        #endif
}

/***** Умножение матрицы на вектор *****/
void IMVL_MultMatrixVector( real* InVector,  // входящий вектор
                            real* OutVector  // полученный вектор
                          )
{
        /* Переопределение указателей и размера матрицы */
        integer *col_ind = IMVL_col_ind,
                *row_ptr = IMVL_row_ptr;
        integer  dim = IMVL_dim;
        real    *val = IMVL_val;

        /* Локальные данные */
        integer i, k, end;
        real    valelem;

        /* Умножение матрицы на вектор */
        for(i=0 ; i<dim ; i++)
        {
            valelem = 0.0;
            end = row_ptr[i+1];
            /* умножение только не нулевых элементов */
            for(k=row_ptr[i] ; k<end ; k++)
            {
               valelem += val[k] * InVector[col_ind[k]];
            }
            OutVector[i] = valelem;
        }
}

/***** Умножение вектора на обратную диагональ матрицы *****/
void IMVL_MultDiagMatrixVector( real* InVector,  // входящий вектор
                                real* OutVector  // полученный вектор
                              )
{
        /* Переопределение размера и указателей матрицы */
        integer dim = IMVL_dim;
        real    *di = IMVL_di;

        /* Локальные перемнные */
        integer i;

        /* Умножение диагонали */
        for(i=0 ; i<dim ; i++)
            OutVector[i] = di[i] * InVector[i];
}

/***** Сложение векторов с домножением X=Y+aZ *****/
void IMVL_SummVectorMultConstant( real* y, real* z,
                                  real al, real* x)
{
        /* Переопределение размера матрицы */
        integer dim = IMVL_dim;

        /* Локальные данные */
        integer i;

        /* Вычислительный цикл */
        for(i=0 ; i<dim ; i++)
             x[i] = y[i] + al * z[i];
}

/***** Скалярное произведение *****/
real IMVL_ScalarMultiply( real* x, real* y)
{
        /* Переопределение размера матрицы */
        integer dim = IMVL_dim;

        /* Локальные данные */
        integer i;
        real dot = 0.0;

        /* Вычислительный цикл */
        for(i=0 ; i<dim ; i++) dot += x[i] * y[i];

        /* Скалярное проивздение */
        return dot;
}

/***** Установка начального приближения *****/
void IMVL_SetStartVector(real *x)
{
       /* установка нулевого начального приближения */
       memset(x,0,IMVL_dim*sizereal); 
}

/***** Копирование векторов *****/
void IMVL_CopyVector(real *dest, real *src)
{
       memcpy(dest,src,IMVL_dim*sizereal);
}

/***** Метод сопряжённых градиентов *****/
integer IMVL_ConjGradientMethod(integer MaxIter)
{
        /* Переопределение размера */
        integer dim = IMVL_dim;

        /* Локальные данные */
        integer NumIter = 0;
        real NormRight, NormR, Residual, 
             alpha, alphach, alphazn,
             beta, betach, betazn;

        /* Необходимые вектора */
        real *x = IMVL_CreateVector(IMVL_dim),
             *q = IMVL_CreateVector(IMVL_dim),
             *r = IMVL_CreateVector(IMVL_dim),
             *f = IMVL_right,
             *z = IMVL_CreateVector(IMVL_dim);

        /* Реализация метода */
        IMVL_PreConditioner();
        NormRight = sqrt(IMVL_ScalarMultiply(f,f));     // |f|
        IMVL_SetStartVector(x);                         // x = x0
        IMVL_MultMatrixVector(x,q);                     // q = Ax
        IMVL_SummVectorMultConstant(f,q,-1,r);          // r = f-q
        IMVL_MultDiagMatrixVector(r,q);                 // q = D^(-1)r
        IMVL_CopyVector(z,q);                           // z = q
        NormR = sqrt(IMVL_ScalarMultiply(r,r));         // |r|
        Residual = NormR/NormRight;                     // |r|/|f|
        alphach = IMVL_ScalarMultiply(q,r);             // alch = (q,r)
	
        /* Пока не конец по итерациям или по невязке */
        while(Residual>IMVL_EPS && NumIter<MaxIter)		
        {
           IMVL_MultMatrixVector(z,q);                  // q = Az
           alphazn = IMVL_ScalarMultiply(q,z);          // alzn = (q,z)	
           alpha = alphach/alphazn;                     // al = alch/alzn
           IMVL_SummVectorMultConstant(x,z,alpha,x);    // x = x + al * z
           IMVL_SummVectorMultConstant(r,q,-alpha,r);   // r = r - al * q		
           IMVL_MultDiagMatrixVector(r,q);              // q = D^(-1)r
           betach = IMVL_ScalarMultiply(q,r);           // btch = (q,r)
           betazn = alphach;                            // btzn = alch
           alphach = betach;                            // alch = betach
           beta = betach/betazn;                        // bt = btch/btzn
           IMVL_SummVectorMultConstant(q,z,beta,z);     // z = q + bt * z
           NormR = sqrt(IMVL_ScalarMultiply(r,r));      // |r|
           Residual = NormR/NormRight;                  // |r|/|f|
           NumIter++;
           printf("Number iteration %d Residual %le\n",NumIter,Residual);
        }

        /* Указатель ответа настраиваем на решение */
        IMVL_Solve = x;

        /* Если выход по невязке */
        if(NumIter<MaxIter) return IMVL_RUN_TO_EPSILON;
        else return IMVL_RUN_TO_MAX_ITER;
}

/***** Сброс значений матрицы и вектора правой части *****/
void IMVL_Reset()
{
        /* Сброс матрицы */
        memset(IMVL_val,0,IMVL_row_ptr[IMVL_dim]*sizereal);
	
        /* Сброс вектора правой части */
        memset(IMVL_right,0,IMVL_dim*sizereal);

        /* Сброс предобуславливающей диагонали */
        memset(IMVL_di,0,IMVL_dim*sizereal);

        /* Переопределение указателя */
        IMVL_start_ptr = IMVL_next_start;
}

/***** Сохранение матрицы и вектора правой части *****/
integer IMVL_Save()
{
        FILE* fout;
	integer i;
        /* Сохранение параметров СЛАУ */
	fout = fopen("param","w");
	if(fout)
	{
		fprintf(fout,"%d\n",IMVL_dim);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* Сохранение указателей на начало строк */
	fout = fopen("row","w");
	if(fout)
	{
		for(i=0 ; i<=IMVL_dim ; i++)
			fprintf(fout,"%d\n",IMVL_row_ptr[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* Сохранение правой части */
	fout = fopen("f","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_dim ; i++)
			fprintf(fout,"%g\n",IMVL_right[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* Сохранение диагонали */
	fout = fopen("di","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_dim ; i++)
			fprintf(fout,"%g\n",IMVL_di[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* Сохранение значений матрицы */
	fout = fopen("val","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_row_ptr[IMVL_dim] ; i++)
			fprintf(fout,"%g\n",IMVL_val[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* Сохранение столбцов элементов */
	fout = fopen("col","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_row_ptr[IMVL_dim] ; i++)
			fprintf(fout,"%d\n",IMVL_col_ind[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	return IMVL_OK;
}
