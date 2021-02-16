/*
        ���������� IMVL (�������� - ��������� ��������)
        �����      ������� �����
        ����       15.12.2020
        ����       imvl.h
        ������     1.3
        �������� ������:	
                        - ������ ������� (IMVL_dim)
                        - ��������� ���������� ������ (IMVL_memory)
                        - ��������� �� ������ ��������� ������� ������ (IMVL_start_ptr)
                        - ��������� ��� ������ ������� (IMVL_next_start)
                        - �������� ��������� ������� (IMVL_val)
                        - ��������������� ������� ��������� (IMVL_di)
                        - ������� �������� ��������� (IMVL_col_ind)
                        - ��������� �� ������ ������ ������ (IMVL_row_ptr)
                        - ������ ������� (IMVL_Solve)
                        - ������ ������ ����� (IMVL_right)
        ��������������� ������:
                        - ��������� �� ����� ������ (IMVL_end_list)
                        - ������ (IMVL_list)
        ������:	
                        - �������� ������� (IMVL_CreateMatrix)
                        - �������� ������� (IMVL_CreateVector)
                        - ����������� ������� (IMVL_DestroyMatrix)
                        - ����� �������� ������� (IMVL_Reset)
                        - �������� �������� ������� (IMVL_CreateElement)
                        - ���������� �������� ������� (IMVL_CreatePortrait)
                        - ���������� �������� � ������� (IMVL_AddElement)
                        - ���������� �������� � ������ (IMVL_AddElementVector)
                        - ��������� ������� �� ������ (IMVL_MultMatrixVector)
                        - ��������� ������������ ������� �� ������ (IMVL_MultDiagMatrixVector)
                        - �������� �������� (IMVL_SummVectorMultConstant)
                        - ����������� �������� (IMVL_CopyVector)
                        - ��������� ������������ (IMVL_ScalarMultiply)
                        - ��������� ���������� ����������� (IMVL_SetStartVector)
                        - ������� ���� ������� ���������� ���������� (IMVL_ConjGradientMethod)
                        - ���������� ������� � CRS ������� (IMVL_Save)
*/

/* ����������� ����������� ������������ ������ */
#include <memory.h>
#include <math.h>

/* ���������������� ����� � �������� ��� ���� */
#define byte            char                    // ���� ������
#define integer         long int                // ����� ���
#define real            double                  // ������������ ���
#define sizeint         sizeof(integer)         // ������ ������
#define sizereal        sizeof(real)            // ������ �������������

/* ����������� ������ � ��������� */
#define IMVL_OK                 0       // �� � �������
#define IMVL_OUT_OF_SPACE       1000    // ������������ ������ 
#define IMVL_ELEMENT_CREATED    1001    // ������� ����������
#define IMVL_EMPTY_ERROR        1002    // ���������� ������
#define IMVL_ELEMENT_NOT_FOUND  1003    // ������� �� ������
#define IMVL_NULL_DIAGONAL      1004    // ���� �� ���������
#define IMVL_RUN_TO_MAX_ITER    1005    // ����� �� ��������
#define IMVL_RUN_TO_EPSILON     1006    // ����� �� �������
#define IMVL_FILE_NOT_OPEN      1007    // ���� �� ������� �������

/* ����������� ��������� ��������*/
real    IMVL_EPS;

/* ����������� ������ */
integer         IMVL_dim;       // ������ ������� 
byte*           IMVL_memory;    // ������
byte*           IMVL_start_ptr; // ������ ��������� ������
byte*           IMVL_next_start;// ������ ������ ��� ������
real*           IMVL_val;       // �������� ���������
real*           IMVL_di;        // ��������������� ���������
integer*        IMVL_col_ind;   // ������ ��������
integer*        IMVL_row_ptr;   // ������ �����
real*	        IMVL_Solve;     // ������� ����
real*	        IMVL_right;     // ������ ������ �����
                                        
integer         IMVL_end_list;  // ����� ������
integer*        IMVL_list;      // ������

/*********************************************************** 
                   ����������� ������� 
 ***********************************************************/

/***** �������� ������� *****/
integer IMVL_CreateMatrix( integer Dimension,   // ������ �������
                           integer DumpMemory,  // ����� ������
                           real Epsilon         // ��������
                         )
{
        /* ��������� ��������� */
        integer i;

        /* ���������� �������� */
        IMVL_EPS = Epsilon;

        /* ������ ������ ��� ��������� */
        #ifdef __cplusplus
             IMVL_memory = new byte[DumpMemory];
        #else
             IMVL_memory = (byte *)malloc(DumpMemory);
        #endif
	
        /* ���� ������ ���� ������� �� ���� */
        if(IMVL_memory == NULL) return IMVL_OUT_OF_SPACE;

        /* ������������� ������ ������� */
        IMVL_dim = Dimension;

        /* �������������� ������ */
        IMVL_list = (integer *)IMVL_memory;
        IMVL_end_list = 2*IMVL_dim;
        for(i=0 ; i<IMVL_end_list ; i++) IMVL_list[i] = -1;

        /* ���������� �������� ���������� */
        return IMVL_OK;
}

/***** �������� ������� *****/
real* IMVL_CreateVector(integer Dim)
{
    real *tptr = (real *)IMVL_start_ptr;
    IMVL_start_ptr += sizereal * Dim;
    return tptr;
}

/***** �������� �������� ������� *****/
integer IMVL_CreateElement( integer i,  // ����� ������
                            integer k   // ����� �������
                          )
{
        /* ��������� ���������� */
        integer NextElem;

        /* ������ ������ ������ */
        integer CurElem = 2*i;
	
        /* ���� �� ����� ����� ������ */
        while(1)
        {
                /* ���� ��������� ������� � ������ �
                   �� ��� �� ���������*/
                if( IMVL_list[CurElem+1] == -1 &&
                    IMVL_list[CurElem] != k)
                {
                    /* ������ ������� � ������ */
                    IMVL_list[IMVL_end_list] = k;	
                    IMVL_list[IMVL_end_list+1] = -1;

                    /* ����������� ��������� �� ��������� ������� */
                    IMVL_list[CurElem+1] = IMVL_end_list;

                    /* ������ ������ ��������������� ������ */
                    IMVL_end_list += 2;

                    /* ������� �������� */
                    return IMVL_OK;
                }
                else
                {
                    /* ���� ��������� ������� */
                    NextElem = IMVL_list[CurElem+1];

                    /* ���� ������� ����� ����� �������� */
                    if( IMVL_list[CurElem] < k && 
                        IMVL_list[NextElem] > k)
                    {
                          /* ������ ������� */
                          IMVL_list[IMVL_end_list] = k;
                          IMVL_list[IMVL_end_list+1] = NextElem;

                          /* ����������� ��������� �� ����� ������� */
                          IMVL_list[CurElem+1] = IMVL_end_list;

                          /* ������ ������ ��������������� ������ */
                          IMVL_end_list += 2;

                          /* ������� �������� */
                          return IMVL_OK;				
                    }
                   else
                   {
                          /* ���� ������� ��� ���������� */
                          if(IMVL_list[CurElem] == k)
                               return IMVL_ELEMENT_CREATED;

                          /* ����� ��������� � ���������� ������� */
                          else CurElem = NextElem;
                   }
                }
        }
        return IMVL_EMPTY_ERROR;
}

/***** �������� �������� ������� *****/
integer IMVL_CreatePortrait()
{
        /* ��������� ������ */
        integer CurElem, i;

        /* ����� �������� �������� */
        integer NumElem = 0;

        /* ��������� ������ ���������� �������� */
        IMVL_start_ptr = (byte *)(IMVL_memory+IMVL_end_list*sizeint);

        /* �������� ����� ��� ������ ���������� �� ������ ����� 
           � �������� ���������*/
        IMVL_row_ptr = (integer *)IMVL_start_ptr;
        IMVL_col_ind = (integer *)(IMVL_start_ptr+(IMVL_dim+1)*sizeint);

        /* ������ ������ ��� � ������ */
        IMVL_row_ptr[0] = NumElem;

        /* ������������ �� ������ �������� */
        for(i=0 ; i<IMVL_dim ; i++)
        {
             /* ��������� �� ������ ������ */
             CurElem = IMVL_list[2*i+1];

             /* ���� �� ����� ������ */
             while(CurElem!=-1)
             {
                 /* ������ �������� ������ */
                 IMVL_col_ind[NumElem] = IMVL_list[CurElem];

                 /* ����������� ���������� ��������� */
                 NumElem++;

                 /* ��������� �� ��������� ������� */
                 CurElem = IMVL_list[CurElem+1];
             }

             /* ���������� ��������� �� ��������� ������ */
             IMVL_row_ptr[i+1] = NumElem;
        }

        /* ������������ � ������ ������ */
        memcpy(IMVL_list,IMVL_row_ptr,(IMVL_dim+NumElem+1)*sizeint);

        /* ����������� ��������� */
        IMVL_row_ptr = IMVL_list;
        IMVL_col_ind = IMVL_list + IMVL_dim +1;

        /* ���������� ������ ������� ��������� ������� � �������� ��� */
        IMVL_val = (real *)(IMVL_row_ptr + IMVL_dim + NumElem + 1);
        memset(IMVL_val,0,NumElem*sizereal);

        /* ���������� ������ ��������� ������� ������ �
           ���������� ���� ��������� ��������� */
        IMVL_start_ptr = (byte *)(IMVL_val + NumElem);
        IMVL_di = (real *)(IMVL_start_ptr);
        memset(IMVL_di,0,IMVL_dim*sizereal);
        IMVL_start_ptr = (byte *)(IMVL_di + IMVL_dim);

        /* ������ ������� ������ ������ ����� */
        IMVL_right = (real *)(IMVL_start_ptr);
        IMVL_start_ptr += IMVL_dim * sizereal;
        memset(IMVL_right,0,IMVL_dim * sizereal);

        /* �������� ������ ���������� ������� ������� */
        IMVL_next_start = IMVL_start_ptr;

        return IMVL_OK;
}

/***** ���������� �������� � ������� *****/
integer IMVL_AddElement( integer i,     // ������
                         integer k,     // �������
                         real value     // ��������
                       )
{
        /* ��������� ������ */
        integer ind, end;

        /* ��������� �� �������������� ��������� */
        if(i==k) IMVL_di[i] += value;

        /* ���� �� ������ */
        end = IMVL_row_ptr[i+1];
        for(ind=IMVL_row_ptr[i] ; ind<end ; ind++)
        {
            /* ���� ������� ������ */
            if(IMVL_col_ind[ind] == k)
            {
               /* ���������� �������� */
               IMVL_val[ind] += value;
               return IMVL_OK;
            }
        }
        
        return IMVL_ELEMENT_NOT_FOUND;
}

/***** ���������� �������� � ������ ������ ����� *****/
void IMVL_AddElementVector( integer i, real value)
{
     IMVL_right[i] += value;
}

/***** ������������������ ������� *****/
integer IMVL_PreConditioner()
{
        /* ��������� ���������� */
        integer i;
	
        /* �������������� ��������� */
        for(i=0 ; i<IMVL_dim ; i++)
        {
            if(fabs(IMVL_di[i])<IMVL_EPS) 
                return IMVL_NULL_DIAGONAL;
            else
                IMVL_di[i] = 1.0/IMVL_di[i];
        }

        return IMVL_OK;
}

/***** ����������� ������� � �������� *****/
void IMVL_Destroy()
{
        #ifdef __cplusplus
           delete IMVL_memory;
        #else
           free(IMVL_memory);
        #endif
}

/***** ��������� ������� �� ������ *****/
void IMVL_MultMatrixVector( real* InVector,  // �������� ������
                            real* OutVector  // ���������� ������
                          )
{
        /* ��������������� ���������� � ������� ������� */
        integer *col_ind = IMVL_col_ind,
                *row_ptr = IMVL_row_ptr;
        integer  dim = IMVL_dim;
        real    *val = IMVL_val;

        /* ��������� ������ */
        integer i, k, end;
        real    valelem;

        /* ��������� ������� �� ������ */
        for(i=0 ; i<dim ; i++)
        {
            valelem = 0.0;
            end = row_ptr[i+1];
            /* ��������� ������ �� ������� ��������� */
            for(k=row_ptr[i] ; k<end ; k++)
            {
               valelem += val[k] * InVector[col_ind[k]];
            }
            OutVector[i] = valelem;
        }
}

/***** ��������� ������� �� �������� ��������� ������� *****/
void IMVL_MultDiagMatrixVector( real* InVector,  // �������� ������
                                real* OutVector  // ���������� ������
                              )
{
        /* ��������������� ������� � ���������� ������� */
        integer dim = IMVL_dim;
        real    *di = IMVL_di;

        /* ��������� ��������� */
        integer i;

        /* ��������� ��������� */
        for(i=0 ; i<dim ; i++)
            OutVector[i] = di[i] * InVector[i];
}

/***** �������� �������� � ����������� X=Y+aZ *****/
void IMVL_SummVectorMultConstant( real* y, real* z,
                                  real al, real* x)
{
        /* ��������������� ������� ������� */
        integer dim = IMVL_dim;

        /* ��������� ������ */
        integer i;

        /* �������������� ���� */
        for(i=0 ; i<dim ; i++)
             x[i] = y[i] + al * z[i];
}

/***** ��������� ������������ *****/
real IMVL_ScalarMultiply( real* x, real* y)
{
        /* ��������������� ������� ������� */
        integer dim = IMVL_dim;

        /* ��������� ������ */
        integer i;
        real dot = 0.0;

        /* �������������� ���� */
        for(i=0 ; i<dim ; i++) dot += x[i] * y[i];

        /* ��������� ����������� */
        return dot;
}

/***** ��������� ���������� ����������� *****/
void IMVL_SetStartVector(real *x)
{
       /* ��������� �������� ���������� ����������� */
       memset(x,0,IMVL_dim*sizereal); 
}

/***** ����������� �������� *****/
void IMVL_CopyVector(real *dest, real *src)
{
       memcpy(dest,src,IMVL_dim*sizereal);
}

/***** ����� ���������� ���������� *****/
integer IMVL_ConjGradientMethod(integer MaxIter)
{
        /* ��������������� ������� */
        integer dim = IMVL_dim;

        /* ��������� ������ */
        integer NumIter = 0;
        real NormRight, NormR, Residual, 
             alpha, alphach, alphazn,
             beta, betach, betazn;

        /* ����������� ������� */
        real *x = IMVL_CreateVector(IMVL_dim),
             *q = IMVL_CreateVector(IMVL_dim),
             *r = IMVL_CreateVector(IMVL_dim),
             *f = IMVL_right,
             *z = IMVL_CreateVector(IMVL_dim);

        /* ���������� ������ */
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
	
        /* ���� �� ����� �� ��������� ��� �� ������� */
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

        /* ��������� ������ ����������� �� ������� */
        IMVL_Solve = x;

        /* ���� ����� �� ������� */
        if(NumIter<MaxIter) return IMVL_RUN_TO_EPSILON;
        else return IMVL_RUN_TO_MAX_ITER;
}

/***** ����� �������� ������� � ������� ������ ����� *****/
void IMVL_Reset()
{
        /* ����� ������� */
        memset(IMVL_val,0,IMVL_row_ptr[IMVL_dim]*sizereal);
	
        /* ����� ������� ������ ����� */
        memset(IMVL_right,0,IMVL_dim*sizereal);

        /* ����� ������������������� ��������� */
        memset(IMVL_di,0,IMVL_dim*sizereal);

        /* ��������������� ��������� */
        IMVL_start_ptr = IMVL_next_start;
}

/***** ���������� ������� � ������� ������ ����� *****/
integer IMVL_Save()
{
        FILE* fout;
	integer i;
        /* ���������� ���������� ���� */
	fout = fopen("param","w");
	if(fout)
	{
		fprintf(fout,"%d\n",IMVL_dim);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* ���������� ���������� �� ������ ����� */
	fout = fopen("row","w");
	if(fout)
	{
		for(i=0 ; i<=IMVL_dim ; i++)
			fprintf(fout,"%d\n",IMVL_row_ptr[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* ���������� ������ ����� */
	fout = fopen("f","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_dim ; i++)
			fprintf(fout,"%g\n",IMVL_right[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* ���������� ��������� */
	fout = fopen("di","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_dim ; i++)
			fprintf(fout,"%g\n",IMVL_di[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* ���������� �������� ������� */
	fout = fopen("val","w");
	if(fout)
	{
		for(i=0 ; i<IMVL_row_ptr[IMVL_dim] ; i++)
			fprintf(fout,"%g\n",IMVL_val[i]);
	}
	else
		return IMVL_FILE_NOT_OPEN;			
	fclose(fout);

	/* ���������� �������� ��������� */
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
