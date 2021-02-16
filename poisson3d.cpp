/* 
   Решение уравнения Пуассона методом конечных элементов 
*/
#include <stdio.h>
#include "imvl.h"

/* Объявление констант */
#define PI 3.141592653589793

/* Параметры задачи */
real     sx, sy,sz;         // размер области
integer  nx, ny, nz;        // количество конечных элементов
real     h, hx, hy, hz;     // шаги по осям

/* Функция индекса */
integer index(integer i, integer k, integer l)
{
	return i*(nz+1)*(ny+1) + k*(nz+1) + l;
}

/* Правая часть */
real right(real x, real y, real z)
{
	real r = sqrt(x*x + y*y + z*z);
	if(r <= 1.0)
		return 4.0 * PI * (2.0*r*r*r - 3.0*r*r + 1.0);
	else
		return 0.0;
}

/* Точное решение */
real exact(real x, real y, real z)
{
	real r;
        return x*x + y*y + z*z;
	r = sqrt(x*x + y*y +z*z);
	if(r <= 1.0) 
		return  4.0*PI*r*r*r*r*r/15.0 - 3.0*PI*r*r*r*r/5.0 + 2.0*PI*r*r/3.0 - 3.0*PI/5.0;
	else
		return -4.0*PI/r/15.0;
}

/* Основная функция */
int main()
{
        /* Локальные переменные */
        integer i, k, l, nppp, nppm, npmp, npmm, nmpp, nmpm, nmmp, nmmm;
        real x, y, z, hx2, hy2, hz2, tmp, tmp2, fppp, fppm, fpmp, fpmm, fmpp, fmpm, fmmp, fmmm;
        real B[8][8], C[8][8], F[8];
        FILE* FSolve;

        /* Чтение параметров задачи */
        sx = 4.0;
        sy = 4.0;
        sz = 4.0;

        nx = 64;
        ny = 64;
        nz = 64;

        hx = sx/(real)nx;
        hy = sy/(real)ny;
        hz = sz/(real)nz;
        h  = hx;

        hx2 = hx*hx; 
        hy2 = hy*hy;
        hz2 = hz*hz;

        /* Создание матрицы */
        IMVL_CreateMatrix((nx+1)*(ny+1)*(nz+1),128*1024*1024,1.0e-40);

        /* Построение портрета матрицы */
        for(i=0 ; i<nx ; i++)
          for(k=0 ; k<ny ; k++)
            for(l=0 ; l<nz ; l++)
            {
               /* Нахождение опорной точки */
               nppp = index(i+1,k+1,l+1);
               nppm = index(i+1,k+1,l);
               npmp = index(i+1,k,l+1);
               npmm = index(i+1,k,l);
               nmpp = index(i,k+1,l+1);
               nmpm = index(i,k+1,l);
               nmmp = index(i,k,l+1);
               nmmm = index(i,k,l);

               /* Создание элементов матрицы */
               IMVL_CreateElement(nppp,nppp);   IMVL_CreateElement(npmp,nppp);   IMVL_CreateElement(nmpp,nppp);   IMVL_CreateElement(nmmp,nppp);
               IMVL_CreateElement(nppp,nppm);   IMVL_CreateElement(npmp,nppm);   IMVL_CreateElement(nmpp,nppm);   IMVL_CreateElement(nmmp,nppm);
               IMVL_CreateElement(nppp,npmp);   IMVL_CreateElement(npmp,npmp);   IMVL_CreateElement(nmpp,npmp);   IMVL_CreateElement(nmmp,npmp);
               IMVL_CreateElement(nppp,npmm);   IMVL_CreateElement(npmp,npmm);   IMVL_CreateElement(nmpp,npmm);   IMVL_CreateElement(nmmp,npmm);
               IMVL_CreateElement(nppp,nmpp);   IMVL_CreateElement(npmp,nmpp);   IMVL_CreateElement(nmpp,nmpp);   IMVL_CreateElement(nmmp,nmpp);
               IMVL_CreateElement(nppp,nmpm);   IMVL_CreateElement(npmp,nmpm);   IMVL_CreateElement(nmpp,nmpm);   IMVL_CreateElement(nmmp,nmpm);
               IMVL_CreateElement(nppp,nmmp);   IMVL_CreateElement(npmp,nmmp);   IMVL_CreateElement(nmpp,nmmp);   IMVL_CreateElement(nmmp,nmmp);
               IMVL_CreateElement(nppp,nmmm);   IMVL_CreateElement(npmp,nmmm);   IMVL_CreateElement(nmpp,nmmm);   IMVL_CreateElement(nmmp,nmmm);
                                                                                                                                                
               IMVL_CreateElement(nppm,nppp);   IMVL_CreateElement(npmm,nppp);   IMVL_CreateElement(nmpm,nppp);   IMVL_CreateElement(nmmm,nppp);
               IMVL_CreateElement(nppm,nppm);   IMVL_CreateElement(npmm,nppm);   IMVL_CreateElement(nmpm,nppm);   IMVL_CreateElement(nmmm,nppm);
               IMVL_CreateElement(nppm,npmp);   IMVL_CreateElement(npmm,npmp);   IMVL_CreateElement(nmpm,npmp);   IMVL_CreateElement(nmmm,npmp);
               IMVL_CreateElement(nppm,npmm);   IMVL_CreateElement(npmm,npmm);   IMVL_CreateElement(nmpm,npmm);   IMVL_CreateElement(nmmm,npmm);
               IMVL_CreateElement(nppm,nmpp);   IMVL_CreateElement(npmm,nmpp);   IMVL_CreateElement(nmpm,nmpp);   IMVL_CreateElement(nmmm,nmpp);
               IMVL_CreateElement(nppm,nmpm);   IMVL_CreateElement(npmm,nmpm);   IMVL_CreateElement(nmpm,nmpm);   IMVL_CreateElement(nmmm,nmpm);
               IMVL_CreateElement(nppm,nmmp);   IMVL_CreateElement(npmm,nmmp);   IMVL_CreateElement(nmpm,nmmp);   IMVL_CreateElement(nmmm,nmmp);
               IMVL_CreateElement(nppm,nmmm);   IMVL_CreateElement(npmm,nmmm);   IMVL_CreateElement(nmpm,nmmm);   IMVL_CreateElement(nmmm,nmmm);
          }

        /* Достроим портрет */
        IMVL_CreatePortrait();

        /* Формирование матрицы жёсткости */
        tmp = h/12;

        B[0][0] =-tmp * 4.0;  B[2][0] = tmp * 0.0;   B[4][0] = tmp * 0.0;   B[6][0] = tmp * 1.0;   
        B[0][1] = tmp * 0.0;  B[2][1] = tmp * 1.0;   B[4][1] = tmp * 1.0;   B[6][1] = tmp * 1.0;   
        B[0][2] = tmp * 0.0;  B[2][2] =-tmp * 4.0;   B[4][2] = tmp * 1.0;   B[6][2] = tmp * 0.0;   
        B[0][3] = tmp * 1.0;  B[2][3] = tmp * 0.0;   B[4][3] = tmp * 1.0;   B[6][3] = tmp * 1.0;   
        B[0][4] = tmp * 0.0;  B[2][4] = tmp * 1.0;   B[4][4] =-tmp * 4.0;   B[6][4] = tmp * 0.0;   
        B[0][5] = tmp * 1.0;  B[2][5] = tmp * 1.0;   B[4][5] = tmp * 0.0;   B[6][5] = tmp * 1.0;   
        B[0][6] = tmp * 1.0;  B[2][6] = tmp * 0.0;   B[4][6] = tmp * 0.0;   B[6][6] =-tmp * 4.0;   
        B[0][7] = tmp * 1.0;  B[2][7] = tmp * 1.0;   B[4][7] = tmp * 1.0;   B[6][7] = tmp * 0.0;   
                                                                                                   
        B[1][0] = tmp * 0.0;  B[3][0] = tmp * 1.0;   B[5][0] = tmp * 1.0;   B[7][0] = tmp * 1.0;   
        B[1][1] =-tmp * 4.0;  B[3][1] = tmp * 0.0;   B[5][1] = tmp * 0.0;   B[7][1] = tmp * 1.0;   
        B[1][2] = tmp * 1.0;  B[3][2] = tmp * 0.0;   B[5][2] = tmp * 1.0;   B[7][2] = tmp * 1.0;   
        B[1][3] = tmp * 0.0;  B[3][3] =-tmp * 4.0;   B[5][3] = tmp * 1.0;   B[7][3] = tmp * 0.0;   
        B[1][4] = tmp * 1.0;  B[3][4] = tmp * 1.0;   B[5][4] = tmp * 0.0;   B[7][4] = tmp * 1.0;   
        B[1][5] = tmp * 0.0;  B[3][5] = tmp * 1.0;   B[5][5] =-tmp * 4.0;   B[7][5] = tmp * 0.0;   
        B[1][6] = tmp * 1.0;  B[3][6] = tmp * 1.0;   B[5][6] = tmp * 1.0;   B[7][6] = tmp * 0.0;   
        B[1][7] = tmp * 1.0;  B[3][7] = tmp * 0.0;   B[5][7] = tmp * 0.0;   B[7][7] =-tmp * 4.0;   

        /* Формирование матрицы масс */
        tmp = (h*h*h)/216.0;

        C[0][0] = tmp * 8.0;  C[2][0] = tmp * 4.0;   C[4][0] = tmp * 4.0;   C[6][0] = tmp * 2.0;   
        C[0][1] = tmp * 4.0;  C[2][1] = tmp * 2.0;   C[4][1] = tmp * 2.0;   C[6][1] = tmp * 1.0;   
        C[0][2] = tmp * 4.0;  C[2][2] = tmp * 8.0;   C[4][2] = tmp * 2.0;   C[6][2] = tmp * 4.0;   
        C[0][3] = tmp * 2.0;  C[2][3] = tmp * 4.0;   C[4][3] = tmp * 1.0;   C[6][3] = tmp * 2.0;   
        C[0][4] = tmp * 4.0;  C[2][4] = tmp * 2.0;   C[4][4] = tmp * 8.0;   C[6][4] = tmp * 4.0;   
        C[0][5] = tmp * 2.0;  C[2][5] = tmp * 1.0;   C[4][5] = tmp * 4.0;   C[6][5] = tmp * 2.0;   
        C[0][6] = tmp * 2.0;  C[2][6] = tmp * 4.0;   C[4][6] = tmp * 4.0;   C[6][6] = tmp * 8.0;   
        C[0][7] = tmp * 1.0;  C[2][7] = tmp * 2.0;   C[4][7] = tmp * 2.0;   C[6][7] = tmp * 4.0;   
                                                                                                   
        C[1][0] = tmp * 4.0;  C[3][0] = tmp * 2.0;   C[5][0] = tmp * 2.0;   C[7][0] = tmp * 1.0;   
        C[1][1] = tmp * 8.0;  C[3][1] = tmp * 4.0;   C[5][1] = tmp * 4.0;   C[7][1] = tmp * 2.0;   
        C[1][2] = tmp * 2.0;  C[3][2] = tmp * 4.0;   C[5][2] = tmp * 1.0;   C[7][2] = tmp * 2.0;   
        C[1][3] = tmp * 4.0;  C[3][3] = tmp * 8.0;   C[5][3] = tmp * 2.0;   C[7][3] = tmp * 4.0;   
        C[1][4] = tmp * 2.0;  C[3][4] = tmp * 1.0;   C[5][4] = tmp * 4.0;   C[7][4] = tmp * 2.0;   
        C[1][5] = tmp * 4.0;  C[3][5] = tmp * 2.0;   C[5][5] = tmp * 8.0;   C[7][5] = tmp * 4.0;   
        C[1][6] = tmp * 1.0;  C[3][6] = tmp * 2.0;   C[5][6] = tmp * 2.0;   C[7][6] = tmp * 4.0;   
        C[1][7] = tmp * 2.0;  C[3][7] = tmp * 4.0;   C[5][7] = tmp * 4.0;   C[7][7] = tmp * 8.0;   

        /* Запуск конечноэлементной постановки */
        for(i=0 ; i<nx ; i++)
          for(k=0 ; k<ny ; k++)
            for(l=0 ; l<nz ; l++)
            {
               /* Нахождение опорной точки */
               nppp = index(i+1,k+1,l+1);
               nppm = index(i+1,k+1,l);
               npmp = index(i+1,k,l+1);
               npmm = index(i+1,k,l);
               nmpp = index(i,k+1,l+1);
               nmpm = index(i,k+1,l);
               nmmp = index(i,k,l+1);
               nmmm = index(i,k,l);

               /* Вычисление координат */
               x = i*hx - sx/2.0;
               y = k*hy - sy/2.0;
               z = l*hz - sz/2.0;

               /* Формирование матрицы правой части */
               fppp = right(x+hx, y+hy, z+hz);
               fppm = right(x+hx, y+hy, z);
               fpmp = right(x+hx, y, z+hz);
               fpmm = right(x+hx, y, z);
               fmpp = right(x, y+hy, z+hz);
               fmpm = right(x, y+hy, z);
               fmmp = right(x, y, z+hz);
               fmmm = right(x, y, z);

               F[0] = C[0][0]*fppp + C[0][1]*fppm + C[0][2]*fpmp + C[0][3]*fpmm + C[0][4]*fmpp + C[0][5]*fmpm + C[0][6]*fmmp + C[0][7]*fmmm;
               F[1] = C[1][0]*fppp + C[1][1]*fppm + C[1][2]*fpmp + C[1][3]*fpmm + C[1][4]*fmpp + C[1][5]*fmpm + C[1][6]*fmmp + C[1][7]*fmmm;
               F[2] = C[2][0]*fppp + C[2][1]*fppm + C[2][2]*fpmp + C[2][3]*fpmm + C[2][4]*fmpp + C[2][5]*fmpm + C[2][6]*fmmp + C[2][7]*fmmm;
               F[3] = C[3][0]*fppp + C[3][1]*fppm + C[3][2]*fpmp + C[3][3]*fpmm + C[3][4]*fmpp + C[3][5]*fmpm + C[3][6]*fmmp + C[3][7]*fmmm;
               F[4] = C[4][0]*fppp + C[4][1]*fppm + C[4][2]*fpmp + C[4][3]*fpmm + C[4][4]*fmpp + C[4][5]*fmpm + C[4][6]*fmmp + C[4][7]*fmmm;
               F[5] = C[5][0]*fppp + C[5][1]*fppm + C[5][2]*fpmp + C[5][3]*fpmm + C[5][4]*fmpp + C[5][5]*fmpm + C[5][6]*fmmp + C[5][7]*fmmm;
               F[6] = C[6][0]*fppp + C[6][1]*fppm + C[6][2]*fpmp + C[6][3]*fpmm + C[6][4]*fmpp + C[6][5]*fmpm + C[6][6]*fmmp + C[6][7]*fmmm;
               F[7] = C[7][0]*fppp + C[7][1]*fppm + C[7][2]*fpmp + C[7][3]*fpmm + C[7][4]*fmpp + C[7][5]*fmpm + C[7][6]*fmmp + C[7][7]*fmmm;

               /* Добавление элемента в вектор правой части */
               IMVL_AddElementVector(nppp,F[0]);
               IMVL_AddElementVector(nppm,F[1]);
               IMVL_AddElementVector(npmp,F[2]);
               IMVL_AddElementVector(npmm,F[3]);
               IMVL_AddElementVector(nmpp,F[4]);
               IMVL_AddElementVector(nmpm,F[5]);
               IMVL_AddElementVector(nmmp,F[6]);
               IMVL_AddElementVector(nmmm,F[7]);
                                       
               /* Добавление элементов в матрицу */
               IMVL_AddElement(nppp,nppp,B[0][0]);   IMVL_AddElement(npmp,nppp,B[2][0]);   IMVL_AddElement(nmpp,nppp,B[4][0]);   IMVL_AddElement(nmmp,nppp,B[6][0]);
               IMVL_AddElement(nppp,nppm,B[0][1]);   IMVL_AddElement(npmp,nppm,B[2][1]);   IMVL_AddElement(nmpp,nppm,B[4][1]);   IMVL_AddElement(nmmp,nppm,B[6][1]);
               IMVL_AddElement(nppp,npmp,B[0][2]);   IMVL_AddElement(npmp,npmp,B[2][2]);   IMVL_AddElement(nmpp,npmp,B[4][2]);   IMVL_AddElement(nmmp,npmp,B[6][2]);
               IMVL_AddElement(nppp,npmm,B[0][3]);   IMVL_AddElement(npmp,npmm,B[2][3]);   IMVL_AddElement(nmpp,npmm,B[4][3]);   IMVL_AddElement(nmmp,npmm,B[6][3]);
               IMVL_AddElement(nppp,nmpp,B[0][4]);   IMVL_AddElement(npmp,nmpp,B[2][4]);   IMVL_AddElement(nmpp,nmpp,B[4][4]);   IMVL_AddElement(nmmp,nmpp,B[6][4]);
               IMVL_AddElement(nppp,nmpm,B[0][5]);   IMVL_AddElement(npmp,nmpm,B[2][5]);   IMVL_AddElement(nmpp,nmpm,B[4][5]);   IMVL_AddElement(nmmp,nmpm,B[6][5]);
               IMVL_AddElement(nppp,nmmp,B[0][6]);   IMVL_AddElement(npmp,nmmp,B[2][6]);   IMVL_AddElement(nmpp,nmmp,B[4][6]);   IMVL_AddElement(nmmp,nmmp,B[6][6]);
               IMVL_AddElement(nppp,nmmm,B[0][7]);   IMVL_AddElement(npmp,nmmm,B[2][7]);   IMVL_AddElement(nmpp,nmmm,B[4][7]);   IMVL_AddElement(nmmp,nmmm,B[6][7]);
                                                                                                                                                                        
               IMVL_AddElement(nppm,nppp,B[1][0]);   IMVL_AddElement(npmm,nppp,B[3][0]);   IMVL_AddElement(nmpm,nppp,B[5][0]);   IMVL_AddElement(nmmm,nppp,B[7][0]);
               IMVL_AddElement(nppm,nppm,B[1][1]);   IMVL_AddElement(npmm,nppm,B[3][1]);   IMVL_AddElement(nmpm,nppm,B[5][1]);   IMVL_AddElement(nmmm,nppm,B[7][1]);
               IMVL_AddElement(nppm,npmp,B[1][2]);   IMVL_AddElement(npmm,npmp,B[3][2]);   IMVL_AddElement(nmpm,npmp,B[5][2]);   IMVL_AddElement(nmmm,npmp,B[7][2]);
               IMVL_AddElement(nppm,npmm,B[1][3]);   IMVL_AddElement(npmm,npmm,B[3][3]);   IMVL_AddElement(nmpm,npmm,B[5][3]);   IMVL_AddElement(nmmm,npmm,B[7][3]);
               IMVL_AddElement(nppm,nmpp,B[1][4]);   IMVL_AddElement(npmm,nmpp,B[3][4]);   IMVL_AddElement(nmpm,nmpp,B[5][4]);   IMVL_AddElement(nmmm,nmpp,B[7][4]);
               IMVL_AddElement(nppm,nmpm,B[1][5]);   IMVL_AddElement(npmm,nmpm,B[3][5]);   IMVL_AddElement(nmpm,nmpm,B[5][5]);   IMVL_AddElement(nmmm,nmpm,B[7][5]);
               IMVL_AddElement(nppm,nmmp,B[1][6]);   IMVL_AddElement(npmm,nmmp,B[3][6]);   IMVL_AddElement(nmpm,nmmp,B[5][6]);   IMVL_AddElement(nmmm,nmmp,B[7][6]);
               IMVL_AddElement(nppm,nmmm,B[1][7]);   IMVL_AddElement(npmm,nmmm,B[3][7]);   IMVL_AddElement(nmpm,nmmm,B[5][7]);   IMVL_AddElement(nmmm,nmmm,B[7][7]);
            }

        /* Учёт краевых условий первого рода*/
        for(i=0 ; i<=nx ; i++)
          for(k=0 ; k<=ny ; k++)
          {
            l = 0;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));

            l = nz;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));
          }

        for(k=0 ; k<=ny ; k++)
          for(l=0 ; l<=nz ; l++)
          {
            i = 0;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));

            i = nx;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));
          }

        for(i=0 ; i<=nx ; i++)
          for(l=0 ; l<=nz ; l++)
          {
            k = 0;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));

            k = ny;
            x = i*hx - sx/2.0;
            y = k*hy - sy/2.0;
            z = l*hz - sz/2.0;

            IMVL_AddElement(index(i,k,l),index(i,k,l),1.0e+20);
            IMVL_AddElementVector(index(i,k,l),1.0e+20 * exact(x,y,z));
          }

	/* Сохранение матрицы */
	IMVL_Save();

        /* Решение СЛАУ */
        IMVL_ConjGradientMethod(5*nx*ny*nz);

        /* Сохранение решения */
        FSolve = fopen("solve.dat","w");

	tmp  = 0.0;
        tmp2 = 0.0;
        for(i=0 ; i<=nx ; i++)
          for(k=0 ; k<=ny ; k++)
            for(l=0 ; l<=nz ; l++)
            {
               x = i*hx - sx/2.0;
               y = k*hy - sy/2.0;
               z = l*hz - sz/2.0;

  	       tmp  += (IMVL_Solve[index(i,k,l)]-exact(x,y,z))*(IMVL_Solve[index(i,k,l)]-exact(x,y,z));
               tmp2 += exact(x,y,z) * exact(x,y,z);
            }
        printf("Residual = %le\n",sqrt(tmp)/sqrt(tmp2));


        for(i=0 ; i<=nx ; i++)
          for(k=0 ; k<=ny ; k++)
            for(l=0 ; l<=nz ; l++)
            {
               x = i*hx - sx/2.0;
               y = k*hy - sy/2.0;
               z = l*hz - sz/2.0;

               fprintf(FSolve,"%lf %lf %lf %lf %lf\n",x,y,z,IMVL_Solve[index(i,k,l)],exact(x,y,z));
            }

        fclose(FSolve);

        return 0;
}
