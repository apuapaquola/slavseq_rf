#ifndef SW_H
#define SW_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-10

int alfabeto[256]={['A']=1, ['C']=2, ['G']=3, ['T']=4, ['N']=5,
		   ['a']=1, ['c']=2, ['g']=3, ['t']=4, ['n']=5};

double matriz[6][6]=
{
  /*          -     A     C     G     T     N  */
  /* - */   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  /* A */   0.0,  1.0, -1.0, -1.0, -1.0, -1.0,
  /* C */   0.0, -1.0,  1.0, -1.0, -1.0, -1.0,
  /* G */   0.0, -1.0, -1.0,  1.0, -1.0, -1.0,
  /* T */   0.0, -1.0, -1.0, -1.0,  1.0, -1.0,
  /* N */   0.0, -1.0, -1.0, -1.0, -1.0, -1.0
};
double abre_gap=-5.0;
double estende_gap=-1.0;

void sw(char *s, char *t, SV ***_sp)
{
  /* Eq. de recorrencia em Gusfield, pag. 244 */
  
  int s_len=strlen(s);
  int t_len=strlen(t);


  /* Initialize the dynamic programming matrix */
  double v[(s_len+2)*(t_len+1)];
  double e[(s_len+2)*(t_len+1)];
  double f[(s_len+2)*(t_len+1)];

  v[0]=0;
  e[0]=f[0]=-1e+200;

  int i;
  for(i=1; i<=s_len+1; i++)
  {
    e[i]=-1e+200; /* -Inf */
    v[i]=f[i]=0; /* abre_gap + i*estende_gap; */
  }
  
  int j;
  for(j=1; j<=t_len; j++)
  {
    f[j*(s_len+2)]=-1e+200; /* -Inf */
    v[j*(s_len+2)]=e[j*(s_len+2)]=0; /* =abre_gap + j*estende_gap; */
  }
  
  for(j=1; j<=t_len; j++)
  {
    f[j*(s_len+2)+(s_len+1)]=-1e+200; /* -Inf */
    v[j*(s_len+2)+(s_len+1)]=e[j*(s_len+2)+(s_len+1)]=0; /* =abre_gap + j*estende_gap; */
  }

  /* Calculate the dynamic programming matrix */

  for(j=1; j<=t_len; j++)
    for(i=1; i<=s_len; i++)
    {
      double g=v[(j-1)*(s_len+2)+(i-1)]+matriz[alfabeto[(int)s[i-1]]][alfabeto[(int)t[j-1]]];
      
      e[(j)*(s_len+2)+(i)]=fmax(e[(j-1)*(s_len+2)+(i)],
				v[(j-1)*(s_len+2)+(i)]+abre_gap)+estende_gap;
	
      f[(j)*(s_len+2)+(i)]=fmax(f[(j)*(s_len+2)+(i-1)],
				v[(j)*(s_len+2)+(i-1)]+abre_gap)+estende_gap;
	
      double tmp=fmax(e[(j)*(s_len+2)+(i)],f[(j)*(s_len+2)+(i)]);
      v[(j)*(s_len+2)+(i)]=fmax(tmp,g);
    }

  /* Calculate the alignment */
  char s_ali[s_len+t_len+1]; s_ali[s_len+t_len]='\0';
  char t_ali[s_len+t_len+1]; t_ali[s_len+t_len]='\0';
  char m_ali[s_len+t_len+1]; m_ali[s_len+t_len]='\0';

  double score=v[(t_len)*(s_len+2)+(1)];
  int i_max=1;
  int j_max=t_len;
  
  for(i=2; i<=s_len+1; i++)
    if(score<v[(t_len)*(s_len+2)+(i)])
      {
	score=v[(t_len)*(s_len+2)+(i)];
	i_max=i;
	j_max=t_len;
      }

  if(i_max==s_len+1) j_max=0;

  for(j=1; j<=t_len-1; j++)
    if(score<v[(j)*(s_len+2)+(s_len)])
      {
	score=v[(j)*(s_len+2)+(s_len)];
	i_max=s_len;
	j_max=j;
      }

   
  i=s_len; 
  j=t_len; 
  int k=s_len+t_len-1;

  while(i>i_max || j>j_max)
    if(i>i_max)
      {
	s_ali[k]=s[i-1];
	t_ali[k]='-';
	m_ali[k]=' ';
	i--; k--;
      }
    else
      {
	s_ali[k]='-';
	t_ali[k]=t[j-1];
	m_ali[k]=' ';
	j--; k--;
      }
  
  int flag=0;

  while(i>0 && j>0)
  {
    double g=v[(j-1)*(s_len+2)+(i-1)]+matriz[alfabeto[(int)s[i-1]]][alfabeto[(int)t[j-1]]];
      
    //      printf("%d %d %3.0f %d\n",i,j,v[(j)*(s_len+2)+(i)],flag);
      
    if(abs(v[(j)*(s_len+2)+(i)] - g)<EPSILON && flag==0)
    {
      s_ali[k]=s[i-1];
      t_ali[k]=t[j-1];
      m_ali[k]=alfabeto[(int)s[i-1]]==alfabeto[(int)t[j-1]]?'|': matriz[alfabeto[(int)s[i-1]]][alfabeto[(int)t[j-1]]] > 0?' ':' ';
      i--; j--; k--;
    }
    else if(abs(v[(j)*(s_len+2)+(i)] - e[(j)*(s_len+2)+(i)])<EPSILON || flag==-1)
    {
      if(abs(e[(j)*(s_len+2)+(i)] - (v[(j-1)*(s_len+2)+(i)]+abre_gap+estende_gap))<EPSILON)
	flag=0;
      else
	flag=-1;
      
      s_ali[k]='-';
      t_ali[k]=t[j-1];
      m_ali[k]=' ';
      j--; k--;
    }
    else if(abs(v[(j)*(s_len+2)+(i)] - f[(j)*(s_len+2)+(i)])<EPSILON || flag==-2)
    {
      if(abs(f[(j)*(s_len+2)+(i)] - (v[(j)*(s_len+2)+(i-1)]+abre_gap+estende_gap))<EPSILON)
	flag=0;
      else
	flag=-2;
	  
      s_ali[k]=s[i-1];
      t_ali[k]='-';
      m_ali[k]=' ';
      i--; k--;
    }
    else
    {
      printf("Erro nas comparações de ponto flutuante: favor aumentar o valor de EPSILON e recompilar.\n");
      exit(1);
    }
  }
  
  
  while(i>0 || j>0)
    if(i>0)
    {
      s_ali[k]=s[i-1];
      t_ali[k]='-';
      m_ali[k]=' ';
      i--; k--;
    }
    else
    {
      s_ali[k]='-';
      t_ali[k]=t[j-1];
      m_ali[k]=' ';
      j--; k--;
    }


  /* Return the alignment as three values on perl stack */ 

  SV **sp=*_sp;
  XPUSHs(sv_2mortal(newSVpv(s_ali+k+1,0)));
  XPUSHs(sv_2mortal(newSVpv(m_ali+k+1,0)));
  XPUSHs(sv_2mortal(newSVpv(t_ali+k+1,0)));
  XPUSHs(sv_2mortal(newSVnv(score)));
  *_sp=sp;
}

#endif
