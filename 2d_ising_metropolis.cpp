#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <limits>
#include <sys/stat.h>
#include <omp.h>

using namespace std;
#define SIZE 100
#define NCELL 10000
#define MAIN_DEFAULT
#define MAKE_PICTURE
#define METROPOLIS

double sigma[SIZE*SIZE];
int ir = 12345;
// ランダム数を返す関数
double myrand(){
  int    add = 2147483647, mul = 48828125;
  double nrm = 2147483648.0, rnd;
  ir *= mul;
  if (ir<0) ir += add+1;
  rnd = ir/nrm;
  return rnd;
}

// 初期化関数 すべて上向きスピンに
void init(){
  #pragma omp parallel for
  for(int i=0; i<NCELL; i++){
    sigma[i] = 1.0;
  }
}
// 乱数振って上向き下向きを変える関数
void sweep(double beta){
  int i, j;
  double mean, val, prob, myrnd;
  for(i=0; i<SIZE; i++){
    for(j=0; j<SIZE; j++){
      // n1,n2の格子点にたいして隣り合う4つの格子点を評価
      mean = sigma[i*SIZE+(j+1)%SIZE] + sigma[i*SIZE+(j-1+SIZE)%SIZE]
      + sigma[((i+1)%SIZE)*SIZE+j] + sigma[((i-1+SIZE)%SIZE)*SIZE+j];
      val = exp(beta*mean);
      prob = val/(val+1.0/val);
      myrnd = myrand();
      if(myrnd<prob){
        sigma[i*SIZE+j] = 1.0;
      }else{
        sigma[i*SIZE+j] = -1.0;
      }
    }
  }
}
// メトロポリス法でのsweep
void sweep_metropolis(double beta){
  int i, j;
  double mean;
  double Old_energy, New_energy, Old_sigma, New_sigma, dE;
  for(i=0; i<SIZE; i++){
    for(j=0; j<SIZE; j++){
      // 今のスピンと逆向きを設定
      Old_sigma =   sigma[i*SIZE+j];
      New_sigma = - Old_sigma;
      // n1,n2の格子点にたいして隣り合う4つの格子点を評価
      mean = sigma[i*SIZE+(j+1)%SIZE] + sigma[i*SIZE+(j-1+SIZE)%SIZE]
             + sigma[((i+1)%SIZE)*SIZE+j] + sigma[((i-1+SIZE)%SIZE)*SIZE+j];
      // エネルギーの計算
      Old_energy = - Old_sigma * mean;
      New_energy = - New_sigma * mean;
      dE = beta * (New_energy - Old_energy);
      // アップデート条件文
      if((dE<=0.0) || (exp(-dE)>myrand())){
        sigma[i*SIZE+j] = New_sigma;
      }
    }
  }
}
// 磁場を計算する
double magnet(){
  int i;
  double mag_temp = 0.0;
  #pragma omp parallel for reduction(+:mag_temp)
  for(i=0; i<NCELL; i++){
    mag_temp += sigma[i];
  }
  mag_temp /= NCELL;
  return mag_temp;
}

#ifdef MAIN_DEFAULT
  int main(void){
    // 定数定義
    double beta;
    int maxiter, iter = 0;
    cout << "# Ising model simulation" << endl;
    cout << "input initial beta:";
    cin >> beta;
    cout << "input the number of iteration:";
    cin >> maxiter;
    cout << "size:" << SIZE << " beta:" << beta << " maxiter" << maxiter << endl;

    //保存用フォルダの作成
    char dirname[50], _dirname[50];
    char filename[50];
    FILE * fpg;
    #ifdef METROPOLIS
      sprintf(dirname, "data_ising_metropolis");
    #else
      sprintf(dirname, "data_ising");
    #endif
    strcpy(_dirname,dirname);
    mkdir(dirname, 0777);
    chmod(dirname, 0777);

    //初期値出力
    init();
    sprintf(filename, "/magnet_%.2f.dat", beta);
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
      outputfile << iter << "\t" << magnet() << endl;
    outputfile .close();
    // iterationの開始
    for(iter=1; iter<=maxiter; iter++){
      #ifdef METROPOLIS
        sweep_metropolis(beta);
      #else
        sweep(beta);
      #endif
      ofstream outputfile(_dirname, ios_base::app);
        outputfile << iter << "\t" << magnet() << endl;
      outputfile .close();
    }
    #ifdef MAKE_PICTURE
      fpg = popen("gnuplot", "w");
      char loadname[50], picname[50];
      #ifdef METROPOLIS
        sprintf(loadname, "'./data_ising_metropolis/magnet_%.2f.dat'\n", beta);
        sprintf(picname, "'./data_ising_metropolis/magnet_%.2f.png'\n", beta);
      #else
        sprintf(loadname, "'./data_ising/magnet_%.1f.dat'\n", beta);
        sprintf(picname, "'./data_ising/magnet_%.1f.png'\n", beta);
      #endif
      if(fpg == NULL) return -1;
        fputs("set output ", fpg);
        fputs(picname, fpg);
        fputs("loadname = ", fpg);
        fputs(loadname, fpg);
        fputs("load './gnu_ising_magnet.gp'\n",fpg);
        fflush(fpg);
      pclose(fpg);
    #endif
    return 0;
  }
#endif
