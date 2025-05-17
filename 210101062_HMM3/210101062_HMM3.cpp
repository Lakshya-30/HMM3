// 210101062_HMM3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;

#define M 32	//Number of observation symbols per state
#define N 5		//Number of states
const int MAX_T = 160;		//Maximum number of observations

//Model parameters A, B and Pi
long double A[N+1][N+1] = {0};
long double B[N+1][M+1] = {0};
long double Pi[N+1] = {0};

int T=1;			//Number of observations
int O[MAX_T+1];		//Observation sequence
long double alpha[MAX_T+1][N+1];
long double beta[MAX_T+1][N+1];
long double delta[MAX_T+1][N+1];
int psi[MAX_T+1][N+1];
long double P_O_given_lambda = 0;
long double p_star = 0;
int q_star[MAX_T+1];	//Optimal state sequence

long double gamma[MAX_T+1][N+1];
long double xi[MAX_T + 1][N+1][N+1];
long double A_bar[N+1][N+1] = {0};			//Reestimated model parameters
long double B_bar[N+1][M+1] = {0};
long double Pi_bar[N+1] = {0};

//Function to calculate alpha variable to find the solution of problem number 1
void forwardProcedure(){
	long double sum;
	int obs = O[1];
	P_O_given_lambda = 0;
	for(int i=1; i<=N; i++){			//initialization
		alpha[1][i] = Pi[i]*B[i][obs];
	}
	
	for(int t=1; t<T; t++){				//induction
		for(int j=1; j<=N; j++){
			sum = 0;
			for(int i=1; i<=N; i++){
				sum += (alpha[t][i] * A[i][j]);
			}
			alpha[t + 1][j] = (sum * B[j][O[t + 1]]);
			//cout<<alpha[t+1][j]<<"    ";
		}
		//cout<<'\n';
	}

	for(int i=1; i<=N; i++){			//termination
		P_O_given_lambda += alpha[T][i];
	}
}

//Function to calculate beta variable.
void backwardProcedure(){
	long double sum;
	for(int i=1; i<=N; i++){		//initialization
		beta[T][i] = 1.0;
	}
	for(int t=T-1; t>=1; t--){		//induction
		int obs = O[t+1];
		for(int i=1; i<=N; i++){
			sum = 0;
			for(int j=1; j<=N; j++){
				sum += B[j][obs]*A[i][j]*beta[t+1][j];
			}
			beta[t][i]=sum;
			//cout<<beta[t][i]<<"    ";
		}
		//cout<<"\n";
	}
}

//Function to find the single best state sequence for the given observation sequence
void viterbiAlgorithm(){
    for(int i=1; i<=N; i++){			//Initialization
        delta[1][i] = Pi[i] * B[i][O[1]];
        psi[1][i] = 0;
    }

	for(int j=1; j<=N; j++){			//Induction
		for(int t=2; t<=T; t++){
            long double maxVal = 0, temp = 0;
            int index = 0;
            for(int i=1; i<=N; i++){
                temp = (delta[t-1][i] * A[i][j]);
                if(temp > maxVal){
					maxVal = temp;
					index = i;
				}
            }
            delta[t][j] = maxVal * B[j][O[t]];
			psi[t][j] = index;
        }
    }

    long double max = 0;
    for(int i=1; i<=N; i++){			//Termination
        if(delta[T][i] > max) {
			max = delta[T][i];
			q_star[T] = i;
		}
        p_star = max;
    }

    for(int t=T-1; t>0; t--){			//Path (state sequence) backtracking
        q_star[t] = psi[t+1][q_star[t+1]];
    }

	cout<<"Probability: "<<p_star<<endl;
	for(int t=1; t<=T; t++){
        cout<<q_star[t]<<" ";
    }
	cout<<endl;
}

//Function to calculate xi variable
void computeXi(){
	long double denominator[MAX_T + 1];
	for(int t=1; t<=T; t++){
		denominator[t] = 0;
		for(int i=1; i<=N; i++){
			for(int j=1; j<=N; j++){
				denominator[t] += alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
			}
		}
		for(int i=1; i<=N; i++){
			long double z;
			for(int j=1; j<=N; j++){
				z = alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
				xi[t][i][j]= z/denominator[t];
			}
		}
	}
}

//Function to compute gamma values
void computeGamma(){
	for (int t=1; t<=T-1; t++){
		for (int i=1; i<=N; i++){
			gamma[t][i] = 0;
			for (int j=1; j<=N; j++){
				gamma[t][i] += xi[t][i][j];
			}
		}
	}
}

//Function to find solution to the reestimation problem
void reestimateParameters(){
	for(int i=1; i<=N; i++){
		Pi_bar[i] = gamma[1][i];			//Reevaluating Pi
	}
	for(int i=1; i<=N; i++){				//Reevaluating A
		for(int j=1; j<=N; j++){
			long double num = 0, denom = 0;
			for(int t=1; t<=T-1; t++){
				num += xi[t][i][j];
				denom += gamma[t][i];
			}
			A_bar[i][j] = num/denom;
		}
	}
	for(int j=1; j<=N; j++){				//Reevaluating B
		for(int k=1; k<=M; k++){
			long double num = 0, denom = 0;
			for(int t=1; t<T; t++){
				denom += gamma[t][j];
				if(O[t]==k){
					num += gamma[t][j];
				}
			}
			B_bar[j][k] = num/denom;
		}
	}
}

int main(){
	FILE *f1 = NULL;    // Read the data from the input file.
	double x;
    int err = fopen_s(&f1, "Files/Initial_Model.txt", "r");
    if(err != NULL)     //error handling for file opening
    {
        printf("\nCannot open the file\n");
        system("pause");
        exit(1);
    }

	int count = 0;
	int i=1, j=1;
	while( !feof(f1) )
	{
		if( fscanf_s( f1, "%lf", &x) == 1)
		{
			if(count<25){
				A[i][j] = x;
				j++;
				if(j==N+1){
					i++,j=1;
					if(i==N+1) i=1,j=1;
				}
				count++;
			}
			else if(count<185){
				B[i][j] = x;
				j++;
				if(j==M+1){
					i++,j=1;
					if(i==N+1) i=1,j=1;
				}
				count++;
			}
			else{
				Pi[i]=x;
				i++;
			}
		}
		else{
			char lines[500];
			fgets(lines, 500, f1);
		}
	}
	fclose(f1);

	FILE* f2= fopen("Files/HMM_OBSERVATION_SEQUENCE_1.txt","r");
    while( !feof(f2) )
	{
		if( fscanf_s( f2, "%lf", &x) == 1){
			O[T++]=x;
		}
	}
	fclose(f2);
	T--;

	for(int iter=1; iter<5; iter++){
		cout<<"-------------------Iteration "<<iter<<"------------------------\n";
		forwardProcedure();
		backwardProcedure();
		viterbiAlgorithm();
		computeXi();
		computeGamma();
		reestimateParameters();				//Overwriting the current model
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				A[i][j]= A_bar[i][j];
			}
			for(int k=0; k<M; k++){
				B[i][k] = B_bar[i][k];
			}
			Pi[i] = Pi_bar[i];
		}
	}

	return 0;
}

