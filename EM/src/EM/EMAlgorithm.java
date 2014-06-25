package EM;

import java.util.List;

public class EMAlgorithm {

	static String inFileAddress = "emdata.txt";
	
	// 最初に２つの、中心と分散が異なる正規分布を重ね合わせたようなデータ分布のデータを作っておく
	// いったん２つの１次元の正規分布で近似するようにしてみる。
	public static void main(String[] args){
		
		// 近似する正規分布の数
		int gNum =2;
		
		// 正規分布に対する係数
		double[] w = new double[gNum];
		
		// 正規分布の平均値
		double[] mu = new double[gNum];
		
		// 正規分布の分散
		double[] sig = new double[gNum];
		
		// 初期データ
		w[0] = 0.2;
		w[1] = 0.8;
		mu[0] = -2.0;
		mu[1] = 30.0;
		sig[0] = 1.0;
		sig[1] = 50.0;
		
		// 読み込んだデータを格納する箱（１次元） xNumはデータの個数
		double[] x;
		int xNum;
		List<double[]> list = FileUtil.getDoublesByQuot(inFileAddress);
		xNum = list.size();
		x = new double[xNum];
		for(int i = 0; i < xNum; i++){
			x[i] = list.get(i)[0];
		}
		
		
		// 媒介変数イータ (eta[i][l] ・・・iは読込データの順、lは正規分布の順)
		double[][] eta = new double[xNum][gNum];
		
		for(int i = 0; i < 10000; i++){
			eta = eStep(w, mu, sig, x);
			w = mStepW(eta, gNum, xNum);
			mu = mStepMu(eta, gNum, xNum, x);
			sig = mStepSig(eta, gNum, xNum, x, mu);
		}
		
		for(int i = 0; i < gNum; i++){
			System.out.println("weight" + i + ": " + w[i]);
			System.out.println("means" + i + ": " + mu[i]);
			System.out.println("sigma" + i + ": " + sig[i]);
		}
		
		
	}
	
	
	
	static double gauss(double x, double mu, double sig){
		return 1.0 / Math.sqrt(2.0 * 3.1415 * sig) * Math.exp((-1.0) * Math.pow((x - mu), 2.0) / (2.0 * sig));
		
	}
	
	static double[][] eStep(double[] w, double[] mu, double[] sig, double[] x){
		
		double[][] returnEta = new double[x.length][w.length];
		
		for(int i = 0; i < x.length; i++){
			
			// イータの分母をまず求める
			double sumA = 0.0;
			for(int k = 0; k < w.length; k++){
				sumA = sumA + w[k] * gauss(x[i], mu[k], sig[k]);
			}
			for(int k = 0; k < w.length; k++){
				returnEta[i][k] = w[k] * gauss(x[i], mu[k], sig[k]) / sumA;
			}
			
			
		}
		
		return returnEta;
		
	}
	
	static double[] mStepW(double[][] eta, int gNum, int xNum){
		double[] returnW = new double[gNum];
		
		for(int i = 0; i < gNum; i++){
			double sumW = 0.0;
			for(int k = 0; k < xNum; k++){
				sumW = sumW + eta[k][i];
			}
			returnW[i] = sumW / (double)(xNum);
		}
		return returnW;
	}
	
	static double[] mStepMu(double eta[][], int gNum, int xNum, double[] x){
		double[] returnMu = new double[gNum];
		for(int i = 0; i < gNum; i++){
			double sumA = 0.0;
			double sumB = 0.0;
			for(int k = 0; k < xNum; k++){
				sumA = sumA + eta[k][i] * x[k];
				sumB = sumB + eta[k][i];
			}
			returnMu[i] = sumA / sumB;
		}
		return returnMu;
	}
	
	static double[] mStepSig(double[][] eta, int gNum, int xNum, double[] x, double[] mu){
		double[] returnSig = new double[gNum];
		for(int i = 0; i < gNum; i++){
			double sumA = 0.0;
			double sumB = 0.0;
			for(int k = 0; k < xNum; k++){
				sumA = sumA + eta[k][i] * Math.pow(x[k]-mu[i], 2);
				sumB = sumB + eta[k][i];
			}
			returnSig[i] = sumA / sumB;
		}
		return returnSig;
	}
	
}
