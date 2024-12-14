package plotd;
//Don't touch this code this belongs to Milan and Mamdouh
//FFTcta: Cooley-Tukey FFT calculation

public class FFTcta {

	private final int MAX_POINTS = 1024, MAX_PTS2 = MAX_POINTS / 2;
	private double si_re[] = new double[MAX_POINTS];
	private double si_im[] = new double[MAX_POINTS];
	private double so_re[] = new double[MAX_POINTS];
	private double so_im[] = new double[MAX_POINTS];
	private double WNn_re[] = new double[MAX_PTS2];
	private double WNn_im[] = new double[MAX_PTS2];
	private int NN, NN2, NNlog2;
	// sample data
	private final double T_s = 1.0 / 8.0E3;
	private final double f1 = 1000.0, f2 = 2000.0;
	private int N_SAMPLES = 0;
	private double xx[] = new double[MAX_POINTS];



	private void TestData() {
		int k;
		System.out.printf("test data(%d samples)\n", N_SAMPLES);
		System.out.printf("  k     data\n");
		System.out.printf("__________________\n");
		for(k=0;k<N_SAMPLES;k++) {
			xx[k] = Math.sin(2.0*Math.PI*f1*T_s*k)+0.5*Math.sin(2.0*Math.PI*f2*T_s*k+3.0*Math.PI/4.0);
			System.out.printf("%5d   %10.4f \n", k,xx[k]);
		}
		System.out.printf("________________\n");
	}
	
	//bit reverse algorithm
	private int BitRev(int k) {
		int i, k_mask, br_mask, k_br;
		
		k_br = 0;
		k_mask = 1;
		br_mask = (1 << NNlog2-1);
		for (i=0 ; i < NNlog2; i++) {
			if ((k & k_mask) != 0) {
				k_br |= br_mask; 
			}
			k_mask = k_mask<<1;
			br_mask = br_mask>>1;
		}
		return k_br;
	}
	
	private boolean FftInit() {
		int k, n;
		NN = N_SAMPLES;
		NN2 = NN /2;
		NNlog2 = 0;
		k=NN;
		while(k!=1) {
			if(k%2==1) return false;
			k /=2;
			NNlog2++;
		}
		System.out.printf("test data(%d samples)\n", N_SAMPLES);
		System.out.printf("stages : %d\n",NNlog2);
		System.out.printf("Wn Factors.....\n");
		System.out.printf("    n    Re{WNn}    Im{WNn}\n");
		System.out.printf("________________________________\n");
		for(n = 0; n <NN2; n++){
			WNn_re[n] = Math.cos(2.0*Math.PI*n/NN);
			WNn_im[n] = -Math.sin(2.0*Math.PI*n/NN);
			System.out.printf("%5d %10.4f %10.4f\n",n, WNn_re[n], WNn_im[n] );
		}
		System.out.printf("________________________________\n");
		System.out.printf("Bit reverse ordering of samples....\n");
		for (k = 0; k < NN ; k++) {
			//System.out.printf("k:%3d k_br :%3d\n", k, BitRev(k));
			so_re[k] = xx[BitRev(k)];
			so_im[k] = 0.0;
			System.out.printf(" so_re[%d]:%10.4f\n", k, so_re[k]);
		}
		return true;
	}
	
	private void FftStage(int stage) {
		// algorithm to compute the fft at different stages
		// for 8 bit fft the stages are 3
		int k = 0;
		for (k =0; k <NN; k++) {
			si_re[k] = so_re[k];
			si_im[k] = so_im[k];
		}
		// Butterfly algorithm
		System.out.printf(" ----------------\n");
		System.out.printf(" ---stage %2d ---\n", stage);
		System.out.printf(" ----------------\n");
		System.out.printf("  k		Re{si}		Im{si}		Re{so}		Im{so}\n");
		System.out.printf("-----------------------------------------------------------------------\n");
		for (k=0;k<NN; k++) {
			System.out.printf("%3d	    %10.4f	    %10.4f	    %10.4f	    %10.4f \n", k, si_re[k], si_im[k], so_re[k], so_im[k]);
		}
		System.out.printf("-----------------------------------------------------------------------\n");
	} 
	
	private void FftComp() {
		int stage;
		for (stage = 0; stage<NNlog2 ; stage++) {
			FftStage(stage);
		}
	}

	void ShowTestData() {
		
	}


	int FftGetLength() {
		return N_SAMPLES;
	}

	double FftGetReal(int k) {
		return so_re[k];
	}

	double FftGetImag(int k) {
		return so_im[k];
	}

	boolean FFTrun(int N_Samples)
	{
		boolean rc_success;
		N_SAMPLES = N_Samples;
		System.out.printf("N %d samples.\n", N_SAMPLES);
		rc_success = false;
		TestData();
		if(!(rc_success = FftInit())) {
			System.out.println(" *** length of power is not power of 2");
		} else{
			FftComp();
		}
		return rc_success;
	}

	public FFTcta() {
		System.out.println("--- FFT Cooley-Tukey Algorithm ---");
	}

}
