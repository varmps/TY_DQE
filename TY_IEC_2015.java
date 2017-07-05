import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.frame.*;
import javax.swing.*;
import java.awt.event.*;
import java.lang.*;
import java.util.*;
//
import ij.process.FHT.*;


public class TY_IEC_2015 extends PlugInFrame implements ActionListener,FocusListener 
{
	ImagePlus imp,impcut,impffttot;
	ImageProcessor ipffttot,ip,ipcut;
	float C[],S[];
	double IMean,mode1,mode2,dprof[],xprof[];
	int imamode,nsamp=2048,cutsize=128;
	int nnpsleave=1,nnpscount=1;
	double a,b,pixsize=0.139;
	double slope = 1;
	double intercept = 0;
	
	//subpixsize=0.144/32; 
	/* these are +/- default values, applicable for the systems tested upto now */
	JTextField TFpixsize,TFsubpixsize,TFcutsize,TFnnpsleave,TFnnpscount;
	JLabel Lpixsize,Lsubpixsize,Lproflength,Lstatus,Lcutsize,Lnnpsleave,Lnnpscount;
	JTextField TFproflength;
	JRadioButton CBlin,CBsqrt,CBlog,CBproc,CBprofile,CBprof2mtf,CBangcorr,CBshowNNPS;
	JButton Bmtfhor,Bmtfvert,Bnps;
	JButton Bcopycol1,Bcopycol2,Bcopycol3,Bcopycol4,Bcopycol5,Bcopycol6,Bcopycol7,Bcopycol8,Bcopy;
	JTable JTa;
	
	public TY_IEC_2015() 
	{	
		
		super("TY_IEC_2015");
	//	setBackground(SystemColor.control);
		setLayout(new BorderLayout());
		JPanel panel1=new JPanel();
		panel1.setLayout(new GridLayout(12,3));
	//	panel1.setBackground(SystemColor.control);
	
	//------First row------
		panel1.add(new JLabel("Actual pixel size (mm):"));
		TFpixsize = new JTextField(Double.toString(pixsize));
		TFpixsize.setToolTipText("actual pixel size in the original image");
		panel1.add(TFpixsize);
		TFpixsize.addActionListener(this);
		TFpixsize.addFocusListener(this);
		Lpixsize=new JLabel("");
		panel1.add(Lpixsize);
		
	//------Second row------
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		
	//------Third row-------
		panel1.add(new JLabel("MTF subpixel size (mm):"));
		TFsubpixsize = new JTextField();
		//panel1.add(new JSeparator());
		TFsubpixsize.setToolTipText("Subsampled pixel size for MTF-calculation");
		TFsubpixsize.addFocusListener(this);
		TFsubpixsize.addActionListener(this);
		panel1.add(TFsubpixsize);
		Lsubpixsize=new JLabel("");
		panel1.add(Lsubpixsize);
	
	//------Forth row-------
		panel1.add(new JLabel("Profile length:"));	
		TFproflength = new JTextField(Integer.toString(nsamp));
		TFproflength.setToolTipText("Number of subpixels in the line spread function");
		TFproflength.addActionListener(this);
		TFproflength.addFocusListener(this);
		panel1.add(TFproflength);
		
	//	Lproflength=new JLabel("");
	//  panel1.add(Lproflength);
	
		CBangcorr=new JRadioButton("Correct MTF for angle",true);
		panel1.add(CBangcorr);
		
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		
		TFcutsize = new JTextField(Integer.toString(cutsize));
		TFcutsize.setToolTipText("Length of LSF in subpixels. Power of 2 is required!");
		panel1.add(new JLabel("NNPS fragment size:")); 
		TFcutsize.addActionListener(this); 
		TFcutsize.addFocusListener(this);
		panel1.add(TFcutsize);
		
		Lcutsize=new JLabel("");
		panel1.add(Lcutsize);
		
		TFnnpsleave = new JTextField(Integer.toString(nnpsleave));
		TFnnpsleave.setToolTipText("don't use ... central lines in center for NNPS");
		panel1.add(new JLabel("Don't use .. lines for NNPS:"));
		TFnnpsleave.addActionListener(this);
		TFnnpsleave.addFocusListener(this);
		panel1.add(TFnnpsleave);
		
		Lnnpsleave=new JLabel("");
		panel1.add(Lnnpsleave);

		TFnnpscount = new JTextField(Integer.toString(nnpscount));
		TFnnpscount.setToolTipText("don't use ... central lines in center for NNPS");
		panel1.add(new JLabel("Use .. lines for NNPS:"));
		TFnnpscount.addActionListener(this);
		TFnnpscount.addFocusListener(this);
		
		panel1.add(TFnnpscount);
	//	Lnnpscount=new JLabel("");
	//	panel1.add(Lnnpscount);
		CBshowNNPS=new JRadioButton("Show 2D NNPS",false);
		panel1.add(CBshowNNPS);
		
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		panel1.add(new JSeparator());
		
		panel1.add(new JLabel("Image type:"));
		ButtonGroup BG=new ButtonGroup();
		CBlin	=new JRadioButton ("Linear",true);
		CBsqrt	=new JRadioButton ("Square Root");
	//	CBlog	=new JRadioButton ("Logarithmic");
		BG.add(CBlin);
		BG.add(CBsqrt);
	//	BG.add(CBlog);
		panel1.add(CBlin);
		panel1.add(CBsqrt);
	//	panel1.add(new JLabel(""));
	//	panel1.add(CBlog);
	//	panel1.add(new JLabel(""));
	
	//	panel1.add(new JSeparator());
	//	panel1.add(new JSeparator());
	//	panel1.add(new JSeparator());
		
//ei tee midagi hetkel:
		CBprofile=new JRadioButton("Calculate LSF");
		CBprof2mtf=new JRadioButton("Calculate MTF from LSF",true);

		panel1.add(new JLabel("MTF-mode:"));
		panel1.add(CBprofile);
		panel1.add(CBprof2mtf);
	//	panel1.add(new JLabel(""));
		
		ButtonGroup BGmtf=new ButtonGroup();
		BGmtf.add(CBprofile);
		BGmtf.add(CBprof2mtf);
		Bmtfhor=new JButton("MTF hor");
		Bmtfhor.setToolTipText("Select a ROI including a vertical edge, angle<10° recommended");
		Bmtfhor.addActionListener(this);
		panel1.add(Bmtfhor);
		Bmtfvert=new JButton("MTF vert");
		Bmtfvert.setToolTipText("Select a ROI including a horizontal edge, angle<10° recommended");
		Bmtfvert.addActionListener(this);
		panel1.add(Bmtfvert);
		Bnps=new JButton("NNPS");
		Bnps.setToolTipText("Select a ROI with uniform density");
		Bnps.addActionListener(this);
		panel1.add(Bnps);
		
		add(panel1,BorderLayout.NORTH);
		
		JTa=new JTable(5000,8);
		JTa.setPreferredScrollableViewportSize(new Dimension(50,100));
		JTa.setEnabled(false);
		JScrollPane sp = new JScrollPane(JTa);
		add(sp,BorderLayout.CENTER);
		JPanel panel6=new JPanel();
		panel6.setLayout(new GridLayout(2,1));
		JPanel panel5=new JPanel();
		panel5.setLayout(new GridLayout(1,7));
		Bcopycol1=new JButton("copy A");Bcopycol1.addActionListener(this);panel5.add(Bcopycol1);
		Bcopycol2=new JButton("copy B");Bcopycol2.addActionListener(this);panel5.add(Bcopycol2);
		Bcopycol3=new JButton("copy C");Bcopycol3.addActionListener(this);panel5.add(Bcopycol3);
		Bcopycol4=new JButton("copy D");Bcopycol4.addActionListener(this);panel5.add(Bcopycol4);
		Bcopycol5=new JButton("copy E");Bcopycol5.addActionListener(this);panel5.add(Bcopycol5);
		Bcopycol6=new JButton("copy F");Bcopycol6.addActionListener(this);panel5.add(Bcopycol6);
		Bcopycol7=new JButton("copy G");Bcopycol7.addActionListener(this);panel5.add(Bcopycol7);
		Bcopycol8=new JButton("copy H");Bcopycol8.addActionListener(this);panel5.add(Bcopycol8);
		panel6.add(panel5);
		JPanel panel4=new JPanel();
		panel4.setLayout(new GridLayout(1,3));
		Lstatus= new JLabel("status");
		panel4.add(Lstatus);
		Bcopy=new JButton("Copy all columns");Bcopy.setToolTipText("Copy the system clipboard, data can be pasted to e.g. Excel");
		Bcopy.addActionListener(this);panel4.add(Bcopy);
		panel6.add(panel4);
		add(panel6,BorderLayout.SOUTH);
		pack();show();
	}
	
	double parseDouble(String s) throws NumberFormatException 
	{
		Double d = new Double(s);
		return(d.doubleValue());
	}	
	
	public double MedPix(ImageProcessor ip,int i0,int j0,int s)
	{
		double mp;
		int w,h;
		w=ip.getWidth();
		h=ip.getHeight();
		int i,j,cnt=0;
		double[] data=new double[(2*s+1)*(2*s+1)];
		
		for(i=(i0-s);i<(i0+s+1);i++)
			for(j=(j0-s);j<(j0+s+1);j++)
				if((i>-1)&&(i<w)&&(j>-1)&&(j<h))
				{
					data[cnt]=slope*ip.getPixelValue(i,j)+intercept;
					cnt++;	
				}
		Arrays.sort(data,0,cnt-1);
		mp=data[cnt/2];
		return mp;	
	}
	
	public void edge_finder(ImageProcessor ip)
	{
		/*
			ip should contain an edge crossing only the upper and the lower edge of the ip:
			|------*--|    |--*------|
			|     *   |    |   *     |
			|    *    | or |    *    |
			|   *     |    |     *   |
			|--*------|    |------*--|
			
			
			0) pre-process image with a median-filter to eliminate dead pixels/dropouts/...
				(radius=2 -> 5x5 region)
				(if this is not applied: one deviating pixel may corrupt the line position detection).
			1) find for each row largest signal diff between 2 adjacent pixels
			2) record parameters for least squares fir method
				count	sx=sx+i
					sy=sy+j
					sxx=sxx+i*i
					sxy=sxy+i*j
					syy=syy+j*j
					n++
			3) apply least square fit method
		
		
		*/
		int medfiltersize=2;
		int regx[]=new int[ip.getHeight()]; 
		double sx=0.0,sxx=0.0,sxy=0.0,syy=0.0,sy=0.0;
		double val,diff,diffmed,valmed,rms=0.0;
		int ijump=0,ijumpmed=0;
		int n=0,i,j;
		for(j=0;j<ip.getHeight();j++)
		{
			valmed=0.0;
			
			for(i=1;i<ip.getWidth();i++)
			{
				diffmed=Math.abs((float)MedPix(ipcut,i,j,medfiltersize)-(float)MedPix(ip,i-1,j,medfiltersize));
				if (diffmed>valmed)
				{
					valmed=diffmed;
					ijumpmed=i;
				}	
				regx[j]=ijumpmed;
			}
			i=ijumpmed;
			sx+=i;sy+=j;sxx+=i*i;sxy+=i*j;syy=j*j;n++;
		}
		b=((double)n*sxy-sx*sy)/((double)n*sxx-sx*sx);
		a=sx/n-sy/n/b;
		
		IJ.write("edge straightness?");
		IJ.write("line: a="+a+"; b="+b);
	
		
		for(j=0;j<ip.getHeight();j++)
			rms+=(a+(double)j/b-(double)regx[j])*(a+(double)j/b-(double)regx[j]);
		rms=Math.sqrt(rms/(double)ip.getHeight());
		IJ.write("rms="+rms);
	}
	
	public double get_leftmode(ImageProcessor ip,double a,double b)
	{
		/*
			measure selection to get left "mode"
			*******--
			*    *  |
			*   *   |
			*  *    |
			* *     |
			**-------
		*/
		
		double mode;
		ImageStatistics imstat=new ImageStatistics();
		Roi oroi,nroi;
		int roix[]=new int[4];
		int roiy[]=new int[4];
		roix[0]=0;					roiy[0]=0;
		roix[1]=0;					roiy[1]=ipcut.getHeight()-1;
		roix[2]=(int)(a+(ipcut.getHeight()-1)/b);	roiy[2]=ipcut.getHeight()-1;
		roix[3]=(int)(a);				roiy[3]=0;
		nroi=new PolygonRoi(roix,roiy,4,impcut,Roi.POLYGON);
		
		impcut.setRoi(nroi);
		
		imstat=impcut.getStatistics();
		double step=(imstat.max-imstat.min)/255.0;
		imstat=impcut.getStatistics(ImageStatistics.MODE,256,imstat.min-step/2.0,imstat.max+step/2.0);
		mode=imstat.dmode;
		return mode;	
	}
	
	public double get_rightmode(ImageProcessor ip,double a,double b)
	{
		/*
			measure selection to get right "mode"
			-------***
			|     *  *
			|    *   *
			|   *    *
			|  *     *
			--********
		*/
		
		double mode;
		ImageStatistics imstat=new ImageStatistics();
		Roi oroi,nroi;
		int roix[]=new int[4];
		int roiy[]=new int[4];
		roix[0]=ipcut.getWidth();			roiy[0]=0;
		roix[1]=ipcut.getWidth();			roiy[1]=ipcut.getHeight()-1;
		roix[2]=(int)(a+(ipcut.getHeight()-1)/b);	roiy[2]=ipcut.getHeight()-1;
		roix[3]=(int)(a);				roiy[3]=0;
		nroi=new PolygonRoi(roix,roiy,4,impcut,Roi.POLYGON);
		impcut.setRoi(nroi);
		ipcut.setHistogramSize(10000);
		imstat=impcut.getStatistics();
		double step=(imstat.max-imstat.min)/255.0;
		imstat=impcut.getStatistics(ImageStatistics.MODE,256,imstat.min-step/2.0,imstat.max+step/2.0);
		mode=imstat.dmode;
		ipcut.setHistogramSize(256);
		return mode;	
	}
	
	public double[] subsample(ImageProcessor ip,double a, double b,int n_sub,int mode,int shift)
	{
		// Number of lines (rows) after which the edge has moved by one pixel (column)
		int N = Math.abs((int)Math.round(b));
		// Number of lines in the ESF, that depends on the n_sub so that that the total number of subpixels on the LSF is n_sub
		int Pikkus = (int)Math.round(n_sub/N)+1;
		// ESF
		double prof[]=new double[Pikkus*N];
		
		int i,j;
		
		//If the edge is counter-clockwise the ESF starts from the right
		if(b>0){
	
			//x-position of the edge
			for(i=0;i<Pikkus;i++)
			{
			 	//y-position of the edge
			  	for(j=shift;j<N+shift;j++)
			  	{
			  		//Pikkus devided by 2 before the edge and after the edge
			  		prof[j-shift+i*N]=slope*ip.getPixelValue((int)Math.round(a)+(int)Math.round(Pikkus/2)-i,j)+intercept;
			  		//IJ.write(""+prof[j-shift+i*N]);
			  		//IJ.write("(int)Math.round(a)+Pikkus/2-i "+((int)Math.round(a)+(int)Math.round(Pikkus/2)-i)+" "+j);
			  	}
	
			}
			
		}
		//If the edge is clockwise the ESF starts from the left
		else{
			//x-position of the edge 
			for(i=0;i<Pikkus;i++)
			{
			  
				//y-position of the edge	
			  	for(j=shift;j<N+shift;j++)
			  	{
			 		//Pikkus devided by 2 before the edge and after the edge
			  		prof[j-shift+i*N]=slope*ip.getPixelValue((int)Math.round(a)-Pikkus/2+i,j)+intercept;
			  	//	IJ.write(""+prof[j-shift+i*N]);
			  	//	IJ.write("(int)Math.round(a)+Pikkus/2+i "+((int)Math.round(a)+(int)Math.round(Pikkus/2)+i)+" "+j);
			  		
			  	}
			}
		}
		
		return prof;
	}
	//Differentiates LSF and drunkates it to power of 2 lenght
	public double[] diff (double[] prof, int n_sub, int mode){
		// LSF
		double dp[]=new double [n_sub];
		int i;
		
		// Differentiating by IEC kernel option 1. [–1, 0, 1] also cuts to n_sub length
		if (mode == 1){
			//IJ.write("dp mode "+mode);
			for(i=0;i<n_sub;i++)
			{
					
				dp[i]=-prof[i]+prof[i+2];
			//	IJ.write(""+dp[i]);
			}
		}
		// Differentiating by IEC kernel option 2. [–0.5, 0, 0.5] also cuts to n_sub length
		else if (mode == 2){
			//IJ.write("dp mode "+mode);
			for(i=0;i<n_sub;i++)
			{
			
				dp[i]=-0.5*prof[i]+0.5*prof[i+2];
			}
		}
		// Differentiating by kernel [–1, 1] also cuts to n_sub length
		else{
			//IJ.write("mode "+mode);
			for(i=0;i<n_sub;i++)
			{
			
				dp[i]=-prof[i]+prof[i+1];
			}
		}
		return dp;
	}
	
	
	double[] AbsDFT(double[] adX) 
	{
		/*
			fourier transform for deriving mtf from an lsf-array
			this algorthm was taken and adapted from the measure_mtf_ plugin available on the ImageJ-site
		*/
		int iNK = adX.length;
		double[] adFTX = new double[iNK];
		double dA, dB, dT, dTheta;
		for (int k=0; k<iNK; k++) 
		{
			dA = 0;
			dB = 0;
			dT = -2.0 * Math.PI * k / iNK;
			for (int x=0; x<iNK; x++) 
			{
				dTheta = dT * x;
				// Use the Euler identity exp(j*x) = cos(x) + j*sin(x)
				dA += adX[x] * Math.cos(dTheta);
				dB += adX[x] * Math.sin(dTheta);
			}
			adFTX[k] = Math.sqrt(dA*dA + dB*dB);
		}

		return adFTX;
	}
/*
	void s2f_cut(ImageProcessor ori,int x0,int y0,ImageProcessor cut)
	{
		for(int i=0;i<cut.getWidth();i++)
		{	
			for(int j=0;j<cut.getHeight();j++)
			{	
				cut.putPixelValue(i,j,slope*ori.getPixelValue(x0+i,y0+j)+intercept);
			}
		}
	}
*/
	public boolean powerOf2Size(ImageProcessor ip) 
	{
		/* 
			check whether ip is a square image and if the dimensions are a power of 2 
			otherwise no 2D-FFT possible
		*/
		
		Rectangle r = ip.getRoi();
		return powerOf2(r.width) && r.width==r.height;
	}

	boolean powerOf2(int n) 
	{		
		/* check whether n is a power of 2*/
		int i=2;
		while(i<n) i *= 2;
		return i==n;
	}

	public void fft(ImageProcessor ip, boolean inverse)  
	{
		/* 
			this method has been copied from the FFT-filter plugin available at the ImageJ site)
			but has been modified in the way that it finishes by adding the square of the resulting 
			2D-FFT to an accumulating ipffttot.
		*/
		int  maxN = ip.getWidth();
		int i,j;  
		double factor; 
		ImageProcessor amp;
		makeSinCosTables(maxN);  
		float[]	fht   =   (float[])ip.getPixels();
		rc2DFHT(fht,   inverse,   maxN);  
		{
			ImageProcessor ps =  calculatePowerSpectrum(fht, maxN);
			amp = calculateAmplitude(fht,maxN);
		//	factor=amp.getWidth()*amp.getHeight()*Math.sqrt(2.0);
			for(i=0;i<amp.getWidth();i++)
				for(j=0;j<amp.getHeight();j++)
					ipffttot.putPixelValue(i,j,(amp.getPixelValue(i,j)*amp.getPixelValue(i,j))+ipffttot.getPixelValue(i,j));
		}
	}

	void makeSinCosTables(int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		int n = maxN/4;
		C = new float[n];
		S = new float[n];
		double theta = 0.0;
		double dTheta = 2.0 * Math.PI/maxN;
		for (int i=0; i<n; i++) {
			C[i] = (float)Math.cos(theta);
			S[i] = (float)Math.sin(theta);
			theta += dTheta;
		}
	}
	
	/** Row-column Fast Hartley Transform */
	void rc2DFHT(float[] x, boolean inverse, int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		//IJ.write("FFT: rc2DFHT (row-column Fast Hartley Transform)");
		for (int row=0; row<maxN; row++)
			dfht3(x, row*maxN, inverse, maxN);
		transposeR(x, maxN);
		for (int row=0; row<maxN; row++)		
			dfht3(x, row*maxN, inverse, maxN);
		transposeR(x, maxN);

		int mRow, mCol;
		float A,B,C,D,E;
		for (int row=0; row<maxN/2; row++) 
		{ // Now calculate actual Hartley transform
			for (int col=0; col<maxN/2; col++) 
			{
				mRow = (maxN - row) % maxN;
				mCol = (maxN - col)  % maxN;
				A = x[row * maxN + col];	//  see Bracewell, 'Fast 2D Hartley Transf.' IEEE Procs. 9/86
				B = x[mRow * maxN + col];
				C = x[row * maxN + mCol];
				D = x[mRow * maxN + mCol];
				E = ((A + D) - (B + C)) / 2;
				x[row * maxN + col] = A - E;
				x[mRow * maxN + col] = B + E;
				x[row * maxN + mCol] = C + E;
				x[mRow * maxN + mCol] = D - E;
			}
		}
	}
	
	/* An optimized real FHT */
	void dfht3 (float[] x, int base, boolean inverse, int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		int i, stage, gpNum, gpIndex, gpSize, numGps, Nlog2;
		int bfNum, numBfs;
		int Ad0, Ad1, Ad2, Ad3, Ad4, CSAd;
		float rt1, rt2, rt3, rt4;

		Nlog2 = log2(maxN);
		BitRevRArr(x, base, Nlog2, maxN);	//bitReverse the input array
		gpSize = 2;     //first & second stages - do radix 4 butterflies once thru
		numGps = maxN / 4;
		for (gpNum=0; gpNum<numGps; gpNum++)  {
			Ad1 = gpNum * 4;
			Ad2 = Ad1 + 1;
			Ad3 = Ad1 + gpSize;
			Ad4 = Ad2 + gpSize;
			rt1 = x[base+Ad1] + x[base+Ad2];   // a + b
			rt2 = x[base+Ad1] - x[base+Ad2];   // a - b
			rt3 = x[base+Ad3] + x[base+Ad4];   // c + d
			rt4 = x[base+Ad3] - x[base+Ad4];   // c - d
			x[base+Ad1] = rt1 + rt3;      // a + b + (c + d)
			x[base+Ad2] = rt2 + rt4;      // a - b + (c - d)
			x[base+Ad3] = rt1 - rt3;      // a + b - (c + d)
			x[base+Ad4] = rt2 - rt4;      // a - b - (c - d)
		 }

		if (Nlog2 > 2) {
			 // third + stages computed here
			gpSize = 4;
			numBfs = 2;
			numGps = numGps / 2;
			//IJ.write("FFT: dfht3 "+Nlog2+" "+numGps+" "+numBfs);
			for (stage=2; stage<Nlog2; stage++) {
				for (gpNum=0; gpNum<numGps; gpNum++) {
					Ad0 = gpNum * gpSize * 2;
					Ad1 = Ad0;     // 1st butterfly is different from others - no mults needed
					Ad2 = Ad1 + gpSize;
					Ad3 = Ad1 + gpSize / 2;
					Ad4 = Ad3 + gpSize;
					rt1 = x[base+Ad1];
					x[base+Ad1] = x[base+Ad1] + x[base+Ad2];
					x[base+Ad2] = rt1 - x[base+Ad2];
					rt1 = x[base+Ad3];
					x[base+Ad3] = x[base+Ad3] + x[base+Ad4];
					x[base+Ad4] = rt1 - x[base+Ad4];
					for (bfNum=1; bfNum<numBfs; bfNum++) {
					// subsequent BF's dealt with together
						Ad1 = bfNum + Ad0;
						Ad2 = Ad1 + gpSize;
						Ad3 = gpSize - bfNum + Ad0;
						Ad4 = Ad3 + gpSize;

						CSAd = bfNum * numGps;
						rt1 = x[base+Ad2] * C[CSAd] + x[base+Ad4] * S[CSAd];
						rt2 = x[base+Ad4] * C[CSAd] - x[base+Ad2] * S[CSAd];

						x[base+Ad2] = x[base+Ad1] - rt1;
						x[base+Ad1] = x[base+Ad1] + rt1;
						x[base+Ad4] = x[base+Ad3] + rt2;
						x[base+Ad3] = x[base+Ad3] - rt2;

					} /* end bfNum loop */
				} /* end gpNum loop */
				gpSize *= 2;
				numBfs *= 2;
				numGps = numGps / 2;
			} /* end for all stages */
		} /* end if Nlog2 > 2 */

		if (inverse)  {
			for (i=0; i<maxN; i++)
			x[base+i] = x[base+i] / maxN;
		}
	}

	void transposeR (float[] x, int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		int   r, c;
		float  rTemp;

		for (r=0; r<maxN; r++)  {
			for (c=r; c<maxN; c++) {
				if (r != c)  {
					rTemp = x[r*maxN + c];
					x[r*maxN + c] = x[c*maxN + r];
					x[c*maxN + r] = rTemp;
				}
			}
		}
	}
	
	int log2 (int x) 
	{
		/* 
			... fft-filter plugin
		*/
		int count = 15;
		while (!btst(x, count))
			count--;
		return count;
	}
	
	private boolean btst (int  x, int bit) 
	{
		/* 
			... fft-filter plugin
		*/
		return ((x & (1<<bit)) != 0);
	}

	void BitRevRArr (float[] x, int base, int bitlen, int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		int    l;
		float[] tempArr = new float[maxN];
		for (int i=0; i<maxN; i++)  {
			l = BitRevX (i, bitlen);  //i=1, l=32767, bitlen=15
			tempArr[i] = x[base+l];
		}
		for (int i=0; i<maxN; i++)
			x[base+i] = tempArr[i];
	}

	
	private int BitRevX (int  x, int bitlen) 
	{
		/* 
			... fft-filter plugin
		*/
		int  temp = 0;
		for (int i=0; i<=bitlen; i++)
			if ((x & (1<<i)) !=0)
				temp  |= (1<<(bitlen-i-1));
		return temp & 0x0000ffff;
	}

	private int bset (int x, int bit) 
	{
		/* 
			... fft-filter plugin
		*/
		x |= (1<<bit);
		return x;
	}

	ImageProcessor calculatePowerSpectrum (float[] fht, int maxN) 
	{
		/* 
			... fft-filter plugin
		*/
		int base;
		float  r, scale;
		float min = Float.MAX_VALUE;
  		float max = Float.MIN_VALUE;
   		float[] fps = new float[maxN*maxN];
 		byte[] ps = new byte[maxN*maxN];

  		for (int row=0; row<maxN; row++) {
			FHTps(row, maxN, fht, fps);
			base = row * maxN;
			for (int col=0; col<maxN; col++) {
				r = fps[base+col];
				if (r<min) min = r;
				if (r>max) max = r;
			}
		}

		if (min<1.0)
			min = 0f;
		else
			min = (float)Math.log(min);
		max = (float)Math.log(max);
		scale = (float)(253.0/(max-min));

		for (int row=0; row<maxN; row++) {
			base = row*maxN;
			for (int col=0; col<maxN; col++) {
				r = fps[base+col];
				if (r<1f)
					r = 0f;
				else
					r = (float)Math.log(r);
				ps[base+col] = (byte)(((r-min)*scale+0.5)+1);
			}
		}
		ImageProcessor ip = new ByteProcessor(maxN, maxN, ps, null);
		swapQuadrants(ip);
		return ip;
	}

	/* Power Spectrum of one row from 2D Hartley Transform. */
 	void FHTps(int row, int maxN, float[] fht, float[] ps) 
 	{
 		/* 
			... fft-filter plugin
		*/
		int base = row*maxN;
		int l;
		for (int c=0; c<maxN; c++) {
			l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
			ps[base+c] = (sqr(fht[base+c]) + sqr(fht[l]))/2f;
 		}
	}

	ImageProcessor calculateAmplitude(float[] fht, int maxN) 
	{
   		/* 
			... fft-filter plugin
		*/
		float[] amp = new float[maxN*maxN];
   		for (int row=0; row<maxN; row++) 
   		{
			amplitude(row, maxN, fht, amp);
		}
		ImageProcessor ip = new FloatProcessor(maxN, maxN, amp, null);
		swapQuadrants(ip);
		return ip;
	}

	/** Amplitude of one row from 2D Hartley Transform. */
 	void amplitude(int row, int maxN, float[] fht, float[] amplitude) 
 	{
 		/* 
			... fft-filter plugin
		*/
		int base = row*maxN;
		int l;
		for (int c=0; c<maxN; c++) 
		{
			l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
			amplitude[base+c] = (float)Math.sqrt(sqr(fht[base+c]) + sqr(fht[l]));
 		}
	}

	float sqr(float x) 
	{
		/* 
			... fft-filter plugin
		*/
			return x*x;
	}

	/**	Swap quadrants 1 and 3 and quadrants 2 and 4 so the power 
		spectrum origin is at the center of the image.
		2 1
		3 4
	*/
 	public void swapQuadrants (ImageProcessor ip) 
 	{
 		/* 
			... fft-filter plugin
		*/
		ImageProcessor t1, t2;
		int size = ip.getWidth()/2;
		ip.setRoi(size,0,size,size);
		t1 = ip.crop();
  		ip.setRoi(0,size,size,size);
		t2 = ip.crop();
		ip.insert(t1,0,size);
		ip.insert(t2,size,0);
		ip.setRoi(0,0,size,size);
		t1 = ip.crop();
  		ip.setRoi(size,size,size,size);
		t2 = ip.crop();
		ip.insert(t1,size,size);
		ip.insert(t2,0,0);
	}
	/*
	public void focusLost(FocusEvent e) 
	{
		long N = Math.abs(Math.round(b));		
		IJ.write("N: "+N);
		if (e.getSource() == TFpixsize) 
		{
			double newpixsize=0.1;
			try 
			{
      				newpixsize = parseDouble(TFpixsize.getText());
    			}
    			catch (NumberFormatException nfe) 
    			{
      				TFpixsize.setText(Double.toString(pixsize));
    			}
    			if(!(newpixsize>0.0))
    				TFpixsize.setText(Double.toString(pixsize));
    			pixsize = parseDouble(TFpixsize.getText());
    			subpixsize=pixsize/N;
		//	Lpixsize.setText(Double.toString(pixsize));
			TFsubpixsize.setText(Double.toString(subpixsize));
		//	Lsubpixsize.setText(Double.toString(subpixsize));
			return;
		}
		
		
		if (e.getSource() == TFsubpixsize) 
		{
			double newsubpixsize=0.01;
			try 
			{
				newsubpixsize = parseDouble(TFsubpixsize.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFsubpixsize.setText(Double.toString(subpixsize));
    			}
    			if((!(newsubpixsize>0.0))||(newsubpixsize>pixsize))
    				TFsubpixsize.setText(Double.toString(subpixsize));
    			subpixsize = parseDouble(TFsubpixsize.getText());
    		//	Lsubpixsize.setText(Double.toString(subpixsize));
			return;
		}
		
		if (e.getSource()==TFproflength)
		{
			try 
			{
				nsamp = Integer.parseInt(TFproflength.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFproflength.setText(Integer.toString(2048));
    			}
    			if(!(nsamp>15))
    				TFproflength.setText(Integer.toString(2048));
    			nsamp = Integer.parseInt(TFproflength.getText());
    		//	Lproflength.setText(Integer.toString(nsamp));
			
			return;	
		}
		
		if (e.getSource()==TFcutsize)
		{
			int cutsizenew=0;
			try 
			{
				cutsizenew = Integer.parseInt(TFcutsize.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFcutsize.setText(Integer.toString(cutsize));
    			}
    			if(!(powerOf2(cutsizenew)&&(cutsizenew>31)))
    				TFcutsize.setText(Integer.toString(cutsize));
    			cutsize = Integer.parseInt(TFcutsize.getText());
    		
			
			return;	
		}
		
		if (e.getSource()==TFnnpsleave)
		{
			int nnpsleavenew=1;
			try 
			{
				nnpsleavenew = Integer.parseInt(TFnnpsleave.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFnnpsleave.setText(Integer.toString(nnpsleave));
    			}
    			if(!(nnpsleavenew>-1))
    				TFnnpsleave.setText(Integer.toString(nnpsleave));
    			nnpsleave = Integer.parseInt(TFnnpsleave.getText());
    		
			
			return;	
		}
		if (e.getSource()==TFnnpscount)
		{
			int nnpscountnew=1;
			try 
			{
				nnpscountnew = Integer.parseInt(TFnnpscount.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFnnpscount.setText(Integer.toString(nnpscount));
    			}
    			if(!(nnpscountnew>0))
    				TFnnpscount.setText(Integer.toString(nnpscount));
    			nnpscount = Integer.parseInt(TFnnpscount.getText());
    		
			
			return;	
		}
		
	}
	*/
	public static int NO_WINDOW=0, HAMMING=1, HANN=2, FLATTOP=3; // fourier1D window function types

	/** Calculates the Fourier amplitudes of an array, based on a 1D Fast Hartley Transform.
	* With no Window function, if the array size is a power of 2, the input function
	* should be either periodic or the data at the beginning and end of the array should
	* approach the same value (the periodic continuation should be smooth).
	* With no Window function, if the array size is not a power of 2, the
	* data should decay towards 0 at the beginning and end of the array.
	* For data that do not fulfill these conditions, a window function can be used to
	* avoid artifacts from the edges. See http://en.wikipedia.org/wiki/Window_function.
	*
	* Supported window functions: Hamming, Hann ("raised cosine"), flat-top. Flat-top
	* refers to the HFT70 function in the report cited below, it is named for its
	* response in the frequency domain: a single-frequency sinewave becomes a peak with
	* a short plateau of 3 roughly equal Fourier amplitudes. It is optimized for
	* measuring amplitudes of signals with well-separated sharp frequencies.
	* All window functions will reduce the frequency resolution; this is especially
	* pronounced for the flat-top window.
	*
	* Normalization is done such that the peak height in the Fourier transform
	* (roughly) corresponds to the RMS amplitude of a sinewave (i.e., amplitude/sqrt(2)),
	* and the first Fourier amplitude corresponds to DC component (average value of
	* the data). If the sine frequency falls between two discrete frequencies of the
	* Fourier transform, peak heights can deviate from the true RMS amplitude by up to
	* approx. 36, 18, 15, and 0.1% for no window function, Hamming, Hann and flat-top
	* window functions, respectively.
	* When calculating the power spectrum from the square of the output, note that the
	* result is quantitative only if the input array size is a power of 2; then the
	* spectral density of the power spectrum must be divided by 1.3628 for the Hamming,
	* 1.5 for the Hann, and 3.4129 for the flat-top window.
	*
	* For more details about window functions, see:
	* G. Heinzel, A. Rdiger, and R. Schilling
	* Spectrum and spectral density estimation by the discrete Fourier transform (DFT),
	* including a comprehensive list of window functions and some new flat-top windows.
	* Technical Report, MPI f. Gravitationsphysik, Hannover, 2002; http://edoc.mpg.de/395068
	*
	* @param data Input array; its size need not be a power of 2. The input is not modified..
	* @param windowType may be NO_WINDOW, then the input array is used as it is.
	* Otherwise, it is multiplied by a window function, which can be HAMMING, HANN or
	* FLATTOP.
	* @return Array with the result, i.e., the RMS amplitudes for each frequency.
	* The output array size is half the size of the 2^n-sized array used for the FHT;
	* array element [0]corresponds to frequency zero (the "DC component"). The first
	* nonexisting array element, result[result.length] would correspond to a frequency
	* of 1 cycle per 2 input points, i.e., the Nyquist frequency. In other words, if
	* the spacing of the input data points is dx, results[i] corresponds to a frequency
	* of i/(2*results.length*dx).
	*/
	public float[] fourier1D(float[] data, int windowType) {
		int n = data.length;
		int size = 2;
		while (size<n) size *= 2; // find power of 2 where the data fit
		//IJ.write("dp mode "+mode);
		float[] y = new float[size]; // leave the original data untouched, work on a copy
		System.arraycopy(data, 0, y, 0, n); // pad to 2^n-size
		double sum = 0;
		if (windowType != NO_WINDOW) {
			for (int x=0; x<n; x++) { //calculate non-normalized window function
				double z = (x + 0.5) * (2 * Math.PI / n);
				double w = 0;
				if (windowType == HAMMING)
					w = 0.54 - 0.46 * Math.cos(z);
				else if (windowType == HANN)
					w = 1. - Math.cos(z);
				else if (windowType == FLATTOP)
					w = 1. - 1.90796 * Math.cos(z) + 1.07349 * Math.cos(2*z) - 0.18199 * Math.cos(3*z);
				else
					throw new IllegalArgumentException("Invalid Fourier Window Type");
				y[x] *= w;
				sum += w;
			}
		} else
			sum = n;
		for (int x=0; x<n; x++) //normalize
			y[x] *= (1./sum);
		transform1D(y); //transform
		float[] result = new float[size/2];
		result[0] = (float)Math.sqrt(y[0]*y[0]);
		for (int x=1; x<size/2; x++)
			result[x] = (float)Math.sqrt(y[x]*y[x]+y[size-x]*y[size-x]);
		return result;
	}
	/** Performs an optimized 1D Fast Hartley Transform (FHT) of an array.
	 *  Array size must be a power of 2.
	 *  Note that all amplitudes in the output 'x' are multiplied by the array length.
	 *  Therefore, to get the power spectrum, for 1 <=i < N/2, use
	 *  ps[i] = (x[i]*x[i]+x[maxN-i]*x[maxN-i])/(maxN*maxN), where maxN is the array length.
	 *  To get the real part of the complex FFT, for i=0 use x[0]/maxN,
	 *  and for i>0, use (x[i]+x[maxN-i])/(2*maxN).
	 *  The imaginary part of the complex FFT, with i>0, is given by (x[i]-x[maxN-i])/(2*maxN)
	 *  The coefficients of cosine and sine are like the real and imaginary values above,
	 *  but you have to divide by maxN instead of 2*maxN.
	 */
	public void transform1D(float[] x) {
		int n = x.length;
		dfht3(x, 0, false, n);
	}
	public void actionPerformed(ActionEvent e) 
	{
		int i,j;
	//	float mtf[];
		double mtf[];
		int rowcounter=0;
		int colcounter;
		
		
		if (e.getSource() == TFpixsize) 
		{
			try 
			{
      				pixsize = parseDouble(TFpixsize.getText());
    			}
    			catch (NumberFormatException nfe) 
    			{
      				TFpixsize.setText("0.1");
    			}
    			if(!(pixsize>0.0))
    				TFpixsize.setText("0.1");
    			pixsize = parseDouble(TFpixsize.getText());
    		//	subpixsize=pixsize/10.0;
		//	Lpixsize.setText(Double.toString(pixsize));
		//	TFsubpixsize.setText(Double.toString(subpixsize));
		//	Lsubpixsize.setText(Double.toString(subpixsize));
			return;
		}
		
		/*
		if (e.getSource() == TFsubpixsize) 
		{
			try 
			{
				subpixsize = parseDouble(TFsubpixsize.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFsubpixsize.setText(Double.toString(pixsize/10.0));
    			}
    			if((!(subpixsize>0.0))||(subpixsize>pixsize))
    				TFsubpixsize.setText(Double.toString(pixsize/10.0));
    			subpixsize = parseDouble(TFsubpixsize.getText());
    		//	Lsubpixsize.setText(Double.toString(subpixsize));
			return;
		}
		*/
		
		if (e.getSource()==TFproflength)
		{
			try 
			{
				nsamp = Integer.parseInt(TFproflength.getText());
			}
    			catch (NumberFormatException nfe) 
    			{
    				TFproflength.setText(Integer.toString(2048));
    			}
    			if(!(nsamp>15))
    				TFproflength.setText(Integer.toString(2048));
    			nsamp = Integer.parseInt(TFproflength.getText());
    		//	Lproflength.setText(Integer.toString(nsamp));
			
			return;	
		}
		
		
		if (e.getSource() == Bcopy)
		{
			JTextArea TAcp=new JTextArea();
			for(j=0;j<JTa.getRowCount();j++)
			{
				for(i=0;i<JTa.getColumnCount();i++)
					if(JTa.getValueAt(j,i)==null)
						TAcp.append("\t");
					else
						TAcp.append((String)JTa.getValueAt(j,i)+"\t");
				TAcp.append("\n");
			}
			TAcp.selectAll();
			TAcp.copy();
			return;
		} 
		
		if ((e.getSource()==Bcopycol1)||(e.getSource()==Bcopycol2)||(e.getSource()==Bcopycol3)||(e.getSource()==Bcopycol4)||(e.getSource()==Bcopycol5)||(e.getSource()==Bcopycol6)||(e.getSource()==Bcopycol7)||(e.getSource()==Bcopycol8))
		{
			int col=-1;
			if (e.getSource()==Bcopycol1) col=0;
			if (e.getSource()==Bcopycol2) col=1;
			if (e.getSource()==Bcopycol3) col=2;
			if (e.getSource()==Bcopycol4) col=3;
			if (e.getSource()==Bcopycol5) col=4;
			if (e.getSource()==Bcopycol6) col=5;
			if (e.getSource()==Bcopycol7) col=6;
			if (e.getSource()==Bcopycol8) col=7;
			JTextArea TAcp=new JTextArea();
			for(j=0;j<JTa.getRowCount();j++)
			{
				for(i=col;i<col+1;i++)
					if(JTa.getValueAt(j,i)==null)
						TAcp.append("\n");
					else
						TAcp.append((String)JTa.getValueAt(j,i)+"\n");
			}
			TAcp.selectAll();
			TAcp.copy();
			return;
		}

		imp = WindowManager.getCurrentImage();
		if (imp == null) 
		{
			IJ.error("There is no active image, Please open an image.");
			return;
		}
		ip=imp.getProcessor();
		imamode=0;
		Rectangle r=ip.getRoi();
		if(CBlin.isSelected())
			imamode=1;
		else if (CBsqrt.isSelected())
			imamode=2;
		else if (CBlog.isSelected())
			imamode=3;
		else if (CBproc.isSelected())
			imamode=4;
		Roi ro=imp.getRoi();
		if (ro==null)
		{
			IJ.error("Rectangular ROI needed");
			return;	
		}
		if (ro.getType()!=Roi.RECTANGLE)
		{
			IJ.error("Rectangular ROI needed");
			return;	
		}
		if ((e.getSource() == Bmtfhor)||(e.getSource()==Bmtfvert)) 
		{
	/*	MTF measurement: horizontal or vertical	*/
			ipcut=ip.crop();
			if (e.getSource()==Bmtfhor)
				colcounter=0;
			else
			{
				colcounter=2;
				ipcut=ipcut.rotateRight();
			}
			impcut=new ImagePlus("ipcut",ipcut);
			JTa.setValueAt("Image name",rowcounter,colcounter);
			JTa.setValueAt(imp.getTitle(),rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Image type",rowcounter,colcounter);
			switch (imamode) 
			{
				case 1:	JTa.setValueAt("linear",rowcounter,colcounter+1);
				break;
				case 2:	JTa.setValueAt("square root",rowcounter,colcounter+1);
				break;
				case 3:	JTa.setValueAt("log",rowcounter,colcounter+1);
				break;
				case 4: JTa.setValueAt("processed",rowcounter,colcounter+1);
				break;
			}
			rowcounter++;
			JTa.setValueAt("MTF area dims",rowcounter,colcounter);
			JTa.setValueAt(""+ipcut.getWidth()+"x"+ipcut.getHeight(),rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Pixel Size",rowcounter,colcounter);
			JTa.setValueAt(""+pixsize,rowcounter,colcounter+1);
			rowcounter++;
			//double fact=N;
			edge_finder(ipcut);
			int N =(int )Math.abs(Math.round(b));		
			IJ.write("N: "+N);
			JTa.setValueAt("Edge angle",rowcounter,colcounter);
			JTa.setValueAt(""+Math.abs(Math.atan(1/b)/3.141592*180.0),rowcounter,colcounter+1);

			
			rowcounter++;
			JTa.setValueAt("Cycles * rows per cycle",rowcounter,colcounter);
			JTa.setValueAt(""+Math.abs(ipcut.getHeight()/b)+" * "+b,rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Length subsampled profile",rowcounter,colcounter);
			JTa.setValueAt(""+nsamp,rowcounter,colcounter+1);
			rowcounter++;
			mode1=get_leftmode(ipcut,a,b);
			mode2=get_rightmode(ipcut,a,b);
			IJ.write("mode1: "+mode1);
			IJ.write("mode2: "+mode2);
			double fprof[];
		
		/*
		 -calculate profile as present in image			|
		 -convert to linear data depending on image mode	|
		 -differentiate esf to obtain lsf
		 -normalize lsf (i.e. so that the sum of all lsf data equals 1, as if the profile varied from 0 to 1
		 ======
		 (those actions are performed in "subsample(...)")
		 -if selected, corrected the f-values for the edge angle
		 -draw a linear ROI to show where the edge was found
		*/
	//		subsample(ImageProcessor ip,double a, double b,int n_sub,int mode,int shift)
		// Number of lines (rows) after which the edge has moved by one pixel (column)
		int read = Math.abs((int)Math.round(b));
		
			fprof=new double[nsamp];
			dprof=new double[nsamp];
//KESKMISTAMINE
//Tuleks lisada tingimused, mitut kasutada ja kas a-le liidetakse või lahutatakse
		if (b>0){
					
			double[] LSF0 = subsample(ipcut,a,b,nsamp,imamode,0);
			double[] LSF1 = subsample(ipcut,a+1,b,nsamp,imamode,read);
			double[] LSF2 = subsample(ipcut,a+2,b,nsamp,imamode,2*read);
			double[] LSF3 = subsample(ipcut,a+3,b,nsamp,imamode,3*read);
			double[] LSF4 = subsample(ipcut,a+4,b,nsamp,imamode,4*read);
			double[] LSF  = new double[LSF0.length]; 
				
			for(i=0;i<LSF0.length;i++){
		
				LSF[i]=(LSF0[i]+LSF1[i]+LSF2[i]+LSF3[i]+LSF4[i])/5;
			//	IJ.write(""+LSF[i]);
			}
			dprof=diff(LSF,nsamp,imamode);
		}
		else 
			{
					
			double[] LSF0 = subsample(ipcut,a,b,nsamp,imamode,0);
			double[] LSF1 = subsample(ipcut,a-1,b,nsamp,imamode,read);
			double[] LSF2 = subsample(ipcut,a-2,b,nsamp,imamode,2*read);
			double[] LSF3 = subsample(ipcut,a-3,b,nsamp,imamode,3*read);
			double[] LSF4 = subsample(ipcut,a-4,b,nsamp,imamode,4*read);
			double[] LSF  = new double[LSF0.length]; 
				
			for(i=0;i<LSF0.length;i++){
		
				LSF[i]=(LSF0[i]+LSF1[i]+LSF2[i]+LSF3[i]+LSF4[i])/5;
				IJ.write(""+LSF[i]);
			}
			dprof=diff(LSF,nsamp,imamode);
		}
			
			
			/* //alternative normalisation: based on surface under dprof
			double sumdprof=0.0;
			
			for(i=0;i<nsamp;i++)
				sumdprof+=dprof[i];
			for(i=0;i<nsamp;i++)
				dprof[i]/=sumdprof;
				
			impcut.killProcessor();
			*/
			if (CBangcorr.isSelected()){
				for(i=0;i<nsamp;i++){
					fprof[i]=(double)i*Math.abs((int)Math.round(b))*(1/(nsamp*pixsize*Math.cos(Math.atan(1/b))));
					//IJ.write("i jagatud sellega:"+(double)1/nsamp/(double)subpixsize/(double)Math.cos(Math.atan(1/b)));
					//IJ.write("nsamp:"+nsamp);
					//IJ.write("subpixsize:"+(double)subpixsize);
					//IJ.write("u:"+Math.abs((int)Math.round(b))*(1/(nsamp*pixsize*Math.cos(Math.atan(1/b)))));
				}
			}
			else{
				for(i=0;i<nsamp;i++){
					fprof[i]=(double)i/nsamp/((double)(pixsize/(Math.abs((int)Math.round(b)))));
					//IJ.write("u:"+nsamp/((double)(pixsize/(Math.abs((int)Math.round(b))))));
				}
			}
			rowcounter++;
			if (CBprof2mtf.isSelected())
			{
//MTF
				mtf=AbsDFT(dprof);
				IJ.write("mtf "+mtf.length );
		/*	
			float[] floatArray = new float[dprof.length];
			for (int p = 0 ; p < dprof.length; p++)
			{
  			 floatArray[p] = (float) dprof[p];
  			//IJ.write("dprof "+dprof[p]+"p "+p );
  			 }
				mtf = fourier1D(floatArray,0);
		*/
				if (e.getSource()==Bmtfhor)
				{
					imp.setRoi(new Line(r.x+(int)a,r.y,(int)(r.x+a+r.height/b),r.y+r.height,imp));
					JTa.setValueAt("f",rowcounter,colcounter);
					JTa.setValueAt("MTF hor",rowcounter,colcounter+1);
					rowcounter++;
				}
				else
				{
					imp.setRoi(new Line(r.x,r.y+r.height-(int)a,r.x+r.width,r.y+r.height-(int)(a+r.width/b),imp));
					JTa.setValueAt("f",rowcounter,colcounter);
					JTa.setValueAt("MTF vert",rowcounter,colcounter+1);
					rowcounter++;
				}
			
			
				for(i=0;fprof[i]<(0.5/pixsize);i++)  //while f<fNyquist
				{
					//IJ.write("dprof "+dprof[i]+"i "+i );
					JTa.setValueAt(""+fprof[i],rowcounter,colcounter);
					JTa.setValueAt(""+mtf[i]/mtf[0],rowcounter,colcounter+1);
					rowcounter++;
				}	
				
				while (rowcounter<JTa.getRowCount())
				{
					JTa.setValueAt("",rowcounter,colcounter);
					JTa.setValueAt("",rowcounter,colcounter+1);
					rowcounter++;
				}
			}
			else
			{
				if (e.getSource()==Bmtfhor)
				{
					imp.setRoi(new Line(r.x+(int)a,r.y,(int)(r.x+a+r.height/b),r.y+r.height,imp));
					JTa.setValueAt("x",rowcounter,colcounter);
					JTa.setValueAt("profile hor",rowcounter,colcounter+1);
					rowcounter++;
			IJ.write("r.x "+r.x+","+r.y+" r.height/b "+r.height/b+" imp "+imp+"r.height"+r.height);
				}
				else
				{
					imp.setRoi(new Line(r.x,r.y+r.height-(int)a,r.x+r.width,r.y+r.height-(int)(a+r.width/b),imp));
					JTa.setValueAt("y",rowcounter,colcounter);
					JTa.setValueAt("prof vert",rowcounter,colcounter+1);
					rowcounter++;
				}
			
				for(i=0;i<nsamp;i++)
				{
					JTa.setValueAt(""+i*(pixsize/(int)Math.abs(Math.round(b))),rowcounter,colcounter);
					JTa.setValueAt(""+dprof[i],rowcounter,colcounter+1);
					rowcounter++;
				}	
				while (rowcounter<JTa.getRowCount())
				{
					JTa.setValueAt("",rowcounter,colcounter);
					JTa.setValueAt("",rowcounter,colcounter+1);
					rowcounter++;
				}
			}
		
			ipcut=null;
			impcut=null;
			fprof=null;
			dprof=null;
			mtf=null;
			imp=null;
			ip=null;
			return;
		}
		
		if (e.getSource() == Bnps) 
		{
			ImageStatistics imstat;
			
			double fnyq=0.5/pixsize;
			colcounter=4;
			int iN=(int)Math.floor(r.width/(cutsize/2))-1,jN=(int)Math.floor(r.height/(cutsize/2))-1;
		//	int i0=(r.x+r.width-(iN+1)*(cutsize/2))/2,j0=(r.y+r.height-(jN+1)*(cutsize/2))/2;
			int i0=r.x,j0=r.y;
			int nzones=0;
			boolean inverse;
			ipffttot=new FloatProcessor(cutsize,cutsize);
			for(i=0;i<cutsize;i++)
				for(j=0;j<cutsize;j++)
					ipffttot.putPixelValue(i,j,0.0);
			for(i=0;i<iN;i++)
			{
				IJ.showStatus("NNPS progress:"+i+"/"+iN);	
				for(j=0;j<jN;j++)
				{
					ip.setRoi(i0+i*cutsize/2,j0+j*cutsize/2,cutsize,cutsize);
					ipcut=ip.crop();
					ipcut=ipcut.convertToFloat();
					switch (imamode) 
					{
						case 1:	
							impcut=new ImagePlus("ipcut",ipcut);
							imstat=impcut.getStatistics();
							if ((imstat.max-imstat.min)<imstat.mean)
							{
								IMean=imstat.mean;
								ipcut.add(-IMean);
								ipcut.multiply(1.0/IMean);
								fft(ipcut, false);
								nzones++;
							}
							else
							{
								IJ.write("Rejected zone "+i+","+j+" ("+(i0+i*(cutsize/2))+","+(j0+j*(cutsize/2))+")"+" mean: "+imstat.mean+" range: "+imstat.min+"-"+imstat.max);	
							}
						break;
						case 2:	
							ipcut.sqr();
							impcut=new ImagePlus("ipcut",ipcut);
							imstat=impcut.getStatistics();
							if ((imstat.max-imstat.min)<imstat.mean)
							{
								IMean=imstat.mean;
								ipcut.add(-IMean);
								ipcut.multiply(1.0/IMean);
								fft(ipcut, false);
								nzones++;
							}
							else
							{
								IJ.write("Rejected zone "+i+","+j+" ("+(i0+i*(cutsize/2))+","+(j0+j*(cutsize/2))+")"+" mean: "+imstat.mean+" range: "+imstat.min+"-"+imstat.max);	
							}
						break;
						case 3:	
							for(i=0;i<ipcut.getHeight();i++)
								for(j=0;j<ipcut.getWidth();j++)
									ipcut.putPixelValue(j,i,Math.exp(ipcut.getPixelValue(i,j)));
							impcut=new ImagePlus("ipcut",ipcut);
							
							imstat=impcut.getStatistics();
							if ((imstat.max-imstat.min)<imstat.mean)
							{
								IMean=imstat.mean;
								ipcut.add(-IMean);
								ipcut.multiply(1.0/IMean);
								fft(ipcut, false);
								nzones++;
							}
							else
							{
								IJ.write("rejected zone "+i+","+j+" ("+(i0+i*(cutsize/2))+","+(j0+j*(cutsize/2))+")"+" mean: "+imstat.mean+" range: "+imstat.min+"-"+imstat.max);	
							}
					
						break;
						case 4: 
							impcut=new ImagePlus("ipcut",ipcut);
							imstat=impcut.getStatistics();
							if ((imstat.max-imstat.min)<imstat.mean)
							{
								IMean=imstat.mean;
								ipcut.add(-IMean);
								fft(ipcut, false);
								nzones++;
							}
							else
							{
								IJ.write("rejected zone "+i+","+j+" ("+(i0+i*(cutsize/2))+","+(j0+j*(cutsize/2))+")"+" mean: "+imstat.mean+" range: "+imstat.min+"-"+imstat.max);	
							}
						break;
					}
				}
			}	
			JTa.setColumnSelectionInterval(colcounter,colcounter+2);
			JTa.clearSelection();
			for(i=0;i<cutsize;i++)
				for(j=0;j<cutsize;j++)
					ipffttot.putPixelValue(i,j,ipffttot.getPixelValue(i,j)/(double)nzones);
			JTa.setValueAt("Image name",rowcounter,colcounter);
			JTa.setValueAt(imp.getTitle(),rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Image type",rowcounter,colcounter);
			switch (imamode) 
			{
				case 1:	JTa.setValueAt("linear",rowcounter,colcounter+1);
				break;
				case 2:	JTa.setValueAt("square root",rowcounter,colcounter+1);
				break;
				case 3:	JTa.setValueAt("log",rowcounter,colcounter+1);
				break;
				case 4: JTa.setValueAt("processed",rowcounter,colcounter+1);
				break;
			}
			rowcounter++;
			JTa.setValueAt("NNPS area dimensions",rowcounter,colcounter);
			JTa.setValueAt(""+r.width+"x"+r.height,rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Pixel Size",rowcounter,colcounter);
			JTa.setValueAt(""+pixsize,rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Number of NNPS-zones",rowcounter,colcounter);
			JTa.setValueAt(""+(iN*jN),rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("Number of excluded NNPS-zones",rowcounter,colcounter);
			JTa.setValueAt(""+(iN*jN-nzones),rowcounter,colcounter+1);
			rowcounter++;
			JTa.setValueAt("NNPS fragment-size",rowcounter,colcounter);
			JTa.setValueAt(""+cutsize,rowcounter,colcounter+1);
			rowcounter++;
			rowcounter++;
			JTa.setValueAt("f",rowcounter,colcounter);
			JTa.setValueAt("NNPS hor",rowcounter,colcounter+1);
			JTa.setValueAt("NNPS vert",rowcounter,colcounter+2);
			JTa.setValueAt("NNPS rad av",rowcounter,colcounter+3);
			rowcounter++;
			ipffttot.multiply((double)(pixsize*pixsize)/(double)(cutsize*cutsize*2.0));
			
			// calculate radially averaged NNPS
			
			int l;
			int lmax=cutsize/2;
			int h=cutsize;
			int b=cutsize;
			
			int x0=lmax;
			int y0=lmax;
			
			double[] accu=new double[lmax];
			int[] count=new int[lmax];
			
			for (i=0;i<lmax;i++)
			{
				accu[i]=0.0;
				count[i]=0;	
			}
			for (j=0;j<h;j++)
				for (i=0;i<b;i++)
				{
					if((Math.abs(i-lmax)>nnpsleave-1)&&(Math.abs(j-lmax)>nnpsleave-1))
					{
						l=(int)Math.sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0));
						if (l<lmax)
						{
							accu[l]+=ipffttot.getPixelValue(i,j);
							count[l]+=1;
						}
					}
				}
			for(i=0;i<lmax;i++)
			{
				
				accu[i]=(accu[i]/count[i]);
			}
						
			
			float nnpsfhor,nnpsfvert;
			for(i=cutsize/2;i<cutsize;i++)
			{
				JTa.setValueAt(""+((double)(i-cutsize/2)/(double)(cutsize/2)*fnyq),rowcounter,colcounter);
				
				nnpsfhor= (float)0.0;
				nnpsfvert= (float)0.0;
				
				for(j=0;j<nnpscount;j++)
				{
//Arvutatakse NPS					
					nnpsfhor  =nnpsfhor+(float) (ipffttot.getPixelValue(i,cutsize/2+(nnpsleave+j))+ipffttot.getPixelValue(i,cutsize/2-(nnpsleave+j)));
					nnpsfvert =nnpsfvert+(float) (ipffttot.getPixelValue(cutsize/2+(nnpsleave+j),i)+ipffttot.getPixelValue(cutsize/2-(nnpsleave+j),i));
				}
				nnpsfhor/=2.0*nnpscount;
				nnpsfvert/=2.0*nnpscount;
				JTa.setValueAt(""+nnpsfhor,rowcounter,colcounter+1);
				JTa.setValueAt(""+nnpsfvert,rowcounter,colcounter+2);
				JTa.setValueAt(""+(float)accu[i-lmax],rowcounter,colcounter+3);
				rowcounter++;
			}	
			while (rowcounter<JTa.getRowCount())
			{
				JTa.setValueAt("",rowcounter,colcounter);
				JTa.setValueAt("",rowcounter,colcounter+1);
				JTa.setValueAt("",rowcounter,colcounter+2);
				JTa.setValueAt("",rowcounter,colcounter+3);
				rowcounter++;
			}
			if (CBshowNNPS.isSelected())
			{
				
				new ImagePlus("NNPS of "+imp.getTitle(),ipffttot).show();
			}
			ipffttot=null;
			ipcut=null;
			impcut=null;
			return;
		}
	}
}
