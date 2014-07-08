package salsa.davectorsponge;

import salsa.general.Box;

public class GoldenExamination
{
	public static int[] GoldenID; // Integer ID of Golden records; -1 if not a Golden Point
	public static String[] GoldenLabel; // String ID of Golden records; "" if not a Golden Point
	public static double[][] PeakPosition; // [GlobalIndex][x,y] Peak Positions read in
	public static int NumberClusteringMethods; // Number of Clustering Methods
	public static int OneDHistogramSize; // Size of 1D Point-Cut Center Histograms including overflow
	public static int TwoDHistogramSize; // Size of 2D Point-Cut Center Histograms including overflow
	public static double OneDHistogramInterval; // Interval of 1D Point-Cut Center Histograms
	public static double TwoDHistogramInterval; // Interval of 2D Point-Cut Center Histograms
	public static double OneDHistogramMaxVal; // Upper Limit of 1D Point-Cut Center Histograms
	public static double TwoDHistogramMaxVal; // Upper Limit of 1D Point-Cut Center Histograms
	public static String[] CoordLabels = {"m/z", "RT "};

	public static int OneD_Center_HistogramSize; // Size of 1D Method Center - Golden Center Histograms including overflow
	public static int TwoD_Center_HistogramSize; // Size of 2D Method Center - Golden Center Histograms including overflow
	public static double OneD_Center_HistogramInterval; // Interval of 1D Method Center - Golden Center Histograms
	public static double TwoD_Center_HistogramInterval; // Interval of 2D Method Center - Golden Center Histograms
	public static double OneD_Center_HistogramMaxVal; // Upper Limit of 1D Method Center - Golden Center Histograms
	public static double TwoD_Center_HistogramMaxVal; // Upper Limit of 1D Method Center - Golden Center Histograms

	public double CutValue; // Cut on this
	public int MinimumforAveraging; // Minimum count for an average to be performed
	public int CutMinimumforAveraging; // Minimum count for an average to be performed
	public int HandleGoldenCluster = 0; // How to set target Center =0 Average of all peaks
	public int NumberGoldenAccumulations = 0;
	public int NumberGoldenAccumulationsCenterTotalCalc = 0;
	public int NumberGoldenAccumulationsCenterCutCalc = 0;
	public int[] NumberOtherAccumulationsCenterTotalCalc;
	public int[] NumberOtherAccumulationsCenterCutCalc;
	public int[][][] Other_CenterOneDHistogram; // [Method][x,y][Bin] Histogram of Centers for 1D Projections
	public int[][] Other_CenterTwoDHistogram; // [Method][Bin] Histogram of Centers for radial distance

	public GoldenClusterRecord Accumulation;

	public GoldenExamination()
	{
		NumberGoldenAccumulations = 0;
		NumberGoldenAccumulationsCenterTotalCalc = 0;
		NumberGoldenAccumulationsCenterCutCalc = 0;
		NumberOtherAccumulationsCenterTotalCalc = new int[NumberClusteringMethods];
		NumberOtherAccumulationsCenterCutCalc = new int[NumberClusteringMethods];
		for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
		{
			NumberOtherAccumulationsCenterTotalCalc[ClusteringMethod] = 0;
			NumberOtherAccumulationsCenterCutCalc[ClusteringMethod] = 0;
		}
		Accumulation = new GoldenClusterRecord();
		Accumulation.Label = "Accumulation";
		Accumulation.ID = -1;
		Accumulation.TotalGoldenPeaks = 0;
		Accumulation.TotalGoldenPeaksinCut = 0;
		Accumulation.TotalGoldenPeaksoutofCut = 0;

		Other_CenterOneDHistogram = new int[NumberClusteringMethods][][];
		Other_CenterTwoDHistogram = new int[NumberClusteringMethods][];

		for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
		{
			Other_CenterOneDHistogram[ClusteringMethod] = new int[DAVectorSponge.ParameterVectorDimension][];
			Other_CenterTwoDHistogram[ClusteringMethod] = new int[TwoD_Center_HistogramSize];
			for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
			{
				Other_CenterOneDHistogram[ClusteringMethod][VectorIndex] = new int[OneD_Center_HistogramSize];
			}
		}
	}

	public static class GoldenClusterRecord
	{
		public String Label;
		public int ID;
		public int TotalGoldenPeaks; // Number of Golden peaks input for this cluster
		public int TotalGoldenPeaksinCut; // Number of Golden Peaks in Cut
		public int TotalGoldenPeaksoutofCut; // Number of Golden Peaks out of cut
		public int[] Other_TotalGoldenPeaks; // [Method] Number of Golden peaks for each Method for this cluster
		public int[] Other_TotalGoldenPeaksinCut; // [Method] Number of Golden Peaks for each Method in Cut
		public int[] Other_TotalGoldenPeaksoutofCut; // [Method] Number of Golden Peaks for each Method out of cut but in other cluster
		public int[] Other_TotalGoldenPeaksoutofCluster; // [Method] Number of Golden Peaks for each Method out of cluster
		public int[] Other_TotalNonGoldenPeaksinCluster; // [Method] Number of Non Golden Peaks for each Method in cluster
		public int[] Other_TotalNonGoldenPeaksinCut; // [Method] Number of Non Golden Peaks for each Method in cut
		public double[] GoldenClusterSigma; // [x,y] Standard Deviation -- NOT squared
		public double[] GoldenClusterCenterTotal; // [x,y] Golden Center averaging all peaks
		public double[] GoldenClusterCenterCut; // [x,y] Golden Center averaging peaks in Cut
		public double[][] Other_ClusterCenterDifferenceTotal; // [Method][x,y] Center for Method in Golden Peaks averaging all peaks -- normalized ABSSOLUTE difference
		public double[][] Other_ClusterCenterDifferenceCut; // [Method][x,y] Center Method in Golden Peaks averaging peaks in Cut -- normalized ABSOLUTE difference
		public int[][] GoldenClusterOneDHistogram; // [x,y][Bin] Histogram of Golden Clusters for 1D Projections
		public int[] GoldenClusterTwoDHistogram; // [Bin] Histogram of Golden Clusters for radial distance
		public int[][][] Other_ClusterOneDHistogram; // [Method][x,y][Bin] Histogram of Golden Clusters for 1D Projections
		public int[][] Other_ClusterTwoDHistogram; // [Method][Bin] Histogram of Golden Clusters for radial distance


		public GoldenClusterRecord()
		{

			Other_TotalGoldenPeaks = new int[NumberClusteringMethods];
			Other_TotalGoldenPeaksinCut = new int[NumberClusteringMethods];
			Other_TotalGoldenPeaksoutofCut = new int[NumberClusteringMethods];
			Other_TotalGoldenPeaksoutofCluster = new int[NumberClusteringMethods];
			Other_TotalNonGoldenPeaksinCluster = new int[NumberClusteringMethods];
			Other_TotalNonGoldenPeaksinCut = new int[NumberClusteringMethods];

			// Next 3 not used in an accumulation
			GoldenClusterSigma = new double[DAVectorSponge.ParameterVectorDimension];
			GoldenClusterCenterTotal = new double[DAVectorSponge.ParameterVectorDimension];
			GoldenClusterCenterCut = new double[DAVectorSponge.ParameterVectorDimension];

			Other_ClusterCenterDifferenceTotal = new double[NumberClusteringMethods][];
			Other_ClusterCenterDifferenceCut = new double[NumberClusteringMethods][];

			GoldenClusterOneDHistogram = new int[DAVectorSponge.ParameterVectorDimension][];
			GoldenClusterTwoDHistogram = new int[TwoDHistogramSize];
			Other_ClusterOneDHistogram = new int[NumberClusteringMethods][][];
			Other_ClusterTwoDHistogram = new int[NumberClusteringMethods][];

			for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
			{
				Other_ClusterCenterDifferenceTotal[ClusteringMethod] = new double[DAVectorSponge.ParameterVectorDimension];
				Other_ClusterCenterDifferenceCut[ClusteringMethod] = new double[DAVectorSponge.ParameterVectorDimension];
				Other_ClusterOneDHistogram[ClusteringMethod] = new int[DAVectorSponge.ParameterVectorDimension][];
				Other_ClusterTwoDHistogram[ClusteringMethod] = new int[TwoDHistogramSize];
				for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
				{
					Other_ClusterOneDHistogram[ClusteringMethod][VectorIndex] = new int[OneDHistogramSize];
				}
			}
			for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
			{
				GoldenClusterOneDHistogram[VectorIndex] = new int[OneDHistogramSize];
			}

		} // End Initializer for GoldenClusterRecord

	} // End GoldenClusterRecord

	// Look at Individual Golden Clusters and compares with 3 Clustering Methods
	//  Outputs results for each cluster and accumulations (averages)
	public final void GoldenComparisons(ArbitraryClustering GoldenBase, ArbitraryClustering[] Methods)
	{
		for (int BaseClusterLoop = 0; BaseClusterLoop < GoldenBase.MaxIndependentClusterIndex; BaseClusterLoop++)
		{ // Loop over Golden Clusters skipping small ones
			if (GoldenBase.ClusterCountsbyIndex[BaseClusterLoop] <= this.MinimumforAveraging)
			{
				continue;
			}

			GoldenClusterRecord GoldenRecordforthisCluster = new GoldenClusterRecord();
			int RawClusterLabel = GoldenBase.ClusterIndextoID[BaseClusterLoop];
			GoldenRecordforthisCluster.TotalGoldenPeaks = GoldenBase.ClusterCountsbyIndex[BaseClusterLoop];
			++NumberGoldenAccumulations;

			//  Extract Golden Points for this cluster
			int[] ListofGoldenPoints = new int[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]];
			int[] AlwaysGolden = new int[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]];
			double[][] GoldenPeakPositions = new double[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]][];
			for (int LocalPointIndex = 0; LocalPointIndex < GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]; LocalPointIndex++)
			{
				GoldenPeakPositions[LocalPointIndex] = new double[DAVectorSponge.ParameterVectorDimension];
			}
			int localcount = 0;
			for (int GlobalPointIndex = 0; GlobalPointIndex < GoldenBase.NumberofPoints; GlobalPointIndex++)
			{
				if (GoldenBase.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
				{
					continue;
				}
				AlwaysGolden[localcount] = 1;
				ListofGoldenPoints[localcount] = GlobalPointIndex;
				GoldenPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
				GoldenPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
				GoldenRecordforthisCluster.ID = GoldenExamination.GoldenID[GlobalPointIndex]; // Same for all entries
				GoldenRecordforthisCluster.Label = GoldenExamination.GoldenLabel[GlobalPointIndex];
				++localcount;
			}

			//  Set Sigmas
			double[] GoldenSigma = new double[DAVectorSponge.ParameterVectorDimension];
			Box<double[]> tempRef_GoldenSigma = new Box<>(GoldenSigma);
			this.SetSigma(GoldenPeakPositions, GoldenRecordforthisCluster.TotalGoldenPeaks, tempRef_GoldenSigma);
			GoldenSigma = tempRef_GoldenSigma.content;

			//  Set Cut for Golden Peaks
			int dum1 = 0;
			int dum2 = 0;
			Box<double[]> tempRef_GoldenSigma2 = new Box<>(GoldenSigma);
			Box<Integer> tempRef_TotalGoldenPeaksinCut = new Box<>(GoldenRecordforthisCluster.TotalGoldenPeaksinCut);
			Box<Integer> tempRef_TotalGoldenPeaksoutofCut = new Box<>(GoldenRecordforthisCluster.TotalGoldenPeaksoutofCut);
			Box<Integer> tempRef_dum1 = new Box<>(dum1);
			Box<Integer> tempRef_dum2 = new Box<>(dum2);
			Box<double[]> tempRef_GoldenClusterCenterTotal = new Box<>(GoldenRecordforthisCluster.GoldenClusterCenterTotal);
			Box<double[]> tempRef_GoldenClusterCenterCut = new Box<>(GoldenRecordforthisCluster.GoldenClusterCenterCut);
			this.SetCutStatus(GoldenPeakPositions, AlwaysGolden, GoldenRecordforthisCluster.TotalGoldenPeaks, tempRef_GoldenSigma2, tempRef_TotalGoldenPeaksinCut, tempRef_TotalGoldenPeaksoutofCut, tempRef_dum1, tempRef_dum2, tempRef_GoldenClusterCenterTotal, tempRef_GoldenClusterCenterCut);
			GoldenSigma = tempRef_GoldenSigma2.content;
			GoldenRecordforthisCluster.TotalGoldenPeaksinCut = tempRef_TotalGoldenPeaksinCut.content;
			GoldenRecordforthisCluster.TotalGoldenPeaksoutofCut = tempRef_TotalGoldenPeaksoutofCut.content;
			dum1 = tempRef_dum1.content;
			dum2 = tempRef_dum2.content;
			GoldenRecordforthisCluster.GoldenClusterCenterTotal = tempRef_GoldenClusterCenterTotal.content;
			GoldenRecordforthisCluster.GoldenClusterCenterCut = tempRef_GoldenClusterCenterCut.content;

			//  Basic Accumulations
			++NumberGoldenAccumulationsCenterTotalCalc;
			this.Accumulation.TotalGoldenPeaks += GoldenRecordforthisCluster.TotalGoldenPeaks;
			this.Accumulation.TotalGoldenPeaksinCut += GoldenRecordforthisCluster.TotalGoldenPeaksinCut;
			this.Accumulation.TotalGoldenPeaksoutofCut += GoldenRecordforthisCluster.TotalGoldenPeaksoutofCut;

			if (GoldenRecordforthisCluster.TotalGoldenPeaksinCut <= this.CutMinimumforAveraging)
			{
				continue;
			}

			Box<double[]> tempRef_GoldenSigma3 = new Box<>(GoldenSigma);
			SetHistograms(GoldenPeakPositions, GoldenRecordforthisCluster.TotalGoldenPeaks, GoldenRecordforthisCluster.GoldenClusterCenterCut, tempRef_GoldenSigma3, GoldenRecordforthisCluster.GoldenClusterOneDHistogram, GoldenRecordforthisCluster.GoldenClusterTwoDHistogram);
			GoldenSigma = tempRef_GoldenSigma3.content;

			//  Golden Histogram Accumulations
			++NumberGoldenAccumulationsCenterCutCalc;
			for (int Histbin = 0; Histbin < TwoDHistogramSize; Histbin++)
			{
				this.Accumulation.GoldenClusterTwoDHistogram[Histbin] += GoldenRecordforthisCluster.GoldenClusterTwoDHistogram[Histbin];
			}
			for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
			{
				for (int Histbin = 0; Histbin < OneDHistogramSize; Histbin++)
				{
					this.Accumulation.GoldenClusterOneDHistogram[VectorIndex][Histbin] += GoldenRecordforthisCluster.GoldenClusterOneDHistogram[VectorIndex][Histbin];
				}
			}

			//  Analyze Other Clustering Methods
			for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
			{
				ArbitraryClustering ClusterMethod = Methods[MethodLoop];
				int[] UsedTargetClusters1 = new int[ClusterMethod.MaxIndependentClusterIndex];
				for (int targetindex = 0; targetindex < ClusterMethod.MaxIndependentClusterIndex; targetindex++)
				{
					UsedTargetClusters1[targetindex] = 0;
				}
				for (int LocalPointIndex = 0; LocalPointIndex < GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]; LocalPointIndex++)
				{
					int GlobalPointIndex = ListofGoldenPoints[LocalPointIndex];
					int RawtargetCluster = ClusterMethod.PointstoClusterIDs[GlobalPointIndex];
					++UsedTargetClusters1[ClusterMethod.ClusterIDtoIndex[RawtargetCluster]];
				}

				//  MaxMatch is number of Matching peaks
				//  TargetClusterMatch is index of Method Cluster matching this Golden Cluster
				int MaxMatch = 0;
				int TargetClusterMatch = -1;
				for (int targetindex = 0; targetindex < ClusterMethod.MaxIndependentClusterIndex; targetindex++)
				{
					if (UsedTargetClusters1[targetindex] == 0)
					{
						continue;
					}
					if (MaxMatch < UsedTargetClusters1[targetindex])
					{
						MaxMatch = UsedTargetClusters1[targetindex];
						TargetClusterMatch = targetindex;
					}
				}
				if (MaxMatch <= this.CutMinimumforAveraging)
				{
					continue;
				}

				// Set Points of Clustering Method in this Cluster
				int OtherRawClusterLabel = ClusterMethod.ClusterIndextoID[TargetClusterMatch];
				int num = ClusterMethod.ClusterCountsbyIndex[TargetClusterMatch];
				int[] ListofOtherPoints = new int[num];
				int[] GoldenorNot = new int[num];
				double[][] OtherPeakPositions = new double[num][];
				for (int LocalPointIndex = 0; LocalPointIndex < num; LocalPointIndex++)
				{
					OtherPeakPositions[LocalPointIndex] = new double[DAVectorSponge.ParameterVectorDimension];
				}
				localcount = 0;
				int WrongGoldenCluster = 0;
				for (int GlobalPointIndex = 0; GlobalPointIndex < GoldenBase.NumberofPoints; GlobalPointIndex++)
				{
					if (ClusterMethod.PointstoClusterIDs[GlobalPointIndex] != OtherRawClusterLabel)
					{
						continue;
					}
					GoldenorNot[localcount] = 1;
					if (GoldenBase.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
					{
						++WrongGoldenCluster;
						GoldenorNot[localcount] = 0;
					}
					ListofOtherPoints[localcount] = GlobalPointIndex;
					OtherPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
					OtherPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
					++localcount;
				}
				GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCluster[MethodLoop] = WrongGoldenCluster;
				GoldenRecordforthisCluster.Other_TotalGoldenPeaks[MethodLoop] = MaxMatch;
				GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCluster[MethodLoop] = GoldenRecordforthisCluster.TotalGoldenPeaks - MaxMatch;

				//  Set Cut for Other Method
				int Totin = 0;
				int Totout = 0;
				double[] FullCenterfromClustering = new double[DAVectorSponge.ParameterVectorDimension];
				double[] CutCenterfromClustering = new double[DAVectorSponge.ParameterVectorDimension];
				Box<double[]> tempRef_GoldenSigma4 = new Box<>(GoldenSigma);
				Box<Integer> tempRef_Totin = new Box<>(Totin);
				Box<Integer> tempRef_Totout = new Box<>(Totout);
				Box<Integer> tempRef_Object = new Box<>(GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop]);
				Box<Integer> tempRef_Object2 = new Box<>(GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCut[MethodLoop]);
				Box<double[]> tempRef_FullCenterfromClustering = new Box<>(FullCenterfromClustering);
				Box<double[]> tempRef_CutCenterfromClustering = new Box<>(CutCenterfromClustering);
				this.SetCutStatus(OtherPeakPositions, GoldenorNot, num, tempRef_GoldenSigma4, tempRef_Totin, tempRef_Totout, tempRef_Object, tempRef_Object2, tempRef_FullCenterfromClustering, tempRef_CutCenterfromClustering);
				GoldenSigma = tempRef_GoldenSigma4.content;
				Totin = tempRef_Totin.content;
				Totout = tempRef_Totout.content;
				GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop] = tempRef_Object.content;
				GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCut[MethodLoop] = tempRef_Object2.content;
				FullCenterfromClustering = tempRef_FullCenterfromClustering.content;
				CutCenterfromClustering = tempRef_CutCenterfromClustering.content;
				GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCut[MethodLoop] = Totin - GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop];

				//  Accumulate
				++NumberOtherAccumulationsCenterTotalCalc[MethodLoop];

				this.Accumulation.Other_TotalGoldenPeaks[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaks[MethodLoop];
				this.Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop];
				this.Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCut[MethodLoop];
				this.Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCluster[MethodLoop];
				this.Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop] += GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCluster[MethodLoop];
				this.Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCut[MethodLoop];

				if (Totin <= this.CutMinimumforAveraging)
				{
					continue;
				}

				double distance = 0.0;
				int idist;
				for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
				{
					double tmp = (CutCenterfromClustering[VectorIndex] - GoldenRecordforthisCluster.GoldenClusterCenterCut[VectorIndex]) / GoldenSigma[VectorIndex];
					GoldenRecordforthisCluster.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] = Math.abs(tmp);
					distance += tmp * tmp;
					idist = (int)Math.floor(Math.abs(tmp) / GoldenExamination.OneD_Center_HistogramInterval);
					idist = Math.min(idist, GoldenExamination.OneD_Center_HistogramSize - 1);
					++this.Other_CenterOneDHistogram[MethodLoop][VectorIndex][idist];

					double tmp1 = (FullCenterfromClustering[VectorIndex] - GoldenRecordforthisCluster.GoldenClusterCenterTotal[VectorIndex]) / GoldenSigma[VectorIndex];
					GoldenRecordforthisCluster.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] = Math.abs(tmp1);
				}
				idist = (int)Math.floor(Math.abs(distance) / GoldenExamination.TwoD_Center_HistogramInterval);
				idist = Math.min(idist, GoldenExamination.TwoD_Center_HistogramSize - 1);
				++this.Other_CenterTwoDHistogram[MethodLoop][idist];

				Box<double[]> tempRef_GoldenSigma5 = new Box<>(GoldenSigma);
				SetHistograms(OtherPeakPositions, num, CutCenterfromClustering, tempRef_GoldenSigma5, GoldenRecordforthisCluster.Other_ClusterOneDHistogram[MethodLoop], GoldenRecordforthisCluster.Other_ClusterTwoDHistogram[MethodLoop]);
				GoldenSigma = tempRef_GoldenSigma5.content;

				//  Final accumulations
				++NumberOtherAccumulationsCenterCutCalc[MethodLoop];
				for (int VectorIndex = 0; VectorIndex < DAVectorSponge.ParameterVectorDimension; VectorIndex++)
				{
					this.Accumulation.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] += GoldenRecordforthisCluster.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex];
					this.Accumulation.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] += GoldenRecordforthisCluster.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex];
					for (int Histbin = 0; Histbin < OneDHistogramSize; Histbin++)
					{
						this.Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][Histbin] += GoldenRecordforthisCluster.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][Histbin];
					}
				}
				for (int Histbin = 0; Histbin < TwoDHistogramSize; Histbin++)
				{
					this.Accumulation.Other_ClusterTwoDHistogram[MethodLoop][Histbin] += GoldenRecordforthisCluster.Other_ClusterTwoDHistogram[MethodLoop][Histbin];
				}


			} // End loop over MethodLoop

		} // End loop over Golden Clusters

	} // End GoldenComparisons

	public final void SetSigma(double[][] PointPositions, int NumberofPoints, Box<double[]> Sigma)
	{
		double[] center = new double[2];
		for (int ClusterLoop = 0; ClusterLoop < NumberofPoints; ClusterLoop++)
		{
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				center[VectorIndex] += PointPositions[ClusterLoop][VectorIndex];
			}
		}
		center[0] = center[0] / NumberofPoints;
		center[1] = center[1] / NumberofPoints;
		Sigma.content[0] = DAVectorSponge.SigmaVectorParameters_i_[0] * center[0];
		Sigma.content[1] = DAVectorSponge.SigmaVectorParameters_i_[1];

	} // End SetSigma

	//  Analyze a cluster
	//  Goldenstatus =1 if point Golden
	public final void SetCutStatus(double[][] PointPositions, int[] GoldenStatus, int NumberofPoints, Box<double[]> Sigma, Box<Integer> NumberinCut, Box<Integer> NumberoutofCut, Box<Integer> NumberGoldeninCut, Box<Integer> NumberGoldenoutofCut, Box<double[]> FullCenter, Box<double[]> CutCenter)
	{
		NumberoutofCut.content = 0;
		NumberinCut.content = 0;
		NumberGoldeninCut.content = 0;
		NumberGoldenoutofCut.content = 0;

		for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
		{
			FullCenter.content[VectorIndex] = 0.0;
		}
		int TotalGoldeninCluster = 0;
		for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
		{
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				FullCenter.content[VectorIndex] += PointPositions[ClusterPointLoop][VectorIndex];
			}
			if (GoldenStatus[ClusterPointLoop] == 1)
			{
				++TotalGoldeninCluster;
			}
		}
		for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
		{
			FullCenter.content[VectorIndex] = FullCenter.content[VectorIndex] / NumberofPoints;
			CutCenter.content[VectorIndex] = FullCenter.content[VectorIndex];
		}

		// Iterate points in cut
		double[] NewCenter = new double[DAVectorSponge.ParameterVectorDimension];

		for (int IterationLoop = 0; IterationLoop < 10; IterationLoop++)
		{
			double distance;
			NumberinCut.content = 0;
			NumberGoldeninCut.content = 0;
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				NewCenter[VectorIndex] = 0.0;
			}
			for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
			{
				double tmp;
				distance = 0.0;
				int iout = 0;
				for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
				{
					tmp = (PointPositions[ClusterPointLoop][VectorIndex] - CutCenter.content[VectorIndex]) / Sigma.content[VectorIndex];
					if (Math.abs(tmp) > this.CutValue)
					{
						iout = 1;
					}
					distance += tmp * tmp;
				}
				if ((this.HandleGoldenCluster == 1) && (iout == 1))
				{
					continue;
				}
				if ((this.HandleGoldenCluster == 0) && (distance > this.CutValue))
				{
					continue;
				}
				++NumberinCut.content;
				if (GoldenStatus[ClusterPointLoop] == 1)
				{
					++NumberGoldeninCut.content;
				}
				for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
				{
					NewCenter[VectorIndex] += PointPositions[ClusterPointLoop][VectorIndex];
				}
			}
			if (NumberinCut.content <= 0)
			{
				break;
			}
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				CutCenter.content[VectorIndex] = NewCenter[VectorIndex] / NumberinCut.content;
			}
		}
		NumberoutofCut.content = NumberofPoints - NumberinCut.content;
		NumberGoldenoutofCut.content = TotalGoldeninCluster - NumberGoldeninCut.content;

    } // End SetCutStatus

	//  Set Histograms -- each point is put in three histograms
	public static void SetHistograms(double[][] PointPositions, int NumberofPoints, double[] CutCenter, Box<double[]> Sigma, int[][] OneDHistogram, int[] TwoDHistogram)
	{
		for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
		{
			double tmp;
			double distance = 0.0;
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				tmp = (PointPositions[ClusterPointLoop][VectorIndex] - CutCenter[VectorIndex]) / Sigma.content[VectorIndex];
				int idist = (int)Math.floor(Math.abs(tmp) / GoldenExamination.OneDHistogramInterval);
				idist = Math.min(idist, GoldenExamination.OneDHistogramSize - 1);
				++OneDHistogram[VectorIndex][idist];
				distance += tmp * tmp;
			}
			int idist2D = (int)Math.floor(distance / GoldenExamination.TwoDHistogramInterval);
			idist2D = Math.min(idist2D, GoldenExamination.TwoDHistogramSize - 1);
			++TwoDHistogram[idist2D];
		}

	} // End SetHistograms

	public final void PrintAccumulation()
	{
		String[] ClusteringLabels = {"DAVS  ", "Medea ", "Mclust", "Golden"};

		DAVectorUtility.SALSAPrint(0, "\n*********************** Golden Cluster Analysis Cut Value " + String.format("%1$4.3f", CutValue) + " Minimum needed to Average " + MinimumforAveraging + " After Cut " + CutMinimumforAveraging + " Handling Option " + HandleGoldenCluster);

		DAVectorUtility.SALSAPrint(0, "1D Hist Size " + OneDHistogramSize + " 1D Hist Interval " + String.format("%1$4.3f", OneDHistogramInterval) + " Max Val " + String.format("%1$4.3f", OneDHistogramMaxVal) + " 2D Hist Size " + TwoDHistogramSize + " 2D Hist Interval " + String.format("%1$4.3f", TwoDHistogramInterval) + " Max Val " + String.format("%1$4.3f", TwoDHistogramMaxVal));

		DAVectorUtility.SALSAPrint(0, " Records Accumulated " + NumberGoldenAccumulations);
		String tempstring = "";
		for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
		{
			tempstring += " " + ClusteringLabels[MethodLoop] + " " + NumberOtherAccumulationsCenterTotalCalc[MethodLoop];
		}
		DAVectorUtility.SALSAPrint(0, " Records with Good Centers and Accumulated " + NumberGoldenAccumulationsCenterTotalCalc + tempstring);
		tempstring = "";
		for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
		{
			tempstring += " " + ClusteringLabels[MethodLoop] + " " + NumberOtherAccumulationsCenterCutCalc[MethodLoop];
		}
		DAVectorUtility.SALSAPrint(0, " Records with Good Cut Centers and Histograms " + NumberGoldenAccumulationsCenterCutCalc + tempstring);

		double divisor1 = 1.0 / ((double)NumberGoldenAccumulations);
		double tmp1 = Accumulation.TotalGoldenPeaksinCut * divisor1;
		double tmp2 = Accumulation.TotalGoldenPeaksoutofCut * divisor1;
		DAVectorUtility.SALSAPrint(0, " Golden Cut Analysis Cut Value " + String.format("%1$4.3f", CutValue) + " Minimum needed to Average " + MinimumforAveraging + " After Cut " + CutMinimumforAveraging + " Handling Option " + HandleGoldenCluster + " Full Set of Golden Cluster Peaks " + Accumulation.TotalGoldenPeaksinCut + " In Cut " + Accumulation.TotalGoldenPeaksinCut + " Per Cluster " + String.format("%1$5.4f", tmp1) + " Out of Cut " + Accumulation.TotalGoldenPeaksoutofCut + " Per Cluster " + String.format("%1$5.4f", tmp2));
		tempstring = "\n2D Pure Golden Histogram ";
		int ActualHistSize = TwoDHistogramSize;
		Box<Integer> tempRef_ActualHistSize = new Box<>(ActualHistSize);
		tempstring += ArbitraryClustering.TrimHistograms(Accumulation.GoldenClusterTwoDHistogram, tempRef_ActualHistSize, 0.0, TwoDHistogramInterval);
		ActualHistSize = tempRef_ActualHistSize.content;
		for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
		{
			tempstring += Accumulation.GoldenClusterTwoDHistogram[HistBin] + ", ";
		}
		DAVectorUtility.SALSAPrint(0, tempstring);
		for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
		{
			tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " Golden Histogram ";
			ActualHistSize = OneDHistogramSize;
			Box<Integer> tempRef_ActualHistSize2 = new Box<>(ActualHistSize);
			tempstring += ArbitraryClustering.TrimHistograms(Accumulation.GoldenClusterOneDHistogram[VectorIndex], tempRef_ActualHistSize2, 0.0, OneDHistogramInterval);
			ActualHistSize = tempRef_ActualHistSize2.content;
			for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
			{
				tempstring += Accumulation.GoldenClusterOneDHistogram[VectorIndex][HistBin] + ", ";
			}
			DAVectorUtility.SALSAPrint(0, tempstring);
		}

		//  Loop over Methods
		for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
		{
			double divisor2 = 1.0 / ((double)NumberOtherAccumulationsCenterTotalCalc[MethodLoop]);
			tmp1 = NumberOtherAccumulationsCenterTotalCalc[MethodLoop] * divisor1;
			tmp2 = NumberOtherAccumulationsCenterCutCalc[MethodLoop] * divisor1;
			DAVectorUtility.SALSAPrint(0, " Golden Cut Analysis Cut Value " + String.format("%1$4.3f", CutValue) + " Minimum needed to Average " + MinimumforAveraging + " After Cut " + CutMinimumforAveraging + " Handling Option " + HandleGoldenCluster + "\nMethod " + ClusteringLabels[MethodLoop] + " Accumulations " + NumberOtherAccumulationsCenterTotalCalc[MethodLoop] + " " + String.format("%1$5.4f", tmp1) + " with Cut and Histogram " + NumberOtherAccumulationsCenterCutCalc[MethodLoop] + " " + String.format("%1$5.4f", tmp2));
			int totalGPeaksinthislot = Accumulation.Other_TotalGoldenPeaks[MethodLoop] + Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop];
			double divisor3 = 1.0 / ((double)totalGPeaksinthislot);
			tmp1 = Accumulation.Other_TotalGoldenPeaks[MethodLoop] * divisor3;
			tmp2 = Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop] * divisor3;
			double tmp3 = Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop] * divisor3;
			DAVectorUtility.SALSAPrint(0, "Method " + ClusteringLabels[MethodLoop] + " Golden Peaks in Clusters " + Accumulation.Other_TotalGoldenPeaks[MethodLoop] + " " + String.format("%1$5.4f", tmp1) + " Golden Peaks NOT in Method Cluster at all " + Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop] + " " + String.format("%1$5.4f", tmp2) + " False Peaks in Cluster " + Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop] + " " + String.format("%1$5.4f", tmp3) + " NOT in divisor");
			tmp1 = Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop] * divisor3;
			tmp2 = Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop] * divisor3;
			tmp3 = Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop] * divisor3;
			DAVectorUtility.SALSAPrint(0, "Method " + ClusteringLabels[MethodLoop] + " Golden Peaks in Cut " + Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop] + " " + String.format("%1$5.4f", tmp1) + " Golden Peaks out of cut " + Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop] + " " + String.format("%1$5.4f", tmp2) + " False Peaks in Cut " + Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop] + " " + String.format("%1$5.4f", tmp3) + " NOT in divisor");

			double divisor4 = 1.0 / ((double)NumberOtherAccumulationsCenterTotalCalc[MethodLoop]);
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				tmp1 = Accumulation.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] * divisor4;
				tmp2 = Accumulation.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] * divisor4;
				DAVectorUtility.SALSAPrint(0, GoldenExamination.CoordLabels[VectorIndex] + " Method " + ClusteringLabels[MethodLoop] + " Average discrepancy for Full Cluster " + String.format("%1$5.4f", tmp1) + " Average discrepancy for Cut Cluster " + String.format("%1$5.4f", tmp2));
			}

			tempstring = "\n2D " + ClusteringLabels[MethodLoop] + " Abs(Point - Center) Histogram ";
			ActualHistSize = TwoDHistogramSize;
			Box<Integer> tempRef_ActualHistSize3 = new Box<>(ActualHistSize);
			tempstring += ArbitraryClustering.TrimHistograms(Accumulation.Other_ClusterTwoDHistogram[MethodLoop], tempRef_ActualHistSize3, 0.0, TwoDHistogramInterval);
			ActualHistSize = tempRef_ActualHistSize3.content;
			for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
			{
				tempstring += Accumulation.Other_ClusterTwoDHistogram[MethodLoop][HistBin] + ", ";
			}
			DAVectorUtility.SALSAPrint(0, tempstring);
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + ClusteringLabels[MethodLoop] + " " + " Abs(Point - Center) Histogram ";
				ActualHistSize = OneDHistogramSize;
				Box<Integer> tempRef_ActualHistSize4 = new Box<>(ActualHistSize);
				tempstring += ArbitraryClustering.TrimHistograms(Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex], tempRef_ActualHistSize4, 0.0, OneDHistogramInterval);
				ActualHistSize = tempRef_ActualHistSize4.content;
				for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
				{
					tempstring += Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][HistBin] + ", ";
				}
				DAVectorUtility.SALSAPrint(0, tempstring);
			}

			tempstring = "\n2D " + ClusteringLabels[MethodLoop] + " Method-Golden Center Histogram ";
			ActualHistSize = TwoD_Center_HistogramSize;
			Box<Integer> tempRef_ActualHistSize5 = new Box<>(ActualHistSize);
			tempstring += ArbitraryClustering.TrimHistograms(this.Other_CenterTwoDHistogram[MethodLoop], tempRef_ActualHistSize5, 0.0, TwoD_Center_HistogramInterval);
			ActualHistSize = tempRef_ActualHistSize5.content;
			for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
			{
				tempstring += this.Other_CenterTwoDHistogram[MethodLoop][HistBin] + ", ";
			}
			DAVectorUtility.SALSAPrint(0, tempstring);
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + ClusteringLabels[MethodLoop] + " " + " Method-Golden Center Histogram ";
				ActualHistSize = TwoD_Center_HistogramSize;
				Box<Integer> tempRef_ActualHistSize6 = new Box<>(ActualHistSize);
				tempstring += ArbitraryClustering.TrimHistograms(this.Other_CenterOneDHistogram[MethodLoop][VectorIndex], tempRef_ActualHistSize6, 0.0, OneD_Center_HistogramInterval);
				ActualHistSize = tempRef_ActualHistSize6.content;
				for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
				{
					tempstring += this.Other_CenterOneDHistogram[MethodLoop][VectorIndex][HistBin] + ", ";
				}
				DAVectorUtility.SALSAPrint(0, tempstring);
			}

		} // End output Loop over Methods

		DAVectorUtility.SALSAPrint(0, "\n***************** End Average Status of Golden Clusters Cut Value " + String.format("%1$4.3f", CutValue) + " Minimum needed to Average " + MinimumforAveraging + " After Cut " + CutMinimumforAveraging + " Handling Option " + HandleGoldenCluster + "\n");

	} // End PrintAccumulation()

} // End GoldenExamination
 // End Namespace salsa.davectorsponge