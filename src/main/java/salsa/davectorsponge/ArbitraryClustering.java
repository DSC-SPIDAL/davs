package salsa.davectorsponge;

import salsa.general.Box;

public class ArbitraryClustering
{
	public int NumberofPoints = 0;
	public int MaxClusterID = 0;
	public int MaxIndependentClusterIndex = 0;
	public String ClusterString = "";

	public int[] PointstoClusterIDs;
	public int[] ClusterIDtoIndex;
	public int[] ClusterIDtoGoldenIndex;
	public int[] ClusterIndextoID;
	public int[] ClusterCountsbyIndex;
	public int[] ClusterGoldenCountsbyIndex;

	public ArbitraryClustering(int NumberofPointsINPUT, String ClusterStringINPUT)
	{
		ClusterString = ClusterStringINPUT;
		NumberofPoints = NumberofPointsINPUT;
		PointstoClusterIDs = new int[NumberofPoints];
	}

	//  Rather Messy Indexing
	//  For each "type" (DAVS Medea Mcluster Golden) there is an array PointstoClusterIDs that maps each Point to a cluster. 
	//  This array is set in public static void ReadLabelsFromFile(string fname) for "foreign types"
	//  This array is set in LCMSCalculateClusterStatus() for "DAVS clusters" and ALL Sponger points are converted to individual clusters with one point each and stored at end of list
	//  PointstoClusterIDs is -1 if no cluster
	//  The Clusters are listed in order they they first appear in file
	//  this.MaxIndependentClusterIndex is number of clusters and index runs from 0 to this.MaxIndependentClusterIndex-1
	//  this.ClusterIndextoID[index] = is RawClusterIndex for this collection. RawClusterIndex will appear in PointstoClusterIDs for each point
	//  this.ClusterCountsbyIndex[index] is number of points in this cluster
	//  this.ClusterGoldenCountsbyIndex[index] is number of Golden points in this cluster
	//
	public final void setup()
	{
		this.MaxClusterID = 0;
		for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
		{
			this.MaxClusterID = Math.max(this.MaxClusterID, this.PointstoClusterIDs[GlobalPointIndex]);
		}
		++this.MaxClusterID;
		this.ClusterIDtoIndex = new int[MaxClusterID];
		this.ClusterIDtoGoldenIndex = new int[MaxClusterID];
		for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
		{
			this.ClusterIDtoIndex[RawClusterIndex] = 0;
			this.ClusterIDtoGoldenIndex[RawClusterIndex] = 0;
		}
		for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
		{
			int RawClusterIndex = this.PointstoClusterIDs[GlobalPointIndex];
			if (RawClusterIndex < 0)
			{
				continue;
			}
			this.ClusterIDtoIndex[RawClusterIndex]++;
			if (Program.GoldenPeaks.PointstoClusterIDs[GlobalPointIndex] >= 0)
			{
				this.ClusterIDtoGoldenIndex[RawClusterIndex]++;
			}
		}
		this.MaxIndependentClusterIndex = 0;
		for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
		{
			int Count = this.ClusterIDtoIndex[RawClusterIndex];
			if (Count == 0)
			{
				this.ClusterIDtoIndex[RawClusterIndex] = -1;
				continue;
			}
			++this.MaxIndependentClusterIndex;
		}

		this.ClusterIndextoID = new int[this.MaxIndependentClusterIndex];
		this.ClusterCountsbyIndex = new int[this.MaxIndependentClusterIndex];
		this.ClusterGoldenCountsbyIndex = new int[this.MaxIndependentClusterIndex];

		int index = 0;
		for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
		{
			int Count = this.ClusterIDtoIndex[RawClusterIndex];
			if (Count < 0)
			{
				continue;
			}
			this.ClusterIndextoID[index] = RawClusterIndex;
			this.ClusterCountsbyIndex[index] = Count;
			this.ClusterGoldenCountsbyIndex[index] = this.ClusterIDtoGoldenIndex[RawClusterIndex];
			this.ClusterIDtoIndex[RawClusterIndex] = index;
			++index;
		}
	} // End setup

	public final void Statistics()
	{ // Individual statistics for clustering

		for (int histogramloops = 0; histogramloops < 2; histogramloops++)
		{
			int[] ListofCounts;
			String HistogramType;
			int CountLimit;
			int[] Limits;
			if (histogramloops == 0)
			{
				ListofCounts = this.ClusterCountsbyIndex;
				CountLimit = 1;
				HistogramType = "Full";
				Limits = new int[] {5, 10, 20, 30};
			}
			else
			{
				ListofCounts = this.ClusterGoldenCountsbyIndex;
				CountLimit = 0;
				HistogramType = "Select Golden Peaks";
				Limits = new int[] {5, 15,40};
			}
			int NumIntervals = Limits.length + 1;
			int MaxCount = 0;
			double PWHammySingle = 0.0;
			double PWHammyReal = 0.0;
			double[] Width_xReal = new double[NumIntervals];
			double[] Width_yReal = new double[NumIntervals];
			int[] NumPeaksinCategory = new int[NumIntervals];
			int[] NumClustersinCategory = new int[NumIntervals];
			for (int Category = 0; Category < NumIntervals; Category++)
			{
				Width_xReal[Category] = 0.0;
				Width_yReal[Category] = 0.0;
				NumPeaksinCategory[Category] = 0;
				NumClustersinCategory[Category] = 0;
			}
			int NSingle = 0;
			int NReal = 0;

			for (int ClusterIndex = 0; ClusterIndex < this.MaxIndependentClusterIndex; ClusterIndex++)
			{
				MaxCount = Math.max(MaxCount, ListofCounts[ClusterIndex]);
			}
			int[] HistCounts = new int[MaxCount + 1];

			int GreaterThanOne = 0;
			double Average = 0.0;
			for (int ClusterIndex = 0; ClusterIndex < this.MaxIndependentClusterIndex; ClusterIndex++)
			{
				int Count = ListofCounts[ClusterIndex];
				if (Count < CountLimit)
				{
					DAVectorUtility.printAndThrowRuntimeException("Illegal Cluster Count " + this.ClusterString + " Index " + ClusterIndex + " Original Label " + this.ClusterIndextoID[ClusterIndex]);
				}
				HistCounts[Count]++;
				if (Count <= 0)
				{
					continue;
				}

				int Category = NumIntervals - 1;
				for (int FindCategory = 0; FindCategory < (NumIntervals - 1); FindCategory++)
				{
					if (Count < Limits[FindCategory])
					{
						Category = FindCategory;
						break;
					}
				}

				int RawCluster = this.ClusterIndextoID[ClusterIndex];
				if (Count > 1)
				{ // Real Clusters
					GreaterThanOne++;
					Average += Count;

					double LocalWidth_x = 0.0;
					double LocalWidth_y = 0.0;
					Box<Double> tempRef_LocalWidth_x = new Box<>(LocalWidth_x);
					Box<Double> tempRef_LocalWidth_y = new Box<>(LocalWidth_y);
					SetWidths(RawCluster, tempRef_LocalWidth_x, tempRef_LocalWidth_y);
					LocalWidth_x = tempRef_LocalWidth_x.content;
					LocalWidth_y = tempRef_LocalWidth_y.content;
					NReal += Count;
					Width_xReal[Category] += LocalWidth_x;
					Width_yReal[Category] += LocalWidth_y;
					NumPeaksinCategory[Category] += Count;
					NumClustersinCategory[Category]++;
					PWHammyReal += 0.5 * (LocalWidth_x + LocalWidth_y);
				}
				else
				{ // Singletons
					++NSingle;
					PWHammySingle += 0.5 * Program.SpongeFactor * Program.SpongeFactor;
				}

			}
			if (GreaterThanOne > 0)
			{
				Average = Average / GreaterThanOne;
			}

			double PWHammy = PWHammyReal + PWHammySingle;

			double Width_xReal_Total = 0.0;
			double Width_yReal_Total = 0.0;
			for (int Category = 0; Category < NumIntervals; Category++)
			{
				Width_xReal_Total += Width_xReal[Category];
				Width_yReal_Total += Width_yReal[Category];
				if (NumPeaksinCategory[Category] == 0)
				{
					continue;
				}
				Width_xReal[Category] = Width_xReal[Category] / NumPeaksinCategory[Category];
				Width_yReal[Category] = Width_yReal[Category] / NumPeaksinCategory[Category];
			}
			if (NReal > 0)
			{
				Width_xReal_Total = Width_xReal_Total / NReal;
				Width_yReal_Total = Width_yReal_Total / NReal;
			}

			String nextline = "\nBasic Statistics of Cluster " + ClusterString + " " + HistogramType + " Max Count " + MaxCount + " Avg Count above 1 " + String.format("%1$4.3f", Average) + " Number of Clusters " + this.MaxIndependentClusterIndex + " Above Count of One " + GreaterThanOne + "\n" + " Total Hamiltonian " + String.format("%1$5.4E", PWHammy) + " Singletons " + NSingle + " Single Hamiltonian " + String.format("%1$5.4E", PWHammySingle) + " > 1 Cluster Particles " + NReal + " >1 Hamiltonian " + String.format("%1$5.4E", PWHammyReal) + " Width-x " + String.format("%1$5.4E", Width_xReal_Total) + " Width-y " + String.format("%1$5.4E", Width_yReal_Total) + " Followed by Count Interval Averages preceded by #Clusters(#Peaks)\n";
			for (int Category = 0; Category < NumIntervals; Category++)
			{
				String start = "Rest ";
				if (Category < (NumIntervals - 1))
				{
					start = "Up to " + Limits[Category] + " : ";
				}
				nextline += start + NumClustersinCategory[Category] + "(" + NumPeaksinCategory[Category] + ") Width-x " + String.format("%1$5.4E", Width_xReal[Category]) + " Width-y " + String.format("%1$5.4E", Width_yReal[Category]) + " * ";
			}
			nextline += "\n";
			int ActualMaxCount = MaxCount + 1;
			Box<Integer> tempRef_ActualMaxCount = new Box<>(ActualMaxCount);
			nextline += TrimHistograms(HistCounts, tempRef_ActualMaxCount, 0, 1);
			ActualMaxCount = tempRef_ActualMaxCount.content;
			for (int histloop = 0; histloop < ActualMaxCount; histloop++)
			{
				nextline += HistCounts[histloop] + ", ";
			}
			DAVectorUtility.SALSAPrint(0, nextline);
		}

	} // End Statistics

	public static String TrimHistograms(int[] CountList, Box<Integer> NumCounts, double start, double interval)
	{ // Remove zeros at end of Histograms

		String summary = "";
		int LocalSize = NumCounts.content;
		for (int CountBackwards = (LocalSize - 1); CountBackwards >= 0; CountBackwards--)
		{
			if (CountList[CountBackwards] > 0)
			{
				break;
			}
			--NumCounts.content;
		}
		if (NumCounts.content == 0)
		{
			summary = "Empty Histogram ";
			NumCounts.content = 1;
			return summary;
		}
		summary = "Start " + String.format("%1$5.4f", start) + " Interval " + String.format("%1$5.4f", interval) + " End " + String.format("%1$5.4f", start + (NumCounts.content - 1) * interval) + " * ";
		return summary;

	} // End TrimHistograms(int[] CountList, ref int NumCounts, double start, double interval)

	public static String TrimHistograms(int[] CountList, Box<Integer> NumCounts, int start, int interval)
	{ // Remove zeros at end of Histograms

		String summary = "";
		int LocalSize = NumCounts.content;
		for (int CountBackwards = (LocalSize - 1); CountBackwards >= 0; CountBackwards--)
		{
			if (CountList[CountBackwards] > 0)
			{
				break;
			}
			--NumCounts.content;
		}
		if (NumCounts.content == 0)
		{
			summary = "Empty Histogram ";
			NumCounts.content = 1;
			return summary;
		}
		summary = "Start " + start + " Interval " + interval + " End " + (start + (NumCounts.content - 1) * interval) + " * ";
		return summary;

	} // End TrimHistograms(int[] CountList, ref int NumCounts, int start, int interval)

	public final void SetWidths(int ClusterNumber, Box<Double> Width_x, Box<Double> Width_y)
	{ // Calculate Contribution to the x and y widths of this cluster. This is NOT divided by Occupation Count

		Width_x.content = 0.0;
		Width_y.content = 0.0;
		double[] Center = new double[2];
		double[] Sigma = new double[2];
		int PointsinCluster = 0;

		for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
		{
			int ClusterforPoint = this.PointstoClusterIDs[GlobalPointIndex];
			if (ClusterforPoint != ClusterNumber)
			{
				continue;
			}
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				Center[VectorIndex] += GoldenExamination.PeakPosition[GlobalPointIndex][VectorIndex];
			}
			++PointsinCluster;
		}
		if (PointsinCluster == 0)
		{
			return;
		}

		Center[0] = Center[0] / PointsinCluster;
		Center[1] = Center[1] / PointsinCluster;
		Sigma[0] = Program.SigmaVectorParameters_i_[0] * Center[0];
		Sigma[1] = Program.SigmaVectorParameters_i_[1];

		for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
		{
			int ClusterforPoint = this.PointstoClusterIDs[GlobalPointIndex];
			if (ClusterforPoint != ClusterNumber)
			{
				continue;
			}
			double tmp = (GoldenExamination.PeakPosition[GlobalPointIndex][0] - Center[0]) / Sigma[0];
			Width_x.content += tmp * tmp;
			tmp = (GoldenExamination.PeakPosition[GlobalPointIndex][1] - Center[1]) / Sigma[1];
			Width_y.content += tmp * tmp;
		}

	} // End SetWidths

	public final void Difference(ArbitraryClustering BaseClusters)
	{
		int Sizelimit = 0;
		for (int skipcount = 0; skipcount < 6; skipcount++)
		{
			int NumberTargetPeaksinMajorClusters = 0;
			int NumberBasePeaks = 0;
			int NumberBaseClustersUsed = 0;
			if (skipcount == 1)
			{
				Sizelimit = 2;
			}
			if (skipcount == 2)
			{
				Sizelimit = 5;
			}
			if (skipcount == 3)
			{
				Sizelimit = 10;
			}
			if (skipcount == 4)
			{
				Sizelimit = 20;
			}
			if (skipcount == 5)
			{
				Sizelimit = 50;
			}

			// UsedTargetClusters = 0 This cluster has no overlap with Golden;  = 1 has overlap but not largest cluster; = 2 matched
			int[] UsedTargetClusters = new int[this.MaxIndependentClusterIndex];
			for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
			{
				UsedTargetClusters[targetindex] = 0;
			}

			//  Loop over Clusters in Base Cluster set
			for (int BaseClusterLoop = 0; BaseClusterLoop < BaseClusters.MaxIndependentClusterIndex; BaseClusterLoop++)
			{
				if (BaseClusters.ClusterCountsbyIndex[BaseClusterLoop] <= Sizelimit)
				{
					continue;
				}
				++NumberBaseClustersUsed;

				//  UsedTargetClusters1 is number of golden peaks of this cluster in given other cluster 
				int RawClusterLabel = BaseClusters.ClusterIndextoID[BaseClusterLoop];
				int[] UsedTargetClusters1 = new int[this.MaxIndependentClusterIndex];
				for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
				{
					UsedTargetClusters1[targetindex] = 0;
				}
				NumberBasePeaks += BaseClusters.ClusterCountsbyIndex[BaseClusterLoop];

				for (int GlobalPointIndex = 0; GlobalPointIndex < BaseClusters.NumberofPoints; GlobalPointIndex++)
				{
					if (BaseClusters.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
					{
						continue;
					}
					int RawtargetCluster = this.PointstoClusterIDs[GlobalPointIndex];
					++UsedTargetClusters1[this.ClusterIDtoIndex[RawtargetCluster]];
				}

				int MaxMatch = 0;
				int TargetClusterMatch = -1;
				for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
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
					UsedTargetClusters[targetindex] = Math.max(1, UsedTargetClusters[targetindex]);
				}
				NumberTargetPeaksinMajorClusters += MaxMatch;
				UsedTargetClusters[TargetClusterMatch] = 2;
			}
			int NumberNotOne = 0;
			int NumberOne = 0;
			int NumberMajor = 0;
			for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
			{
				if (UsedTargetClusters[targetindex] == 0)
				{
					continue;
				}
				if (UsedTargetClusters[targetindex] == 2)
				{
					++NumberMajor;
				}
				if (UsedTargetClusters[targetindex] == 1)
				{
					if (this.ClusterCountsbyIndex[targetindex] == 1)
					{
						++NumberOne;
					}
					else
					{
						++NumberNotOne;
					}
				}
			}

			DAVectorUtility.SALSAPrint(0, "Base= " + BaseClusters.ClusterString + " Compared to " + this.ClusterString + " Count Limit " + Sizelimit + " Base Peaks " + NumberBasePeaks + " In Major Cluster " + NumberTargetPeaksinMajorClusters + " Base Clusters " + NumberBaseClustersUsed + " " + " Mapped to Major Clusters " + NumberMajor + " # Stray >1 " + NumberNotOne + " # Stray Clusters 1 member " + NumberOne);
		}

	} // end Difference

	public int[][] GeneralClusterOneDHistogram; // [x,y][Bin] Histogram of General Clusters for 1D Projections
	public int[] GeneralClusterTwoDHistogram; // [Bin] Histogram of General Clusters for radial distance

	public final void HistogramPeaks(int ClusterCountCut)
	{
		GeneralClusterOneDHistogram = new int[Program.ParameterVectorDimension][];
		GeneralClusterTwoDHistogram = new int[GoldenExamination.TwoDHistogramSize];
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			GeneralClusterOneDHistogram[VectorIndex] = new int[GoldenExamination.OneDHistogramSize];
		}

		int NumberHistogrammed = 0;
		for (int ThisClusterLoop = 0; ThisClusterLoop < this.MaxIndependentClusterIndex; ThisClusterLoop++)
		{ // Loop over all clusters skipping small ones
			if (this.ClusterCountsbyIndex[ThisClusterLoop] <= ClusterCountCut)
			{
				continue;
			}

			int RawClusterLabel = this.ClusterIndextoID[ThisClusterLoop];
			int TotalThisPeaks = this.ClusterCountsbyIndex[ThisClusterLoop];
			++NumberHistogrammed;

			//  Extract Points for this cluster
			double[][] ThisPeakPositions = new double[this.ClusterCountsbyIndex[ThisClusterLoop]][];
			for (int LocalPointIndex = 0; LocalPointIndex < this.ClusterCountsbyIndex[ThisClusterLoop]; LocalPointIndex++)
			{
				ThisPeakPositions[LocalPointIndex] = new double[Program.ParameterVectorDimension];
			}
			int localcount = 0;
			for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
			{
				if (this.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
				{
					continue;
				}
				ThisPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
				ThisPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
				++localcount;
			}

			//  Set Sigmas
			double[] Sigma = new double[Program.ParameterVectorDimension];
			double[] Center = new double[Program.ParameterVectorDimension];
			Box<double[]> tempRef_Center = new Box<>(Center);
			Box<double[]> tempRef_Sigma = new Box<>(Sigma);
			this.SetCenter_Sigma(ThisPeakPositions, TotalThisPeaks, tempRef_Center, tempRef_Sigma);
			Center = tempRef_Center.content;
			Sigma = tempRef_Sigma.content;

			Box<double[]> tempRef_Sigma2 = new Box<>(Sigma);
			GoldenExamination.SetHistograms(ThisPeakPositions, TotalThisPeaks, Center, tempRef_Sigma2, GeneralClusterOneDHistogram, GeneralClusterTwoDHistogram);
			Sigma = tempRef_Sigma2.content;

		} // End loop over clusters

		DAVectorUtility.SALSAPrint(0, "\n" + this.ClusterString + " Point-Center Histograms for Occupation Count Cuts: 1D Hist Size " + GoldenExamination.OneDHistogramSize + " 1D Hist Interval " + String.format("%1$4.3f", GoldenExamination.OneDHistogramInterval) + " Max Val " + String.format("%1$4.3f", GoldenExamination.OneDHistogramMaxVal) + " 2D Hist Size " + GoldenExamination.TwoDHistogramSize + " 2D Hist Interval " + String.format("%1$4.3f", GoldenExamination.TwoDHistogramInterval) + " Max Val " + String.format("%1$4.3f", GoldenExamination.TwoDHistogramMaxVal));
		String tempstring = this.ClusterString + " " + NumberHistogrammed + " Clusters with occupation count cut Greater Than " + ClusterCountCut + "\n";
		int ActualHistSize = GoldenExamination.TwoDHistogramSize;
		Box<Integer> tempRef_ActualHistSize = new Box<>(ActualHistSize);
		tempstring += " 2D " + this.ClusterString + " Histogram " + TrimHistograms(GeneralClusterTwoDHistogram, tempRef_ActualHistSize, 0.0, GoldenExamination.TwoDHistogramInterval);
		ActualHistSize = tempRef_ActualHistSize.content;
		for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
		{
			tempstring += GeneralClusterTwoDHistogram[HistBin] + ", ";
		}
		DAVectorUtility.SALSAPrint(0, tempstring);
		for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
		{
			tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + this.ClusterString + " " + " Histogram ";
			ActualHistSize = GoldenExamination.OneDHistogramSize;
			Box<Integer> tempRef_ActualHistSize2 = new Box<>(ActualHistSize);
			tempstring += TrimHistograms(GeneralClusterOneDHistogram[VectorIndex], tempRef_ActualHistSize2, 0.0, GoldenExamination.OneDHistogramInterval);
			ActualHistSize = tempRef_ActualHistSize2.content;
			for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
			{
				tempstring += GeneralClusterOneDHistogram[VectorIndex][HistBin] + ", ";
			}
			DAVectorUtility.SALSAPrint(0, tempstring);
		}

	} // End HistogramPeaks(int ClusterCountCut)

	public final void SetCenter_Sigma(double[][] PointPositions, int NumberofPoints, Box<double[]> Center, Box<double[]> Sigma)
	{
		for (int ClusterLoop = 0; ClusterLoop < NumberofPoints; ClusterLoop++)
		{
			for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
			{
				Center.content[VectorIndex] += PointPositions[ClusterLoop][VectorIndex];
			}
		}
		Center.content[0] = Center.content[0] / NumberofPoints;
		Center.content[1] = Center.content[1] / NumberofPoints;
		Sigma.content[0] = Program.SigmaVectorParameters_i_[0] * Center.content[0];
		Sigma.content[1] = Program.SigmaVectorParameters_i_[1];

	} // End SetCenterSigma

} // End ArbitraryClustering