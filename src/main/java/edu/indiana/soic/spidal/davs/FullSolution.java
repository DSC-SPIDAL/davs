package edu.indiana.soic.spidal.davs;

import edu.indiana.soic.spidal.general.Box;

public class FullSolution
{
	public int SpongeCluster; // Sponge Cluster Location
	public double[][] CenterPosition; // Center Position
	public int[] CreatedIndex; // Created Index of Cluster
	public int[] CurrentNode; // Current Node holding Cluster
	public int[] OccupationCount; // Number of points in Cluster
	public int NumberofCenters; // Number of Centers
	public int MaximumNumberofCenters; // maximum Number of Centers
	public int IterationSetAt; // Iteration where this instance set at

	public static int PointsinHistogram = 50; // Number of Points in Histograms

	public FullSolution(int MaxCenters)
	{
		SpongeCluster = -1;
		CenterPosition = new double[MaxCenters][];
		for (int CenterIndex = 0; CenterIndex < MaxCenters; CenterIndex++)
		{
			CenterPosition[CenterIndex] = new double[Program.ParameterVectorDimension];
		}
		CreatedIndex = new int[MaxCenters];
		CurrentNode = new int[MaxCenters];
		OccupationCount = new int[MaxCenters];
		NumberofCenters = 0;
		MaximumNumberofCenters = MaxCenters;
		IterationSetAt = -1;
	}


	public final void HistogramClusterProperties()
	{ // Histogram Occupation Counts

		if (DAVectorUtility.MPI_Rank != 0)
		{
			return;
		}
		int maxcount = 0;
		double InterCenterDistcemax = 0.0;
		for (int FullCenterIndex1 = 0; FullCenterIndex1 < this.NumberofCenters; FullCenterIndex1++)
		{
			if (FullCenterIndex1 == SpongeCluster)
			{
				continue;
			}
			maxcount = Math.max(maxcount, this.OccupationCount[FullCenterIndex1]);
			for (int FullCenterIndex2 = FullCenterIndex1; FullCenterIndex2 < this.NumberofCenters; FullCenterIndex2++)
			{
				if (FullCenterIndex2 == SpongeCluster)
				{
					continue;
				}
				InterCenterDistcemax = Math.max(InterCenterDistcemax, DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterPosition[FullCenterIndex1], CenterPosition[FullCenterIndex2]));
			}
		}
		double hist1min = 0.0;
		double hist1max = maxcount;
		double hist2min = 0.0;
		double hist2max = InterCenterDistcemax;
		Box<Double> tempRef_hist1min = new Box<>(hist1min);
		Box<Double> tempRef_hist1max = new Box<>(hist1max);
		DAVectorUtility.SetUpHistogramRange(FullSolution.PointsinHistogram, tempRef_hist1min, tempRef_hist1max);
		hist1min = tempRef_hist1min.content;
		hist1max = tempRef_hist1max.content;
		Box<Double> tempRef_hist2min = new Box<>(hist2min);
		Box<Double> tempRef_hist2max = new Box<>(hist2max);
		DAVectorUtility.SetUpHistogramRange(FullSolution.PointsinHistogram, tempRef_hist2min, tempRef_hist2max);
		hist2min = tempRef_hist2min.content;
		hist2max = tempRef_hist2max.content;
		int [] Hist1 = new int[2 + FullSolution.PointsinHistogram];
		int [] Hist2 = new int[2 + FullSolution.PointsinHistogram];
		for (int FullCenterIndex1 = 0; FullCenterIndex1 < this.NumberofCenters; FullCenterIndex1++)
		{
			if (FullCenterIndex1 == SpongeCluster)
			{
				continue;
			}
			int thiscount = FullSolution.getHistposition((double) this.OccupationCount[FullCenterIndex1], hist1min, hist1max, FullSolution.PointsinHistogram);
			++Hist1[thiscount];
			for (int FullCenterIndex2 = FullCenterIndex1; FullCenterIndex2 < this.NumberofCenters; FullCenterIndex2++)
			{
				if (FullCenterIndex2 == SpongeCluster)
				{
					continue;
				}
				double tmp = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterPosition[FullCenterIndex1], CenterPosition[FullCenterIndex2]);
				int CorrelCount = FullSolution.getHistposition(tmp, hist2min, hist2max, FullSolution.PointsinHistogram);
				++Hist2[CorrelCount];
			}
		}
		DAVectorUtility.SALSAPrint(0, "\nOccupation Count Histogram " + HistOutput(hist1min, hist1max, FullSolution.PointsinHistogram, Hist1));
		DAVectorUtility.SALSAPrint(0, "\nCenter-Center Distance Histogram " + HistOutput(hist2min, hist2max, FullSolution.PointsinHistogram, Hist2));
	}

	public static int getHistposition(double value, double histmin, double histmax, int NumberofPoints)
	{
		if (value < histmin)
		{
			return 0;
		}
		if (value > histmax)
		{
			return NumberofPoints + 1;
		}
		int index = 1 + (int)Math.floor((value - histmin) * NumberofPoints / (histmax - histmin));
		if (index > NumberofPoints)
		{
			index = NumberofPoints;
		}
		return index;
	}

	public static String HistOutput(double Histmin, double Histmax, int NumberofPoints, int[] Count)
	{
		String message = "Min " + String.format("%1$4.3f", Histmin) + " Max " + String.format("%1$4.3f", Histmax) + " Bin Size " + String.format("%1$4.3f", (Histmax - Histmin) / NumberofPoints) + " Num " + NumberofPoints + "\n";
		for (int histloop = 0; histloop < (NumberofPoints + 2); histloop++)
		{
			message += Count[histloop] + " ";
		}
		return message;
	}

} // End FullSolution