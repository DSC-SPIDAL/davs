package edu.indiana.soic.spidal.davs;

import mpi.MPIException;

public class LCMSAnalyze
{
	public static void InitializeLCMS()
	{
		Program.CalculateCorrelationMatrix = true;
		Program.CalculateIndividualWidths = true;
		Program.Printeigenvectors = true;
		Program.CalculateEigenvaluesfromMatrix = true;

	} // End InitializeLCMS()

	public static void JustAnalyze()
	{
		SpongePointAnalysis();

		//  Set up to read full file
		Program.SelectedInputLabel = -100000000;
		Program.InputFileType = 1;
		Program.RestartTemperature = 0.025;

		Program.CompareSolution = 1;
		Program.ComparisonSelectedInputLabel = 2;
		Program.ComparisonInputFileType = 0;
		//Program.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

		Program.Refinement = false;
		Program.UseSponge = true;

		//  Set Final Sigma[0]
//		Program.SigmaVectorParameters_i_[0] = 0.000217;
		//Program.SigmaMethod = 2;


		//  Please set these
		Program.SelectedInputLabel = -1000000; // Old Data Style
		Program.RestartSelectedInputLabel = -1000000; // Old Data Style

		// Modern Style
		Program.SelectedInputLabel = -100000000; // New Data Style
		Program.RestartSelectedInputLabel = -100000000; // New Data Style

		//  Set Sponge
//		Program.SpongeFactor1 = 12.0;
//		Program.SpongeFactor2 = 2.0;
//		Program.SpongeFactor1 = 12.0;
//		Program.SpongeFactor2 = 2.0;

		//  No Sponge
	   /*
	    Program.SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
	    Program.SpongeTemperature2 = -1.0;
	    Program.NearbySpongePointLimit = 4.0;
	    Program.UseSponge = false;
	    */

	}

	public static void JoinEquilibriate()
	{ // Join 2 Files Together (typically original and analysis of its sponge) and equilibriate
		//  Set Distance Matrix File to be output (ClusterFinal) of Previous job i.e. basic file
		//  Set Label File to be second (sponge analysis) or update

		SpongePointAnalysis();

		//  Set up to read full files
		Program.SelectedInputLabel = -100000000;
		Program.InputFileType = 1;
		Program.RestartSelectedInputLabel = -100000000;
		Program.RestartInputFileType = 1;
		Program.RestartClusterFile = Program.config.DistanceMatrixFile;

		Program.Refinement = true;
		Program.RestartTemperature = 0.025; // was 0.2
		Program.Tminimum = 0.005;


		Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);
		Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);

		//  Set Final Sigma[0] for low temperature restart
		Program.SigmaVectorParameters_i_[0] = 0.00000598;
		Program.SigmaMethod = 2;

		Program.CompareSolution = 1;
		Program.ComparisonSelectedInputLabel = 2;
		Program.ComparisonInputFileType = 0;
		Program.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

		//  Please set these
		Program.UseSponge = true;
		Program.SpongeFactor1 = 2.0;
		Program.SpongeFactor2 = 2.0;

		// Old Sponge 1
		Program.SpongeFactor1 = 1.0;
		Program.SpongeFactor2 = 1.0;
		Program.RestartTemperature = 1.0; // was 0.025
		Program.SelectedInputLabel = -1000000;
		Program.RestartSelectedInputLabel = -1000000;

		// Redo DA2D

		//  No Sponge
		Program.SpongeTemperature1 = -1.0; // Minimum Temperature where Sponge Introduced
		Program.SpongeTemperature2 = -1.0;
		Program.NearbySpongePointLimit = 4.0;
		Program.UseSponge = false;
		Program.RestartTemperature = 0.4;
		Program.SelectedInputLabel = -100000000;
		Program.RestartSelectedInputLabel = -100000000;

	} // End JoinEquilibriate()

	public static void ChangeSpongeFactor()
	{

		JoinEquilibriate();

		//  Reduce Sponge
		Program.SpongeFactor1 = 2.0;
		Program.SpongeFactor2 = 1.0;
		Program.SpongeTemperature1 = 5.0;
		Program.SpongeTemperature2 = 2.0;
		Program.RestartTemperature = 5.0;

		Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);

	} // End ChangeSpongeFactor()

	public static void SpongePointAnalysis()
	{ // Analysis of Sponge Points from a Previous run
		//  Set Distance Matrix File to be output (ClusterFinal) of Previous job

		FreshFullAnalysis();
		Program.Refinement = true;
		Program.MaxNumberSplitClusters = 32;

		// Cooling
		Program.InitialCoolingFactor1 = 0.9875;
		Program.FineCoolingFactor1 = 0.9975;
		Program.InitialCoolingFactor2 = 0.99375;
		Program.FineCoolingFactor2 = 0.999375;

		Program.InputFileType = 1; // Output File of Previous Job
		Program.SelectedInputLabel = 0; // Select Sponge Points
		Program.RestartTemperature = -1.0;
		Program.CompareSolution = -1;

		//  Analyze Sponge 1
		Program.SpongeFactor2 = 1.0;

	} // End SpongePointAnalysis()

	public static void FreshFullAnalysis()
	{ // Complete Fresh Analysis

		/*//  Center Parameters
		Program.maxNcentperNode = 200000;
		Program.maxNcentTOTAL = 200000;
		Program.targetNcentperPoint = 8;
		Program.maxNcentperPoint = 17;
		Program.targetMinimumNcentperPoint = 2;
		Program.maxNcentCreated = 5 * Program.maxNcentperNode;

		//  Ending Parameters
		Program.Iterationatend = 2000;
		Program.Malpha_MaxChange = 0.005;
		Program.Tminimum = 0.025;

		//  Split Parameters
		Program.MaxNumberSplitClusters = 256;
		Program.ClusterLimitforDistribution = 256;

		Program.Waititerations = 1;
		Program.Waititerations_Converge = 4;

		//  Set Print Options
		Program.PrintInterval = 50;
		Program.ClusterPrintNumber = 5;
		DAVectorUtility.DebugPrintOption = 1;

		//  Magic Temperatures
		Program.MagicTemperatures[0] = -1.0;
		Program.MagicTemperatures[1] = -1.0;
		Program.MagicTemperatures[2] = -1.0;
		Program.MagicTemperatures[3] = -1.0;
		Program.MagicTemperatures[4] = -1.0;

		// Cooling
		Program.InitialCoolingFactor1 = 0.9875;
		Program.FineCoolingFactor1 = 0.9975;
		Program.InitialCoolingFactor2 = 0.99375;
		Program.FineCoolingFactor2 = 0.999375;

		//  Slow Cooling
		Program.CoolingTemperatureSwitch = 30.0;
		Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);
		Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);

		//  Exponential Cuts
		Program.ExpArgumentCut1 = 20.0;
		Program.ExpArgumentCut2 = 40.0;

		// Correct center sizing parameters
		Program.MinimumScaledWidthsquaredtosplit = 0.5;
		Program.ScaledSquaredDistanceatClosenessTest = 0.35;
		Program.MinimumCountforCluster_C_kwithSponge = 1.5;
		Program.MinimumCountforCluster_C_k = 0.5;
		Program.MinimumCountforCluster_Points = -1;

		Program.TemperatureforClosenessTest = 4.0;
		Program.TemperatureforClosenessTest = 9.0;

		// Set Data Reading parameters including Charge
		//  DAVectorSpongeSection.DistanceMatrixFile is basic
		//  DAVectorSpongeSection.LabelFile is update
		Program.SelectedInputLabel = -60;
		Program.Replicate = 1;
		Program.ComparisonInputFileType = 0;
		Program.ComparisonClusterFile = "E:\\Sali\\InCloud\\IUBox\\My Box Files\\Sponge\\input\\DarTBFULLSorted" +
                ".txt";
		Program.RestartSelectedInputLabel = -100000000;
		Program.RestartInputFileType = 1;
		Program.RestartClusterFile = Program.config.DistanceMatrixFile;
		Program.ComparisonSelectedInputLabel = 2;
		Program.SelectedInputLabel = 2; // Charge 2
		Program.InputFileType = 0; //  Original Harvard format
		Program.RestartTemperature = -1.0;
		Program.CompareSolution = -1;

		//  Sigma Annealing
		Program.SigmaVectorParameters_i_[0] = 0.000598;
		Program.SigmaVectorParameters_i_[1] = 2.35;
		//       Program.SigmaMethod = 2;
		Program.SigmaMethod = 3;
		Program.FinalTargetSigma0 = 0.00000598;
		Program.FinalTargetTemperature = 12.0;
		Program.FinalTargetTemperature = 30.0;

		//  Sponge Parameters
		Program.UseSponge = false;
		Program.SpongePWeight = 0.1;
		Program.CreateSpongeScaledSquaredWidth = 10.0;

		Program.SpongeFactor1 = 12.0;
		Program.SpongeFactor1 = 45.0;
		Program.SpongeFactor2 = 2.0;
		Program.SpongeTemperature1 = 9.0;
		Program.SpongeTemperature1 = 30.0;
		Program.SpongeTemperature2 = 2.5;
		Program.SpongeTemperature2 = 2.0;

		// No Sponge
		Program.CreateSpongeScaledSquaredWidth = -1.0; // Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches this
		Program.SpongeTemperature1 = -1.0; // Minimum Temperature where Sponge Introduced
		Program.SpongeTemperature2 = -1.0;

		// Run F7
		Program.SigmaVectorParameters_i_[0] = 0.000598;
		Program.SigmaMethod = 3;
		Program.FinalTargetSigma0 = 0.00000598;
		Program.FinalTargetTemperature = 12.0;

		Program.ScaledSquaredDistanceatClosenessTest = 0.5;
		Program.TemperatureforClosenessTest = 4.0;
		Program.MaxNumberSplitClusters = 256;
		Program.MinimumScaledWidthsquaredtosplit = 1.5;
		Program.ClusterLimitforDistribution = 256;
		Program.Tminimum = 0.1;
		Program.RestartTemperature = -1.0;
		Program.Iterationatend = 400;

		Program.UseSponge = false;
		Program.SpongeFactor1 = 12.0;
		Program.SpongeFactor2 = 3.0;
		Program.SpongeTemperature1 = 9.0;
		Program.SpongeTemperature2 = 3.0;
		Program.SpongePWeight = 0.1;
		Program.CreateSpongeScaledSquaredWidth = 10.0;

		Program.CoolingTemperatureSwitch = 12.0;
		Program.InitialCoolingFactor1 = 0.9875000;
		Program.FineCoolingFactor1 = 0.9975000;
		Program.FineCoolingFactor2 = 0.9993750;
		Program.InitialCoolingFactor2 = 0.9937500;

		Program.Malpha_MaxChange = 0.005;
		Program.Malpha_MaxChange1 = 0.005;
		Program.MinimumCountforCluster_C_kwithSponge = 1.5;
		Program.MinimumCountforCluster_C_k = 0.5;
		Program.MinimumCountforCluster_Points = -1;

		//  F7 PLUS
		*//*
		Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
		Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
		Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
		Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);
		*/

	} // End FreshFullAnalysis()

    // TODO - Test Method to Setup SetupCharge5_From_TempestF7_1_1_8_2June2013
    // The only change is SelectedInputLabel is set to 5 instead of 2
    private static void SetupCharge5_From_TempestF7_1_1_8_2June2013_w_noduplicate_assignments(){
        Program.maxNcentperNode = 200000;
        Program.maxNcentTOTAL = 200000;
        Program.targetNcentperPoint = 8;
        Program.maxNcentperPoint = 17;
        Program.targetMinimumNcentperPoint = 2;
        Program.maxNcentCreated = 5 * Program.maxNcentperNode;
        Program.Waititerations = 1;
        Program.Waititerations_Converge = 4;
        Program.PrintInterval = 50;
        Program.ClusterPrintNumber = 5;
        DAVectorUtility.DebugPrintOption = 1;
        Program.MagicTemperatures[0] = -1.0;
        Program.MagicTemperatures[1] = -1.0;
        Program.MagicTemperatures[2] = -1.0;
        Program.MagicTemperatures[3] = -1.0;
        Program.MagicTemperatures[4] = -1.0;
        Program.ExpArgumentCut1 = 20.0;
        Program.ExpArgumentCut2 = 40.0;
        Program.Replicate = 1;
        Program.ComparisonInputFileType = 0;
        Program.ComparisonClusterFile = "/home/saliya/sali/csharpToJava/sponge/input/DarTBFULLSorted.txt";
        Program.RestartInputFileType = 1;
        Program.RestartClusterFile = Program.config.DistanceMatrixFile;
        Program.ComparisonSelectedInputLabel = 2;
        Program.SelectedInputLabel = 5; // Charge 5
        Program.InputFileType = 0;
        Program.CompareSolution = -1;
        Program.SigmaVectorParameters_i_[1] = 2.35;
        Program.SigmaVectorParameters_i_[0] = 0.000598;
        Program.SigmaMethod = 3;
        Program.FinalTargetSigma0 = 5.98E-06;
        Program.FinalTargetTemperature = 12.0;
        Program.ScaledSquaredDistanceatClosenessTest = 0.5;
        Program.TemperatureforClosenessTest = 4.0;
        Program.MaxNumberSplitClusters = 256;
        Program.MinimumScaledWidthsquaredtosplit = 1.5;
        Program.ClusterLimitforDistribution = 256;
        Program.Tminimum = 0.1;
        Program.RestartTemperature = -1.0;
        Program.Iterationatend = 400;
        Program.UseSponge = false;
        Program.SpongeFactor1 = 12.0;
        Program.SpongeFactor2 = 3.0;
        Program.SpongeTemperature1 = 9.0;
        Program.SpongeTemperature2 = 3.0;
        Program.SpongePWeight = 0.1;
        Program.CreateSpongeScaledSquaredWidth = 10.0;
        Program.CoolingTemperatureSwitch = 12.0;
        Program.InitialCoolingFactor1 = 79.0 / 80.0;
        Program.FineCoolingFactor1 = 399.0 / 400.0;
        Program.FineCoolingFactor2 = 0.999375;
        Program.InitialCoolingFactor2 = 159.0 / 160.0;
        Program.Malpha_MaxChange = 0.005;
        Program.Malpha_MaxChange1 = 0.005;
        Program.MinimumCountforCluster_C_kwithSponge = 1.5;
        Program.MinimumCountforCluster_C_k = 0.5;
        Program.MinimumCountforCluster_Points = -1;
    }

    /*// TODO - Test Method to Setup SetupCharge5_From_TempestF7_1_1_8_2June2013
    // The only change is SelectedInputLabel is set to 5 instead of 2
    private static void SetupCharge5_From_TempestF7_1_1_8_2June2013()
    {
        Program.maxNcentperNode = 200000;
        Program.maxNcentTOTAL = 200000;
        Program.targetNcentperPoint = 8;
        Program.maxNcentperPoint = 17;
        Program.targetMinimumNcentperPoint = 2;

        Program.maxNcentCreated = 5 * Program.maxNcentperNode;
        Program.Iterationatend = 2000;
        Program.Malpha_MaxChange = 0.005;
        Program.Tminimum = 0.025;
        Program.MaxNumberSplitClusters = 256;
        Program.ClusterLimitforDistribution = 256;
        Program.Waititerations = 1;

        Program.Waititerations_Converge = 4;

        Program.PrintInterval = 50;
        Program.ClusterPrintNumber = 5;
        DAVectorUtility.DebugPrintOption = 1;
        Program.MagicTemperatures[0] = -1.0;
        Program.MagicTemperatures[1] = -1.0;
        Program.MagicTemperatures[2] = -1.0;
        Program.MagicTemperatures[3] = -1.0;
        Program.MagicTemperatures[4] = -1.0;
        Program.InitialCoolingFactor1 = 79.0 / 80.0;
        Program.FineCoolingFactor1 = 399.0 / 400.0;
        Program.InitialCoolingFactor2 = 159.0 / 160.0;
        Program.FineCoolingFactor2 = 0.999375;
        Program.CoolingTemperatureSwitch = 30.0;
        Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
        Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
        Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
        Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);
        Program.InitialCoolingFactor1 = Math.sqrt(Program.InitialCoolingFactor1);
        Program.FineCoolingFactor1 = Math.sqrt(Program.FineCoolingFactor1);
        Program.FineCoolingFactor2 = Math.sqrt(Program.FineCoolingFactor2);
        Program.InitialCoolingFactor2 = Math.sqrt(Program.InitialCoolingFactor2);
        Program.ExpArgumentCut1 = 20.0;
        Program.ExpArgumentCut2 = 40.0;
        Program.MinimumScaledWidthsquaredtosplit = 0.5;
        Program.ScaledSquaredDistanceatClosenessTest = 0.35;
        Program.MinimumCountforCluster_C_kwithSponge = 1.5;
        Program.MinimumCountforCluster_C_k = 0.5;
        Program.MinimumCountforCluster_Points = -1;
        Program.TemperatureforClosenessTest = 4.0;
        Program.TemperatureforClosenessTest = 9.0;
        Program.SelectedInputLabel = -60;
        Program.Replicate = 1;
        Program.ComparisonInputFileType = 0;
        Program.ComparisonClusterFile = "/home/saliya/sali/csharpToJava/sponge/input/DarTBFULLSorted.txt";
        Program.RestartSelectedInputLabel = -100000000;
        Program.RestartInputFileType = 1;
        Program.RestartClusterFile = Program.config.DistanceMatrixFile;
        Program.ComparisonSelectedInputLabel = 2;
        Program.SelectedInputLabel = 5; // Charge 5
        Program.InputFileType = 0;
        Program.RestartTemperature = -1.0;
        Program.CompareSolution = -1;
        Program.SigmaVectorParameters_i_[0] = 0.000598;
        Program.SigmaVectorParameters_i_[1] = 2.35;
        Program.SigmaMethod = 3;
        Program.FinalTargetSigma0 = 5.98E-06;
        Program.FinalTargetTemperature = 12.0;
        Program.FinalTargetTemperature = 30.0;
        Program.UseSponge = false;
        Program.SpongePWeight = 0.1;
        Program.CreateSpongeScaledSquaredWidth = 10.0;
        Program.SpongeFactor1 = 12.0;
        Program.SpongeFactor1 = 45.0;
        Program.SpongeFactor2 = 2.0;
        Program.SpongeTemperature1 = 9.0;
        Program.SpongeTemperature1 = 30.0;
        Program.SpongeTemperature2 = 2.5;
        Program.SpongeTemperature2 = 2.0;
        Program.CreateSpongeScaledSquaredWidth = -1.0;
        Program.SpongeTemperature1 = -1.0;
        Program.SpongeTemperature2 = -1.0;
        Program.SigmaVectorParameters_i_[0] = 0.000598;
        Program.SigmaMethod = 3;
        Program.FinalTargetSigma0 = 5.98E-06;
        Program.FinalTargetTemperature = 12.0;
        Program.ScaledSquaredDistanceatClosenessTest = 0.5;
        Program.TemperatureforClosenessTest = 4.0;
        Program.MaxNumberSplitClusters = 256;
        Program.MinimumScaledWidthsquaredtosplit = 1.5;
        Program.ClusterLimitforDistribution = 256;
        Program.Tminimum = 0.1;
        Program.RestartTemperature = -1.0;
        Program.Iterationatend = 400;
        Program.UseSponge = false;
        Program.SpongeFactor1 = 12.0;
        Program.SpongeFactor2 = 3.0;
        Program.SpongeTemperature1 = 9.0;
        Program.SpongeTemperature2 = 3.0;
        Program.SpongePWeight = 0.1;
        Program.CreateSpongeScaledSquaredWidth = 10.0;
        Program.CoolingTemperatureSwitch = 12.0;
        Program.InitialCoolingFactor1 = 79.0 / 80.0;
        Program.FineCoolingFactor1 = 399.0 / 400.0;
        Program.FineCoolingFactor2 = 0.999375;
        Program.InitialCoolingFactor2 = 159.0 / 160.0;
        Program.Malpha_MaxChange = 0.005;
        Program.Malpha_MaxChange1 = 0.005;
        Program.MinimumCountforCluster_C_kwithSponge = 1.5;
        Program.MinimumCountforCluster_C_k = 0.5;
        Program.MinimumCountforCluster_Points = -1;
    }*/



    public static void LCMSCalculateClusterStatus() throws MPIException {
		double cut = Program.NearbySpongePointLimit;
		if (cut < 0.0)
		{
			cut = Program.SpongeFactor;
		}
		if (cut < 0.0)
		{
			cut = 1.0;
		}
		ClusteringSolution.SetGlobalClusterNumbers();
		VectorAnnealIterate.OutputClusteringResults(""); // Set Point -- Cluster links

		Program.ClusterStatus = new ClusterQuality(ClusteringSolution.TotalClusterSummary.NumberofCenters, ClusteringSolution.TotalClusterSummary.SpongeCluster, Program.NumberNearbyClusters, cut);

		Program.ClusterStatus.SetClusterStatistics();
		Program.ClusterStatus.SetPointStatistics();
		Program.ClusterStatus.SetNearbyClusters();
		Program.ClusterStatus.OutputStatus();
		Program.ClusterStatus.ExperimentAnalysis();

		//  Set up Cluster Comparisons
		//  Convert Sponge to singleton clusters
		//  Note FullClusterNumber will end up as 1 more than total of clusters if there is a sponge
		if (Program.CompareSolution <= 0)
		{
			return;
		}
		int FullClusterNumber = ClusteringSolution.TotalClusterSummary.NumberofCenters;
		int SpongeCluster = ClusteringSolution.TotalClusterSummary.SpongeCluster;
		for (int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
		{
			int NewClusterIndex = Program.ClusterAssignments[GlobalPointIndex];
			if (NewClusterIndex == SpongeCluster)
			{
				NewClusterIndex = FullClusterNumber;
				++FullClusterNumber;
			}
			Program.OurClusters.PointstoClusterIDs[GlobalPointIndex] = NewClusterIndex;
		}
		Program.OurClusters.setup();
		LCMSAnalyze.ClusterComparison();

	} // End LCMSCalculateClusterStatus()

	public static void ClusterComparison()
	{
		if (Program.CompareSolution <= 0)
		{
			return;
		}


		//  GoldenExamination of Cuts and Differences
		GoldenExamination.NumberClusteringMethods = 4;

		GoldenExamination.OneDHistogramInterval = 0.1;
		GoldenExamination.TwoDHistogramInterval = 0.1;
		GoldenExamination.OneDHistogramMaxVal = 8.0;
		GoldenExamination.TwoDHistogramMaxVal = 8.0;
		GoldenExamination.OneDHistogramSize = (int)Math.floor((GoldenExamination.OneDHistogramMaxVal + 0.0001) / GoldenExamination.OneDHistogramInterval) + 1;
		GoldenExamination.TwoDHistogramSize = (int)Math.floor((GoldenExamination.TwoDHistogramMaxVal + 0.0001) / GoldenExamination.TwoDHistogramInterval) + 1;

		GoldenExamination.OneD_Center_HistogramInterval = 0.025;
		GoldenExamination.TwoD_Center_HistogramInterval = 0.025;
		GoldenExamination.OneD_Center_HistogramMaxVal = 8.0;
		GoldenExamination.TwoD_Center_HistogramMaxVal = 8.0;
		GoldenExamination.OneD_Center_HistogramSize = (int)Math.floor((GoldenExamination.OneD_Center_HistogramMaxVal + 0.0001) / GoldenExamination.OneD_Center_HistogramInterval) + 1;
		GoldenExamination.TwoD_Center_HistogramSize = (int)Math.floor((GoldenExamination.TwoD_Center_HistogramMaxVal + 0.0001) / GoldenExamination.TwoD_Center_HistogramInterval) + 1;

		if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
		{
			//  Handleloop = 0 Cut on 2D distance
			//  Handleloop = 1 Cut on each 1D distance (Removed)
			for (int handleloop = 0; handleloop < 1; handleloop++)
			{
				for (int minloop = 0; minloop < 3; minloop++)
				{
					GoldenExamination GoldenExaminationContainer = new GoldenExamination(); // This sets Accumulation Record
					GoldenExaminationContainer.HandleGoldenCluster = handleloop;
					GoldenExaminationContainer.MinimumforAveraging = 3;
					if (minloop == 1)
					{
						GoldenExaminationContainer.MinimumforAveraging = 15;
					}
					if (minloop == 2)
					{
						GoldenExaminationContainer.MinimumforAveraging = 40;
					}
					GoldenExaminationContainer.CutValue = 0.7;
					GoldenExaminationContainer.CutMinimumforAveraging = GoldenExaminationContainer.MinimumforAveraging;
					if (GoldenExaminationContainer.CutMinimumforAveraging > 10)
					{
						GoldenExaminationContainer.CutMinimumforAveraging = GoldenExaminationContainer.MinimumforAveraging - 4;
					}


					ArbitraryClustering[] ThreeMethods = new ArbitraryClustering[GoldenExamination.NumberClusteringMethods];
					ThreeMethods[0] = Program.OurClusters;
					ThreeMethods[1] = Program.MedeaClusters;
					ThreeMethods[2] = Program.MclustClusters;
					ThreeMethods[3] = Program.GoldenPeaks;
					GoldenExaminationContainer.GoldenComparisons(Program.GoldenPeaks, ThreeMethods);
					GoldenExaminationContainer.PrintAccumulation();
				}
			}
		}

		//  Point Position Histograms
		GoldenExamination.OneDHistogramInterval = 0.05;
		GoldenExamination.TwoDHistogramInterval = 0.025;
		GoldenExamination.OneDHistogramMaxVal = 8.0;
		GoldenExamination.TwoDHistogramMaxVal = 8.0;
		GoldenExamination.OneDHistogramSize = (int)Math.floor((GoldenExamination.OneDHistogramMaxVal + 0.0001) / GoldenExamination.OneDHistogramInterval) + 1;
		GoldenExamination.TwoDHistogramSize = (int)Math.floor((GoldenExamination.TwoDHistogramMaxVal + 0.0001) / GoldenExamination.TwoDHistogramInterval) + 1;

		int ClusterCountcut = 5;
		for (int minloop = 0; minloop < 4; minloop++)
		{
			if (minloop == 2)
			{
				ClusterCountcut = 20;
			}
			if (minloop == 3)
			{
				ClusterCountcut = 50;
			}

			Program.OurClusters.HistogramPeaks(ClusterCountcut);

			if(Program.MedeaClusters.MaxIndependentClusterIndex > 0){
				Program.MedeaClusters.HistogramPeaks(ClusterCountcut);
			}

			if(Program.MclustClusters.MaxIndependentClusterIndex > 0){
				Program.MclustClusters.HistogramPeaks(ClusterCountcut);
			}

			if(Program.GoldenPeaks.MaxIndependentClusterIndex > 0){
				Program.GoldenPeaks.HistogramPeaks(ClusterCountcut);
			}

			ClusterCountcut += 5;
		}

		//  Comparison of Basic Clusters
		DAVectorUtility.SALSAPrint(0, "\n**************** Statistics of Clustering in each method and selection to Golden Clusters\nWith means versus occupation count and 1D/2D Point-Center Histograms");
		Program.OurClusters.Statistics();

		if(Program.MedeaClusters.MaxIndependentClusterIndex > 0) {
			Program.MedeaClusters.Statistics();
		}

		if(Program.MclustClusters.MaxIndependentClusterIndex > 0) {
			Program.MclustClusters.Statistics();
		}

		if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
		{
			Program.GoldenPeaks.Statistics();
		}

		if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
		{
			DAVectorUtility.SALSAPrint(0, "\n Comparisons with Golden Clusters");
			Program.OurClusters.Difference(Program.GoldenPeaks);
			Program.MedeaClusters.Difference(Program.GoldenPeaks);
			Program.MclustClusters.Difference(Program.GoldenPeaks);
			Program.GoldenPeaks.Difference(Program.GoldenPeaks);
		}

		DAVectorUtility.SALSAPrint(0, "\n Comparisons with DAVector Clusters");

		if(Program.MedeaClusters.MaxIndependentClusterIndex > 0) {
			Program.MedeaClusters.Difference(Program.OurClusters);
		}
		if(Program.MclustClusters.MaxIndependentClusterIndex > 0) {
			Program.MclustClusters.Difference(Program.OurClusters);
		}

		if(Program.MedeaClusters.MaxIndependentClusterIndex > 0) {
			Program.OurClusters.Difference(Program.MedeaClusters);
		}
		if(Program.MclustClusters.MaxIndependentClusterIndex > 0) {
			Program.OurClusters.Difference(Program.MclustClusters);
		}

	}

} // End LCMSAnalyze