package salsa.davectorsponge;

import mpi.MPIException;

public class LCMSAnalyze
{
	public static void InitializeLCMS()
	{
		DAVectorSponge.CalculateCorrelationMatrix = true;
		DAVectorSponge.CalculateIndividualWidths = true;
		DAVectorSponge.Printeigenvectors = true;
		DAVectorSponge.CalculateEigenvaluesfromMatrix = true;

	} // End InitializeLCMS()

	public static void JustAnalyze()
	{
		SpongePointAnalysis();

		//  Set up to read full file
		DAVectorSponge.SelectedInputLabel = -100000000;
		DAVectorSponge.InputFileType = 1;
		DAVectorSponge.RestartTemperature = 0.025;

		DAVectorSponge.CompareSolution = 1;
		DAVectorSponge.ComparisonSelectedInputLabel = 2;
		DAVectorSponge.ComparisonInputFileType = 0;
		DAVectorSponge.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

		DAVectorSponge.Refinement = false;
		DAVectorSponge.UseSponge = true;

		//  Set Final Sigma[0]
		DAVectorSponge.SigmaVectorParameters_i_[0] = 0.00000598;
		DAVectorSponge.SigmaMethod = 2;


		//  Please set these
		DAVectorSponge.SelectedInputLabel = -1000000; // Old Data Style
		DAVectorSponge.RestartSelectedInputLabel = -1000000; // Old Data Style

		// Modern Style
		DAVectorSponge.SelectedInputLabel = -100000000; // New Data Style
		DAVectorSponge.RestartSelectedInputLabel = -100000000; // New Data Style

		//  Set Sponge
		DAVectorSponge.SpongeFactor1 = 2.0;
		DAVectorSponge.SpongeFactor2 = 2.0;
		DAVectorSponge.SpongeFactor1 = 1.0;
		DAVectorSponge.SpongeFactor2 = 1.0;

		//  No Sponge
	   /*
	    DAVectorSponge.SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
	    DAVectorSponge.SpongeTemperature2 = -1.0;
	    DAVectorSponge.NearbySpongePointLimit = 4.0;
	    DAVectorSponge.UseSponge = false;
	    */

	}

	public static void JoinEquilibriate()
	{ // Join 2 Files Together (typically original and analysis of its sponge) and equilibriate
		//  Set Distance Matrix File to be output (ClusterFinal) of Previous job i.e. basic file
		//  Set Label File to be second (sponge analysis) or update

		SpongePointAnalysis();

		//  Set up to read full files
		DAVectorSponge.SelectedInputLabel = -100000000;
		DAVectorSponge.InputFileType = 1;
		DAVectorSponge.RestartSelectedInputLabel = -100000000;
		DAVectorSponge.RestartInputFileType = 1;
		DAVectorSponge.RestartClusterFile = DAVectorSponge.config.DistanceMatrixFile;

		DAVectorSponge.Refinement = true;
		DAVectorSponge.RestartTemperature = 0.025; // was 0.2
		DAVectorSponge.Tminimum = 0.005;


		DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);

		//  Set Final Sigma[0] for low temperature restart
		DAVectorSponge.SigmaVectorParameters_i_[0] = 0.00000598;
		DAVectorSponge.SigmaMethod = 2;

		DAVectorSponge.CompareSolution = 1;
		DAVectorSponge.ComparisonSelectedInputLabel = 2;
		DAVectorSponge.ComparisonInputFileType = 0;
		DAVectorSponge.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

		//  Please set these
		DAVectorSponge.UseSponge = true;
		DAVectorSponge.SpongeFactor1 = 2.0;
		DAVectorSponge.SpongeFactor2 = 2.0;

		// Old Sponge 1
		DAVectorSponge.SpongeFactor1 = 1.0;
		DAVectorSponge.SpongeFactor2 = 1.0;
		DAVectorSponge.RestartTemperature = 1.0; // was 0.025
		DAVectorSponge.SelectedInputLabel = -1000000;
		DAVectorSponge.RestartSelectedInputLabel = -1000000;

		// Redo DA2D

		//  No Sponge
		DAVectorSponge.SpongeTemperature1 = -1.0; // Minimum Temperature where Sponge Introduced
		DAVectorSponge.SpongeTemperature2 = -1.0;
		DAVectorSponge.NearbySpongePointLimit = 4.0;
		DAVectorSponge.UseSponge = false;
		DAVectorSponge.RestartTemperature = 0.4;
		DAVectorSponge.SelectedInputLabel = -100000000;
		DAVectorSponge.RestartSelectedInputLabel = -100000000;

	} // End JoinEquilibriate()

	public static void ChangeSpongeFactor()
	{

		JoinEquilibriate();

		//  Reduce Sponge
		DAVectorSponge.SpongeFactor1 = 2.0;
		DAVectorSponge.SpongeFactor2 = 1.0;
		DAVectorSponge.SpongeTemperature1 = 5.0;
		DAVectorSponge.SpongeTemperature2 = 2.0;
		DAVectorSponge.RestartTemperature = 5.0;

		DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);

	} // End ChangeSpongeFactor()

	public static void SpongePointAnalysis()
	{ // Analysis of Sponge Points from a Previous run
		//  Set Distance Matrix File to be output (ClusterFinal) of Previous job

		FreshFullAnalysis();
		DAVectorSponge.Refinement = true;
		DAVectorSponge.MaxNumberSplitClusters = 32;

		// Cooling
		DAVectorSponge.InitialCoolingFactor1 = 0.9875;
		DAVectorSponge.FineCoolingFactor1 = 0.9975;
		DAVectorSponge.InitialCoolingFactor2 = 0.99375;
		DAVectorSponge.FineCoolingFactor2 = 0.999375;

		DAVectorSponge.InputFileType = 1; // Output File of Previous Job
		DAVectorSponge.SelectedInputLabel = 0; // Select Sponge Points
		DAVectorSponge.RestartTemperature = -1.0;
		DAVectorSponge.CompareSolution = -1;

		//  Analyze Sponge 1
		DAVectorSponge.SpongeFactor2 = 1.0;

	} // End SpongePointAnalysis()

	public static void FreshFullAnalysis()
	{ // Complete Fresh Analysis

		/*//  Center Parameters
		DAVectorSponge.maxNcentperNode = 200000;
		DAVectorSponge.maxNcentTOTAL = 200000;
		DAVectorSponge.targetNcentperPoint = 8;
		DAVectorSponge.maxNcentperPoint = 17;
		DAVectorSponge.targetMinimumNcentperPoint = 2;
		DAVectorSponge.maxNcentCreated = 5 * DAVectorSponge.maxNcentperNode;

		//  Ending Parameters
		DAVectorSponge.Iterationatend = 2000;
		DAVectorSponge.Malpha_MaxChange = 0.005;
		DAVectorSponge.Tminimum = 0.025;

		//  Split Parameters
		DAVectorSponge.MaxNumberSplitClusters = 256;
		DAVectorSponge.ClusterLimitforDistribution = 256;

		DAVectorSponge.Waititerations = 1;
		DAVectorSponge.Waititerations_Converge = 4;

		//  Set Print Options
		DAVectorSponge.PrintInterval = 50;
		DAVectorSponge.ClusterPrintNumber = 5;
		DAVectorUtility.DebugPrintOption = 1;

		//  Magic Temperatures
		DAVectorSponge.MagicTemperatures[0] = -1.0;
		DAVectorSponge.MagicTemperatures[1] = -1.0;
		DAVectorSponge.MagicTemperatures[2] = -1.0;
		DAVectorSponge.MagicTemperatures[3] = -1.0;
		DAVectorSponge.MagicTemperatures[4] = -1.0;

		// Cooling
		DAVectorSponge.InitialCoolingFactor1 = 0.9875;
		DAVectorSponge.FineCoolingFactor1 = 0.9975;
		DAVectorSponge.InitialCoolingFactor2 = 0.99375;
		DAVectorSponge.FineCoolingFactor2 = 0.999375;

		//  Slow Cooling
		DAVectorSponge.CoolingTemperatureSwitch = 30.0;
		DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);

		//  Exponential Cuts
		DAVectorSponge.ExpArgumentCut1 = 20.0;
		DAVectorSponge.ExpArgumentCut2 = 40.0;

		// Correct center sizing parameters
		DAVectorSponge.MinimumScaledWidthsquaredtosplit = 0.5;
		DAVectorSponge.ScaledSquaredDistanceatClosenessTest = 0.35;
		DAVectorSponge.MinimumCountforCluster_C_kwithSponge = 1.5;
		DAVectorSponge.MinimumCountforCluster_C_k = 0.5;
		DAVectorSponge.MinimumCountforCluster_Points = -1;

		DAVectorSponge.TemperatureforClosenessTest = 4.0;
		DAVectorSponge.TemperatureforClosenessTest = 9.0;

		// Set Data Reading parameters including Charge
		//  DAVectorSpongeSection.DistanceMatrixFile is basic
		//  DAVectorSpongeSection.LabelFile is update
		DAVectorSponge.SelectedInputLabel = -60;
		DAVectorSponge.Replicate = 1;
		DAVectorSponge.ComparisonInputFileType = 0;
		DAVectorSponge.ComparisonClusterFile = "E:\\Sali\\InCloud\\IUBox\\My Box Files\\Sponge\\input\\DarTBFULLSorted" +
                ".txt";
		DAVectorSponge.RestartSelectedInputLabel = -100000000;
		DAVectorSponge.RestartInputFileType = 1;
		DAVectorSponge.RestartClusterFile = DAVectorSponge.config.DistanceMatrixFile;
		DAVectorSponge.ComparisonSelectedInputLabel = 2;
		DAVectorSponge.SelectedInputLabel = 2; // Charge 2
		DAVectorSponge.InputFileType = 0; //  Original Harvard format
		DAVectorSponge.RestartTemperature = -1.0;
		DAVectorSponge.CompareSolution = -1;

		//  Sigma Annealing
		DAVectorSponge.SigmaVectorParameters_i_[0] = 0.000598;
		DAVectorSponge.SigmaVectorParameters_i_[1] = 2.35;
		//       DAVectorSponge.SigmaMethod = 2;
		DAVectorSponge.SigmaMethod = 3;
		DAVectorSponge.FinalTargetSigma0 = 0.00000598;
		DAVectorSponge.FinalTargetTemperature = 12.0;
		DAVectorSponge.FinalTargetTemperature = 30.0;

		//  Sponge Parameters
		DAVectorSponge.UseSponge = false;
		DAVectorSponge.SpongePWeight = 0.1;
		DAVectorSponge.CreateSpongeScaledSquaredWidth = 10.0;

		DAVectorSponge.SpongeFactor1 = 12.0;
		DAVectorSponge.SpongeFactor1 = 45.0;
		DAVectorSponge.SpongeFactor2 = 2.0;
		DAVectorSponge.SpongeTemperature1 = 9.0;
		DAVectorSponge.SpongeTemperature1 = 30.0;
		DAVectorSponge.SpongeTemperature2 = 2.5;
		DAVectorSponge.SpongeTemperature2 = 2.0;

		// No Sponge
		DAVectorSponge.CreateSpongeScaledSquaredWidth = -1.0; // Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches this
		DAVectorSponge.SpongeTemperature1 = -1.0; // Minimum Temperature where Sponge Introduced
		DAVectorSponge.SpongeTemperature2 = -1.0;

		// Run F7
		DAVectorSponge.SigmaVectorParameters_i_[0] = 0.000598;
		DAVectorSponge.SigmaMethod = 3;
		DAVectorSponge.FinalTargetSigma0 = 0.00000598;
		DAVectorSponge.FinalTargetTemperature = 12.0;

		DAVectorSponge.ScaledSquaredDistanceatClosenessTest = 0.5;
		DAVectorSponge.TemperatureforClosenessTest = 4.0;
		DAVectorSponge.MaxNumberSplitClusters = 256;
		DAVectorSponge.MinimumScaledWidthsquaredtosplit = 1.5;
		DAVectorSponge.ClusterLimitforDistribution = 256;
		DAVectorSponge.Tminimum = 0.1;
		DAVectorSponge.RestartTemperature = -1.0;
		DAVectorSponge.Iterationatend = 400;

		DAVectorSponge.UseSponge = false;
		DAVectorSponge.SpongeFactor1 = 12.0;
		DAVectorSponge.SpongeFactor2 = 3.0;
		DAVectorSponge.SpongeTemperature1 = 9.0;
		DAVectorSponge.SpongeTemperature2 = 3.0;
		DAVectorSponge.SpongePWeight = 0.1;
		DAVectorSponge.CreateSpongeScaledSquaredWidth = 10.0;

		DAVectorSponge.CoolingTemperatureSwitch = 12.0;
		DAVectorSponge.InitialCoolingFactor1 = 0.9875000;
		DAVectorSponge.FineCoolingFactor1 = 0.9975000;
		DAVectorSponge.FineCoolingFactor2 = 0.9993750;
		DAVectorSponge.InitialCoolingFactor2 = 0.9937500;

		DAVectorSponge.Malpha_MaxChange = 0.005;
		DAVectorSponge.Malpha_MaxChange1 = 0.005;
		DAVectorSponge.MinimumCountforCluster_C_kwithSponge = 1.5;
		DAVectorSponge.MinimumCountforCluster_C_k = 0.5;
		DAVectorSponge.MinimumCountforCluster_Points = -1;

		//  F7 PLUS
		*//*
		DAVectorSponge.InitialCoolingFactor1 = Math.Sqrt(DAVectorSponge.InitialCoolingFactor1);
		DAVectorSponge.FineCoolingFactor1 = Math.Sqrt(DAVectorSponge.FineCoolingFactor1);
		DAVectorSponge.FineCoolingFactor2 = Math.Sqrt(DAVectorSponge.FineCoolingFactor2);
		DAVectorSponge.InitialCoolingFactor2 = Math.Sqrt(DAVectorSponge.InitialCoolingFactor2);
		*/

	} // End FreshFullAnalysis()

    // TODO - Test Method to Setup SetupCharge5_From_TempestF7_1_1_8_2June2013
    // The only change is SelectedInputLabel is set to 5 instead of 2
    private static void SetupCharge5_From_TempestF7_1_1_8_2June2013_w_noduplicate_assignments(){
        DAVectorSponge.maxNcentperNode = 200000;
        DAVectorSponge.maxNcentTOTAL = 200000;
        DAVectorSponge.targetNcentperPoint = 8;
        DAVectorSponge.maxNcentperPoint = 17;
        DAVectorSponge.targetMinimumNcentperPoint = 2;
        DAVectorSponge.maxNcentCreated = 5 * DAVectorSponge.maxNcentperNode;
        DAVectorSponge.Waititerations = 1;
        DAVectorSponge.Waititerations_Converge = 4;
        DAVectorSponge.PrintInterval = 50;
        DAVectorSponge.ClusterPrintNumber = 5;
        DAVectorUtility.DebugPrintOption = 1;
        DAVectorSponge.MagicTemperatures[0] = -1.0;
        DAVectorSponge.MagicTemperatures[1] = -1.0;
        DAVectorSponge.MagicTemperatures[2] = -1.0;
        DAVectorSponge.MagicTemperatures[3] = -1.0;
        DAVectorSponge.MagicTemperatures[4] = -1.0;
        DAVectorSponge.ExpArgumentCut1 = 20.0;
        DAVectorSponge.ExpArgumentCut2 = 40.0;
        DAVectorSponge.Replicate = 1;
        DAVectorSponge.ComparisonInputFileType = 0;
        DAVectorSponge.ComparisonClusterFile = "/home/saliya/sali/csharpToJava/sponge/input/DarTBFULLSorted.txt";
        DAVectorSponge.RestartInputFileType = 1;
        DAVectorSponge.RestartClusterFile = DAVectorSponge.config.DistanceMatrixFile;
        DAVectorSponge.ComparisonSelectedInputLabel = 2;
        DAVectorSponge.SelectedInputLabel = 5; // Charge 5
        DAVectorSponge.InputFileType = 0;
        DAVectorSponge.CompareSolution = -1;
        DAVectorSponge.SigmaVectorParameters_i_[1] = 2.35;
        DAVectorSponge.SigmaVectorParameters_i_[0] = 0.000598;
        DAVectorSponge.SigmaMethod = 3;
        DAVectorSponge.FinalTargetSigma0 = 5.98E-06;
        DAVectorSponge.FinalTargetTemperature = 12.0;
        DAVectorSponge.ScaledSquaredDistanceatClosenessTest = 0.5;
        DAVectorSponge.TemperatureforClosenessTest = 4.0;
        DAVectorSponge.MaxNumberSplitClusters = 256;
        DAVectorSponge.MinimumScaledWidthsquaredtosplit = 1.5;
        DAVectorSponge.ClusterLimitforDistribution = 256;
        DAVectorSponge.Tminimum = 0.1;
        DAVectorSponge.RestartTemperature = -1.0;
        DAVectorSponge.Iterationatend = 400;
        DAVectorSponge.UseSponge = false;
        DAVectorSponge.SpongeFactor1 = 12.0;
        DAVectorSponge.SpongeFactor2 = 3.0;
        DAVectorSponge.SpongeTemperature1 = 9.0;
        DAVectorSponge.SpongeTemperature2 = 3.0;
        DAVectorSponge.SpongePWeight = 0.1;
        DAVectorSponge.CreateSpongeScaledSquaredWidth = 10.0;
        DAVectorSponge.CoolingTemperatureSwitch = 12.0;
        DAVectorSponge.InitialCoolingFactor1 = 79.0 / 80.0;
        DAVectorSponge.FineCoolingFactor1 = 399.0 / 400.0;
        DAVectorSponge.FineCoolingFactor2 = 0.999375;
        DAVectorSponge.InitialCoolingFactor2 = 159.0 / 160.0;
        DAVectorSponge.Malpha_MaxChange = 0.005;
        DAVectorSponge.Malpha_MaxChange1 = 0.005;
        DAVectorSponge.MinimumCountforCluster_C_kwithSponge = 1.5;
        DAVectorSponge.MinimumCountforCluster_C_k = 0.5;
        DAVectorSponge.MinimumCountforCluster_Points = -1;
    }

    /*// TODO - Test Method to Setup SetupCharge5_From_TempestF7_1_1_8_2June2013
    // The only change is SelectedInputLabel is set to 5 instead of 2
    private static void SetupCharge5_From_TempestF7_1_1_8_2June2013()
    {
        DAVectorSponge.maxNcentperNode = 200000;
        DAVectorSponge.maxNcentTOTAL = 200000;
        DAVectorSponge.targetNcentperPoint = 8;
        DAVectorSponge.maxNcentperPoint = 17;
        DAVectorSponge.targetMinimumNcentperPoint = 2;

        DAVectorSponge.maxNcentCreated = 5 * DAVectorSponge.maxNcentperNode;
        DAVectorSponge.Iterationatend = 2000;
        DAVectorSponge.Malpha_MaxChange = 0.005;
        DAVectorSponge.Tminimum = 0.025;
        DAVectorSponge.MaxNumberSplitClusters = 256;
        DAVectorSponge.ClusterLimitforDistribution = 256;
        DAVectorSponge.Waititerations = 1;

        DAVectorSponge.Waititerations_Converge = 4;

        DAVectorSponge.PrintInterval = 50;
        DAVectorSponge.ClusterPrintNumber = 5;
        DAVectorUtility.DebugPrintOption = 1;
        DAVectorSponge.MagicTemperatures[0] = -1.0;
        DAVectorSponge.MagicTemperatures[1] = -1.0;
        DAVectorSponge.MagicTemperatures[2] = -1.0;
        DAVectorSponge.MagicTemperatures[3] = -1.0;
        DAVectorSponge.MagicTemperatures[4] = -1.0;
        DAVectorSponge.InitialCoolingFactor1 = 79.0 / 80.0;
        DAVectorSponge.FineCoolingFactor1 = 399.0 / 400.0;
        DAVectorSponge.InitialCoolingFactor2 = 159.0 / 160.0;
        DAVectorSponge.FineCoolingFactor2 = 0.999375;
        DAVectorSponge.CoolingTemperatureSwitch = 30.0;
        DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
        DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
        DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
        DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);
        DAVectorSponge.InitialCoolingFactor1 = Math.sqrt(DAVectorSponge.InitialCoolingFactor1);
        DAVectorSponge.FineCoolingFactor1 = Math.sqrt(DAVectorSponge.FineCoolingFactor1);
        DAVectorSponge.FineCoolingFactor2 = Math.sqrt(DAVectorSponge.FineCoolingFactor2);
        DAVectorSponge.InitialCoolingFactor2 = Math.sqrt(DAVectorSponge.InitialCoolingFactor2);
        DAVectorSponge.ExpArgumentCut1 = 20.0;
        DAVectorSponge.ExpArgumentCut2 = 40.0;
        DAVectorSponge.MinimumScaledWidthsquaredtosplit = 0.5;
        DAVectorSponge.ScaledSquaredDistanceatClosenessTest = 0.35;
        DAVectorSponge.MinimumCountforCluster_C_kwithSponge = 1.5;
        DAVectorSponge.MinimumCountforCluster_C_k = 0.5;
        DAVectorSponge.MinimumCountforCluster_Points = -1;
        DAVectorSponge.TemperatureforClosenessTest = 4.0;
        DAVectorSponge.TemperatureforClosenessTest = 9.0;
        DAVectorSponge.SelectedInputLabel = -60;
        DAVectorSponge.Replicate = 1;
        DAVectorSponge.ComparisonInputFileType = 0;
        DAVectorSponge.ComparisonClusterFile = "/home/saliya/sali/csharpToJava/sponge/input/DarTBFULLSorted.txt";
        DAVectorSponge.RestartSelectedInputLabel = -100000000;
        DAVectorSponge.RestartInputFileType = 1;
        DAVectorSponge.RestartClusterFile = DAVectorSponge.config.DistanceMatrixFile;
        DAVectorSponge.ComparisonSelectedInputLabel = 2;
        DAVectorSponge.SelectedInputLabel = 5; // Charge 5
        DAVectorSponge.InputFileType = 0;
        DAVectorSponge.RestartTemperature = -1.0;
        DAVectorSponge.CompareSolution = -1;
        DAVectorSponge.SigmaVectorParameters_i_[0] = 0.000598;
        DAVectorSponge.SigmaVectorParameters_i_[1] = 2.35;
        DAVectorSponge.SigmaMethod = 3;
        DAVectorSponge.FinalTargetSigma0 = 5.98E-06;
        DAVectorSponge.FinalTargetTemperature = 12.0;
        DAVectorSponge.FinalTargetTemperature = 30.0;
        DAVectorSponge.UseSponge = false;
        DAVectorSponge.SpongePWeight = 0.1;
        DAVectorSponge.CreateSpongeScaledSquaredWidth = 10.0;
        DAVectorSponge.SpongeFactor1 = 12.0;
        DAVectorSponge.SpongeFactor1 = 45.0;
        DAVectorSponge.SpongeFactor2 = 2.0;
        DAVectorSponge.SpongeTemperature1 = 9.0;
        DAVectorSponge.SpongeTemperature1 = 30.0;
        DAVectorSponge.SpongeTemperature2 = 2.5;
        DAVectorSponge.SpongeTemperature2 = 2.0;
        DAVectorSponge.CreateSpongeScaledSquaredWidth = -1.0;
        DAVectorSponge.SpongeTemperature1 = -1.0;
        DAVectorSponge.SpongeTemperature2 = -1.0;
        DAVectorSponge.SigmaVectorParameters_i_[0] = 0.000598;
        DAVectorSponge.SigmaMethod = 3;
        DAVectorSponge.FinalTargetSigma0 = 5.98E-06;
        DAVectorSponge.FinalTargetTemperature = 12.0;
        DAVectorSponge.ScaledSquaredDistanceatClosenessTest = 0.5;
        DAVectorSponge.TemperatureforClosenessTest = 4.0;
        DAVectorSponge.MaxNumberSplitClusters = 256;
        DAVectorSponge.MinimumScaledWidthsquaredtosplit = 1.5;
        DAVectorSponge.ClusterLimitforDistribution = 256;
        DAVectorSponge.Tminimum = 0.1;
        DAVectorSponge.RestartTemperature = -1.0;
        DAVectorSponge.Iterationatend = 400;
        DAVectorSponge.UseSponge = false;
        DAVectorSponge.SpongeFactor1 = 12.0;
        DAVectorSponge.SpongeFactor2 = 3.0;
        DAVectorSponge.SpongeTemperature1 = 9.0;
        DAVectorSponge.SpongeTemperature2 = 3.0;
        DAVectorSponge.SpongePWeight = 0.1;
        DAVectorSponge.CreateSpongeScaledSquaredWidth = 10.0;
        DAVectorSponge.CoolingTemperatureSwitch = 12.0;
        DAVectorSponge.InitialCoolingFactor1 = 79.0 / 80.0;
        DAVectorSponge.FineCoolingFactor1 = 399.0 / 400.0;
        DAVectorSponge.FineCoolingFactor2 = 0.999375;
        DAVectorSponge.InitialCoolingFactor2 = 159.0 / 160.0;
        DAVectorSponge.Malpha_MaxChange = 0.005;
        DAVectorSponge.Malpha_MaxChange1 = 0.005;
        DAVectorSponge.MinimumCountforCluster_C_kwithSponge = 1.5;
        DAVectorSponge.MinimumCountforCluster_C_k = 0.5;
        DAVectorSponge.MinimumCountforCluster_Points = -1;
    }*/



    public static void LCMSCalculateClusterStatus() throws MPIException {
		double cut = DAVectorSponge.NearbySpongePointLimit;
		if (cut < 0.0)
		{
			cut = DAVectorSponge.SpongeFactor;
		}
		if (cut < 0.0)
		{
			cut = 1.0;
		}
		ClusteringSolution.SetGlobalClusterNumbers();
		VectorAnnealIterate.OutputClusteringResults(""); // Set Point -- Cluster links

		DAVectorSponge.ClusterStatus = new ClusterQuality(ClusteringSolution.TotalClusterSummary.NumberofCenters, ClusteringSolution.TotalClusterSummary.SpongeCluster, DAVectorSponge.NumberNearbyClusters, cut);

		DAVectorSponge.ClusterStatus.SetClusterStatistics();
		DAVectorSponge.ClusterStatus.SetPointStatistics();
		DAVectorSponge.ClusterStatus.SetNearbyClusters();
		DAVectorSponge.ClusterStatus.OutputStatus();

		//  Set up Cluster Comparisons
		//  Convert Sponge to singleton clusters
		//  Note FullClusterNumber will end up as 1 more than total of clusters if there is a sponge
		if (DAVectorSponge.CompareSolution <= 0)
		{
			return;
		}
		int FullClusterNumber = ClusteringSolution.TotalClusterSummary.NumberofCenters;
		int SpongeCluster = ClusteringSolution.TotalClusterSummary.SpongeCluster;
		for (int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
		{
			int NewClusterIndex = DAVectorSponge.ClusterAssignments[GlobalPointIndex];
			if (NewClusterIndex == SpongeCluster)
			{
				NewClusterIndex = FullClusterNumber;
				++FullClusterNumber;
			}
			DAVectorSponge.OurClusters.PointstoClusterIDs[GlobalPointIndex] = NewClusterIndex;
		}
		DAVectorSponge.OurClusters.setup();
		LCMSAnalyze.ClusterComparison();

	} // End LCMSCalculateClusterStatus()

	public static void ClusterComparison()
	{
		if (DAVectorSponge.CompareSolution <= 0)
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

		if (DAVectorSponge.GoldenPeaks.MaxIndependentClusterIndex > 0)
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
					ThreeMethods[0] = DAVectorSponge.OurClusters;
					ThreeMethods[1] = DAVectorSponge.MedeaClusters;
					ThreeMethods[2] = DAVectorSponge.MclustClusters;
					ThreeMethods[3] = DAVectorSponge.GoldenPeaks;
					GoldenExaminationContainer.GoldenComparisons(DAVectorSponge.GoldenPeaks, ThreeMethods);
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
			DAVectorSponge.OurClusters.HistogramPeaks(ClusterCountcut);
			DAVectorSponge.MedeaClusters.HistogramPeaks(ClusterCountcut);
			DAVectorSponge.MclustClusters.HistogramPeaks(ClusterCountcut);
			DAVectorSponge.GoldenPeaks.HistogramPeaks(ClusterCountcut);
			ClusterCountcut += 5;
		}

		//  Comparison of Basic Clusters
		DAVectorUtility.SALSAPrint(0, "\n**************** Statistics of Clustering in each method and selection to Golden Clusters\nWith means versus occupation count and 1D/2D Point-Center Histograms");
		DAVectorSponge.OurClusters.Statistics();
		DAVectorSponge.MedeaClusters.Statistics();
		DAVectorSponge.MclustClusters.Statistics();
		if (DAVectorSponge.GoldenPeaks.MaxIndependentClusterIndex > 0)
		{
			DAVectorSponge.GoldenPeaks.Statistics();
		}

		if (DAVectorSponge.GoldenPeaks.MaxIndependentClusterIndex > 0)
		{
			DAVectorUtility.SALSAPrint(0, "\n Comparisons with Golden Clusters");
			DAVectorSponge.OurClusters.Difference(DAVectorSponge.GoldenPeaks);
			DAVectorSponge.MedeaClusters.Difference(DAVectorSponge.GoldenPeaks);
			DAVectorSponge.MclustClusters.Difference(DAVectorSponge.GoldenPeaks);
			DAVectorSponge.GoldenPeaks.Difference(DAVectorSponge.GoldenPeaks);
		}

		DAVectorUtility.SALSAPrint(0, "\n Comparisons with DAVector Clusters");
		DAVectorSponge.MedeaClusters.Difference(DAVectorSponge.OurClusters);
		DAVectorSponge.MclustClusters.Difference(DAVectorSponge.OurClusters);
		DAVectorSponge.OurClusters.Difference(DAVectorSponge.MedeaClusters);
		DAVectorSponge.OurClusters.Difference(DAVectorSponge.MclustClusters);

	}

} // End LCMSAnalyze