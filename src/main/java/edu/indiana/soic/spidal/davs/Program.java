package edu.indiana.soic.spidal.davs;

// MPI Import
//import MPI.*;
import com.google.common.base.Optional;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.sections.DAVectorSpongeSection;
import edu.indiana.soic.spidal.davs.timing.*;
import edu.indiana.soic.spidal.general.Box;
import mpi.MPI;
import mpi.MPIException;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Date;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

public class Program
{
    private static Options programOptions = new Options();
    static {
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_C),Constants.CMD_OPTION_LONG_C, true,
                Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_N),Constants.CMD_OPTION_LONG_N,true,
                Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_T),Constants.CMD_OPTION_LONG_T,true,
                Constants.CMD_OPTION_DESCRIPTION_T);
    }

	//  Specify Clustering Algorithm
	public static boolean DoKmeans = false; //  If true use simple Kmeans
	public static boolean DoLCMS = false; // If true solve 2D LCMS problem
	public static boolean DoBasicDA = false; // If true do basic Vector DA
	public static boolean FollowUpKmeans = false; // If true follow up BasicDA with Kmeans

	//  Specify Nature of Vector Data
	public static int ParameterVectorDimension = 2; // Vector Dimension of Points

	//  Set Fixed errors in each dimension
	public static int SigmaMethod = 0; //   if 0, sigma of clusters unspecified. Otherwise specifies calculation method
	public static double[] SigmaVectorParameters_i_; // Control calculation of Cluster Sigma with dimension ParameterVectorDimension
	public static double FinalTargetTemperature = 3.0; // Used in SigmaMethod =3
	public static double FinalTargetSigma0 = 0.0; // Used in SigmaMethod =3
	public static double InitialSigma0 = 0.0; // Initial value

	//  Set Dataset properties
	public static int NumberDataPoints = -1; // if -1 then find number of points from input
	public static int SelectedInputLabel = 6; // Charge to select; if negative veto this label
	public static int RestartSelectedInputLabel = 6; // Charge to select in Restart; if negative veto this label
	public static int InputFileType = 0; // If 0 raw data; =1 Output file
	public static int RestartInputFileType = 1; // If 0 raw data; =1 Output file for Restart
	public static boolean Refinement = true; // If False do no refinement
	public static int StartPointPositiononInputLine = 1; // Starting position of Point Data on Input File
	public static int ClusterIndexonInputLine = -1; // Position of Cluster Index on Input File for K means
	public static int FirstClusterValue = 1; // Index for first cluster in Kmeans initialization
	public static double[][]  PointPosition; // Position of Point x is PointPosition[x, i] where i runs over vector space directions 0 ... ParameterVectorDimension-1
	public static int[] PointOriginalIndex; // Point Index on File
	public static int[] PointOriginalExperimentNumber; // Point experiment Number
	public static int[] PointLabel; // Label on File
	public static int Replicate = 1; // Replicate points for scaling tests
	public static double[][] FullPoint3DPosition; // 3D Position of Points -- all stored
	public static int RW3DData = -1; // If positive read and write Plotviz files of this dimension

	//  Control Execution Model in general
	public static String ControlFileName = ""; // Control File Name
	public static int ProcessingOption = 0; // Processing Option
	public static int cachelinesize = 0; // Specify Cache line minimum size
	public static boolean CalculateCorrelationMatrix = false; // If True calculate correlations
	public static boolean CalculateIndividualWidths = false; // If True calculate individual widths
	public static boolean CalculateEigenvaluesfromMatrix = false; // If false use iteration

	//  Describe mode for parallelism over clusters (as opposed to normal parallelism over points)
	public static int ClusterLimitforDistribution = -1; // Go into distributed mode if Number of Clusters per node greater than this
	public static int ActualClusterNumberforDistribution = -1; // Number of Clusters per node when go into distributed mode
	public static double TemperatureLimitforDistribution = -1.0; // Go into distributed mode if Temperature less than this
	public static double ActualTemperatureforDistribution = -1.0; // Temperature when go into distributed mode
	public static double TimeatDistribution = 0.0; // Time when distributed execution started

	public static int MaxIntegerComponents = 2; // Maximum Number of Integer Components in Transport (for MPI arrays)
	public static int MaxDoubleComponents = 3; // Maximum Number of Double Components in Transport
	public static int MaxMPITransportBuffer = 500; // Maximum Number of Clusters in MPI move routines
	public static int ActualMaxMPITransportBuffer = 0; // Needed Maximum Number of Clusters in MPI move routines
	public static int MaxNumberAccumulationsperNode = 30000; // Maximum Number of Clusters for a single node in Distributed Mode
	public static int ActualMaxNumberAccumulationsperNode = 0; // Needed Maximum Number of Clusters for a single node in Distributed Mode
	public static int MaxTransportedClusterStorage = 500; // Maximum Number of Remote Clusters stored on a node
	public static int ActualMaxTransportedClusterStorage = 0; // Needed Maximum Number of Remote Clusters stored on a node

	//  Specify a restart from existing file
	public static double RestartTemperature = -1.0; // If positive, Restart based on previous run(s)
	public static String RestartClusterFile = ""; // First file for restart

	//  Specify sizes of arrays dependent on cluster sizes
	public static int maxNcentperNode = 0; // maximum number of cluster centers needed in any one node
	public static int maxNcentTOTAL = 0; // maximum number of cluster centers needed in full arrays
	public static int targetNcentperPoint = 20; // Target Maximum number of cluster centers for each point (includes Sponge)
	public static int targetMinimumNcentperPoint = 1; // Target Minimum number of cluster centers for each point (includes Sponge)
	public static int maxNcentperPoint = 25; // Actual Maximum number of cluster centers for each point (includes Sponge)
	public static int maxNcentCreated = 200; // Maximum Number of clusters that can be created

	public static int InitialNcent = 1; // Initial value of Ncent

	//  Specify Annealing Schedule
	public static double InitialCoolingFactor = 0.9; // InitialCooling Factor in Annealing
	public static double FineCoolingFactor = 0.99; // Refined Cooling Factor in Annealing
	public static double InitialCoolingFactor1 = 0.9; // InitialCooling Factor in Annealing
	public static double FineCoolingFactor1 = 0.99; // Refined Cooling Factor in Annealing
	public static double InitialCoolingFactor2 = 0.9; // InitialCooling Factor in Annealing
	public static double FineCoolingFactor2 = 0.99; // Refined Cooling Factor in Annealing
	public static double CoolingTemperatureSwitch = 12.0; // Temperature when one should switch Cooling Factors
	public static int NumberTemperatureSteps = 0; // Count Temperature Steps

	public static double Tminimum = -1000.0; // If Positive this is minmum Temperature, If negative Divide maximum temperature by this to get target minimum
	public static double ActualStartTemperature; // Actual Starting Temperature;
	public static double ActualEndTemperatureafterconverging; // Final Temperature after Converging
	public static double ActualEndTemperature; // Temperature where we decided to stop
	public static double TimeatSplittingStop; // Time when splitting switched off
	public static double TargetEndTemperature; // Temperature we aimed at

	//  Parameters to specify terms that are too small
	public static double ExpArgumentCut1 = 20.0; // For numerical stability in EM Iteration DAVectorEMIterate
	public static double ExpArgumentCut2 = 40.0; // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut2 in Point Collection (used in public void SetClustersforaPoint)
	public static double ExpArgumentCut3 = 50.0; // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut3 in distributed broadcast (only calculation 0 component) in ClusterHostlinkage in DistributedClusteringSolution
	public static double ExpArgumentCut4 = 80.0; // For inclusion in Triangle Inequality

	//  Convergence and Iteration parameters
	public static int Waititerations = 1; // Wait this number of Temperature iterations before splitting
	public static int Waititerations_Converge = 4; // Wait this number of Temperature iterations before splitting after Cluster Removal

	public static int Iterationatend = 2000; // Finish up EM Loop with at most this number of iterations
	public static int ConvergenceLoopLimit = 20; // Limit on EM Convergence for each step

	public static double FreezingLimit = 0.002; // In finish stop when all freezing measures are < FreezingLimit
	public static double Malpha_MaxChange = 0.005; // Change in average M per point to define convergence in a step (replaces epsilon limit used in Pairwise Case)
	public static double Malpha_MaxChange1 = 0.0005; // Change in average M per point to define convergence in a step for final steps
	public static double YChangeSquared = 0.000001; // Change compared to width of sum of y squared changes

	public static int PowerIterationLimit = 200; //   Limit for Power Iterations
	public static double eigenvaluechange = 0.001; // Limit 1 on Eigenvalue Changes
	public static double eigenvectorchange = 0.001; // Limit 2 on Eigenvalue Changes

	//  Specify Nature of Splitting
	public static int MaxNumberSplitClusters = 3; // System will split upto this number per node (in distributed mode) simultaneously
	public static double ToosmalltoSplit = 4.0; // Size of Cluster measured by C_k_ that should not be split
	public static double MinimumScaledWidthsquaredtosplit = 6.0; // Do not split Cluster whose scaled width is less than this

	//  Specify annealing clean up and way it deals with close clusters
	public static double MinimumCountforCluster_C_k = 0.5; // Remove clusters with fewer points than this (based on C_k_) (This only used for converged clusters)
	public static double MinimumCountforCluster_C_kwithSponge = 1.5;
	public static int MinimumCountforCluster_Points = 2; // Remove clusters with fewer points than this (based on assigned points)
	public static double CountforCluster_C_ktobezero = 0.001; // This limit used to define a zero size cluster while system running. Positions Y are not evolved and Widths set to zero in this case

	public static double TemperatureforClosenessTest = 6.0; // Validity Check removes Close Clusters below this
	public static double ScaledSquaredDistanceatClosenessTest = 1.0; // Coalesce Clusters separated by this
	public static double[] MagicTemperatures = {4.0, 3.5, 2.5, 2.0, 1.75}; // Temperatures for clean up
	public static int magicindex = 0; // Index magic temperatures

	//  Specify Sponge Features
	public static boolean UseSponge = false; // If True use a single Sponge Factor
	public static double SpongeFactor1 = 3.0; //  Initial Sponge Factor
	public static double SpongeFactor2 = 3.0; // Final Sponge Factor
	public static double SpongeFactor = SpongeFactor1; // Sponge Factor has an exp(-SpongeFactor^2/2T)
	public static int SpongePoption = 1; // Set method to calculate sponge normalization P(SpongeCluster) = 0 default; =1 Set to average of others
	public static double SpongePWeight = 0.1;
	public static double CreateSpongeScaledSquaredWidth = -1.0; // Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches this
	public static double SpongeTemperature1 = -1.0; // Minimum Temperature where Sponge Introduced
	public static double ActualSpongeTemperature = -1.0; // Temperature where Sponge Introduced
	public static double ActualWidth = -1.0; // Width when Sponge introduced
	public static double TimeatSponge = 0.0; // Time when Sponge added
	public static double SpongeTemperature2 = -1.0; // Temperature where Final Sponge Factor introduced
	public static double AddSpongeScaledWidthSquared = -1.0; // If Positive add sponge cluster when averaged scaled width less than this

	//  Specify Deterministic Annealing Style
	public static boolean ContinuousClustering = true; // If true use the Ken Rose Continuous Clustering
	public static boolean ConvergeIntermediateClusters = false; // NOT used

	//  Specify Diagnostic Print Out
	public static int ClusterPrintNumber = 5; // Print the last number of these clusters
	public static int PrintInterval = 3; // Output at this interval in iterations
	public static boolean RemovalDiagnosticPrint = false; // If true output histograms on clusters per point
	public static int MaxNumberClustersToPrint = 200; // Maximum Number of Clusters to Print
	public static boolean Printeigenvectors = false; // If true print eigenvectors in Shouldwesplit()

	//  Set Quantities that record properties of clustering
	public static double MdiffSum = 0.0; // Accumulate Mdiff values for average
	public static int NumberMdiffSums = 0; // Number of MdiffAvg's in Sum
	public static int IterationsperStepSum = 0; // Accumulate number of iterations per step
	public static int NumberIterationSteps = 0; // Number of steps
	public static int NumberMajorSynchs1 = 0; // Number of calls to Major Synchromization in Global mode
	public static int NumberMajorSynchs2 = 0; // Number of calls to Major Synchromization in Distributed mode
	public static int CountClusters2 = 0; //  Count clusters while distributed
	public static int NumberMinorSynchs = 0; // Number of calls to Major Synchromization in Distributed mode
	public static int NumberPipelineSteps = 0; // Number of Pipeline steps
	public static int NumberofPipelineClusters = 0; // Number of Cluster-Steps in Pipelines
	public static int NumberPipelineGroups = 0; // Number of Pipeline group calls
	public static String FinalReason = ""; // Final reason to stop


	public static double SumUsefulCalcs = 0.0; // Number of Useful distance calcs
	public static double SumUselessCalcs = 0.0; // Number of Useless distance calcs
	public static double SumIgnoredCalcs = 0.0; // Number of Ignored distance calcs
	public static double SumEigenSPCalcs = 0.0; // Sum eigenvector scalar product computations

	public static int TotalClustersDeleted_CSmall = 0; // Number of Clusters deleted as C_k_ too small
	public static int TotalClustersDeleted_OccCount = 0; // Number of Clusters deleted as Occupation Count too small
	public static int TotalClustersDeleted_Close = 0; // Number of Clusters deleted as too close

	public static double NumberMsuccesses = 0.0; // Number of steps where M succeeds
	public static double NumberYsuccesses = 0.0; // Number of steps where Y succeeds
	public static double NumberMfailures = 0.0; // Number of steps where M fails
	public static double NumberYfailures = 0.0; // Number of steps where Y fails
	public static double AccumulateMvalues = 0.0; // Number of iterations where M succeeeds
	public static double AccumulateYvalues = 0.0; // Number of iterations where Y succeeds

	//  Kmeans Paramters
	public static double KmeansCenterChangeStop = 0.001; // Stop when relative to radius center change reaches this limit
	public static int KmeansIterationLimit = 1000; // Stop at this iteration limit

	//  Triangle Inequality for Kmeans Parameters
	public static int maxNcentTOTALforParallelism_Kmeans = 50; // Use Center parallelism when this limit hit
	public static int OldCenterOption_Kmeans = 0; // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations
	public static boolean DoBackwardFacingTests_Kmeans = true; // If True do backward facing tests
	public static int UseTriangleInequality_Kmeans = 0; // if 0 do NOT use triangle inequality; > 0 use it in a way specified by integer
	public static double TriangleInequality_Delta1_old_Kmeans = 0.1; // Test for center change and old lower bounds (normalized by radius)
	public static double TriangleInequality_Delta1_current_Kmeans = 0.1; // Test for Center change and current lower bounds (normalized by radius)
	public static int MaxClusterLBsperPoint_Kmeans = 50; // Maximum number of Lower Bound values
	public static int MaxCentersperCenter_Kmeans = 49; // Maximum number of Centers in Center Difference Array

	//  Triangle Inequality for DA Parameters
	public static int maxNcentTOTALforParallelism_DA = 50; // Use Center parallelism when this limit hit
	public static int OldCenterOption_DA = 0; // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations
	public static int UseTriangleInequality_DA = 0; // if 0 do NOT use triangle inequality; > 0 use it in a way specified by integer
	public static double TriangleInequality_Delta1_old_DA = 0.1; // Test for center change and old lower bounds (normalized by radius)
	public static double TriangleInequality_Delta1_current_DA = 0.1; // Test for Center change and current lower bounds (normalized by radius)
	public static int MaxClusterLBsperPoint_DA = 50; // Maximum number of Lower Bound values
	public static int MaxCentersperCenter_DA = 49; // Maximum number of Centers in Center Difference Array
	public static double[] ClustersperCenterDiagnostics = new double[8]; // Diagnostics

	//  Very Specific LC-MS Analysis ("print out") of Results
	public static ArbitraryClustering GoldenPeaks; // Certified Peaks
	public static ArbitraryClustering MedeaClusters; // Medusa Clusters
	public static ArbitraryClustering MclustClusters; // Mclust Clusters
	public static ArbitraryClustering OurClusters; // Our Clusters

	public static int CompareSolution = -1; // Control solution comparison
	public static String ComparisonClusterFile = ""; // File for Comparison Data
	public static int ComparisonSelectedInputLabel = 6; // Charge to select in Comparison; if negative veto this label
	public static int ComparisonInputFileType = 0; // If 0 raw data; =1 Output file

	public static int ClusterCountOutput = 0; // Control Label Output = -1 not at all, = 0 at end only, = 1 at each count
	public static int[] ClusterAssignments; // This gives for all points their cluster assignments (set in OutputClusterLabels)
	public static int[] ExperimentNumberAssigments; // This gives for all points their cluster assignments (set in OutputClusterLabels)
	public static ClusterQuality ClusterStatus;
	public static int ClusterNumberOutput = -1;
	public static int NumberNearbyClusters = 5; // specify number of nearby clusters to output
	public static double NearbySpongePointLimit = -1.0; // Multiplier for sigma used in counting sponge points near a cluster; if negative use spongefactor

	//Cluster Quality Parameters
	public static int HistogramOccupationMax;
	public static int HistogramDistancesMax;

	// Config Settings
	public static DAVectorSpongeSection config;
    public static String ClusterFile;
    public static String DistanceMatrixFile;
    public static String LabelFile;
    public static String TimingFile;
    public static String SummaryFile;

    /**
     * Parse command line arguments
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(String [] args, Options opts){

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        } catch (ParseException e) {
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }

	//  User routine to specify Sigma as a function of Cluster Position given in Centre
	public static void CalculateSigma(double[] Centre, Box<double[]> Sigma)
	{

		if (Program.SigmaMethod == 0)
		{
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				Sigma.content[VectorIndex] = 1.0;
			}
			return;
		}
		if (Program.SigmaMethod == 1)
		{
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				Sigma.content[VectorIndex] = Program.SigmaVectorParameters_i_[VectorIndex] * Program.SigmaVectorParameters_i_[VectorIndex];
			}
			return;
		}
		if ((Program.SigmaMethod == 2) || (Program.SigmaMethod == 3))
		{
            System.arraycopy(Program.SigmaVectorParameters_i_, 0, Sigma.content, 0,
                    Program.ParameterVectorDimension);
			Sigma.content[0] = Sigma.content[0] * Centre[0];
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				Sigma.content[VectorIndex] = Sigma.content[VectorIndex] * Sigma.content[VectorIndex];
			}
			return;
		}

    } // End CalculateSigmaSquared

	public static boolean ChangeClusterSigmas(double TemperatureRatioChange, ClusteringSolution Solution) throws MPIException {
		if (Program.SigmaMethod != 3)
		{
			return false;
		}

		boolean asymptotic = (Program.SigmaVectorParameters_i_[0] <= Program.FinalTargetSigma0);
        asymptotic  = DAVectorUtility.SynchronizeMPIBoolean(asymptotic);

		if (asymptotic)
		{
			return false;
		}
		boolean justsetit = Solution.Temperature <= Program.FinalTargetTemperature;
        justsetit = DAVectorUtility.SynchronizeMPIBoolean(justsetit);
		if (justsetit)
		{
			Program.SigmaVectorParameters_i_[0] = Program.FinalTargetSigma0;
		}
		else
		{
			double LogTemperatureChange = Math.min(Math.log(Program.FinalTargetTemperature) - Math.log(Solution.Temperature), Math.log(TemperatureRatioChange));
			double ActualRatio = Math.exp(Math.log(TemperatureRatioChange) * (Math.log(Program.FinalTargetSigma0) - Math.log(
                    Program.SigmaVectorParameters_i_[0])) / LogTemperatureChange);
			Program.SigmaVectorParameters_i_[0] = Math.max(Program.FinalTargetSigma0, Program.SigmaVectorParameters_i_[0] * ActualRatio);
		}

		for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			Box<double[]> tempRef_Object = new Box<>(Solution.Sigma_k_i_[RealClusterIndex]);
			Program.CalculateSigma(Solution.Y_k_i_[RealClusterIndex], tempRef_Object);
			Solution.Sigma_k_i_[RealClusterIndex] = tempRef_Object.content;
		}
		return true;
	} // End ChangeClusterSigmas

	public static boolean ChangeSpongeFactor(double TemperatureRatioChange, ClusteringSolution Solution) throws MPIException {

		if (Solution.SpongeCluster < 0)
		{
			return false;
		}
		if ((Program.SpongeTemperature2 < 0.0) || (Program.SpongeFactor2 < 0.0))
		{
			return false;
		}
		boolean asymptotic = (Program.SpongeFactor <= Program.SpongeFactor2);
		asymptotic = DAVectorUtility.SynchronizeMPIBoolean(asymptotic);
		if (asymptotic)
		{
			return false;
		}

		boolean justsetit = Solution.Temperature <= Program.SpongeTemperature2;
		justsetit = DAVectorUtility.SynchronizeMPIBoolean(justsetit);
		if (justsetit)
		{
			Program.SpongeFactor = Program.SpongeFactor2;
		}
		else
		{
			double LogTemperatureChange = Math.min(Math.log(Program.SpongeTemperature2) - Math.log(Solution.Temperature), Math.log(TemperatureRatioChange));
			double ActualRatio = Math.exp(Math.log(TemperatureRatioChange) * (Math.log(Program.SpongeFactor2) - Math.log(
                    Program.SpongeFactor)) / LogTemperatureChange);
			Program.SpongeFactor = Program.SpongeFactor * ActualRatio;
			Program.SpongeFactor = Program.SpongeFactor2 + (Program.SpongeFactor1 - Program.SpongeFactor2) * (Solution.Temperature - Program.SpongeTemperature2) / (Program.SpongeTemperature1 - Program.SpongeTemperature2);
			Program.SpongeFactor = Math.max(Program.SpongeFactor2, Program.SpongeFactor);
		}
		return true;

	} // End ChangeSpongeFactor

    /**
     * Parallel Vector clustering based on Kmeans or Deterministic Annealing algorithm
     * @param args command line arguments to the program, which should include
     *             -c "path to config file" -t "number of threads" -n "number of nodes"
     *             The options may also be given as longer names
     *             --configFile, --threadCount, and --nodeCount respectively
     */
	public static void main(String[] args) throws MPIException {
		GeneralTiming.startTiming(GeneralTiming.TimingTask.TOTAL);
		SetupTiming.startTiming(SetupTiming.TimingTask.SETUP);
        Optional<CommandLine> parserResult = parseCommandLineArguments(args, programOptions);
        if (!parserResult.isPresent()){
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_T))){
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

		//  Read Metadata using this as source of other metadata
		ReadControlFile(cmd);

        // TODO - Unnecessary as reading this from config as an array.
		//  Initialize Clusters
		/*SigmaVectorParameters_i_ = new double[Program.ParameterVectorDimension];*/

		Program.DoKmeans = false;
		Program.DoLCMS = true;
		Program.DoBasicDA = false;
		Program.FollowUpKmeans = true;

		if (Program.DoLCMS)
		{
			Program.DoBasicDA = false;
		}
		if (Program.DoBasicDA)
		{
			Program.DoKmeans = false;
		}
		if (!Program.DoBasicDA)
		{
			Program.FollowUpKmeans = false;
		}

		Program.SetupLCMS();
		if (!Program.DoLCMS)
		{
			Program.CompareSolution = -1;
		}
		Program.SetupKmeans();
		Program.SetupDA();

		if (Program.UseTriangleInequality_DA > 0)
		{
			for (int DiagnosticLoop = 0; DiagnosticLoop < 6; DiagnosticLoop++)
			{
				Program.ClustersperCenterDiagnostics[DiagnosticLoop] = 0.0;
			}
		}
		if (!Program.CalculateEigenvaluesfromMatrix)
		{
			Program.CalculateCorrelationMatrix = false;
		}

		//  Don't Change anything after this
		if (!Program.Refinement)
		{
			Program.RestartTemperature = Tminimum;
			RestartInputFileType = 1;
			RestartClusterFile = config.DistanceMatrixFile;
			config.LabelFile = "";
		}

		// if (Program.RestartTemperature > 0.0)
		//   Program.UseSponge = true;

		Program.InitialCoolingFactor2 = Math.max(Program.InitialCoolingFactor2, Program.InitialCoolingFactor1);
		Program.FineCoolingFactor2 = Math.max(Program.FineCoolingFactor2, Program.FineCoolingFactor1);
		Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
		Program.FineCoolingFactor = Program.FineCoolingFactor1;

		// Set Initial Values of Annealed Parameters
		Program.InitialSigma0 = Program.SigmaVectorParameters_i_[0];
		Program.SpongeFactor = Program.SpongeFactor1;


		//  Set up MPI and threads parallelism
        try {
            DAVectorParallelism.SetupParallelism(args);
        } catch (MPIException e) {
            DAVectorUtility.printException(e);
            return; // End program on error
        }

        // Note - Logging
        /*initializeLogging();*/

        //  Restrict number of splits per node
		Program.MaxNumberSplitClusters = Math.min(Program.MaxNumberSplitClusters, DAVectorUtility.MPI_Size * 32);
		if (!Program.DoKmeans)
		{
			Program.maxNcentperPoint = Math.max(Program.maxNcentperPoint, Program.targetNcentperPoint + Program.MaxNumberSplitClusters + 3);
		}
		int SaveSplitNumber = Program.MaxNumberSplitClusters;

		//  Find number of points if requested
		if (Program.NumberDataPoints < 1)
		{
			if (Program.DoLCMS)
			{
				if (DAVectorUtility.MPI_Rank == 0)
				{
					DAVectorReadData.AnalyzeDataFromFile(config.DistanceMatrixFile);
				}
                // Note - MPI Call - Broadcast - int
                if (DAVectorUtility.MPI_Size > 1){
//                    DAVectorUtility.MPI_communicator.<Integer>Broadcast(PointCount_Global, 0);
                    DAVectorUtility.PointCount_Global = DAVectorUtility.mpiOps.broadcast(DAVectorUtility.PointCount_Global, 0);
                }
			}
			else
			{
				DAVectorUtility.printAndThrowRuntimeException("Must Set Number of Input Data Points");
			}
		}
		else
		{
			DAVectorUtility.PointCount_Global = Program.NumberDataPoints;
		}

		DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global + " Selection " + Program.SelectedInputLabel + " Input File Type " + Program.InputFileType + " Comparison Selection " + Program.ComparisonSelectedInputLabel + " Comparison Input File Type " + Program.ComparisonInputFileType + " Restart Selection " + Program.RestartSelectedInputLabel + " Restart Input File Type " + Program.RestartInputFileType);
		DAVectorUtility.SALSAPrint(1, "Files " + config.DistanceMatrixFile + " Comparison " + Program.ComparisonClusterFile + " Restarts " + Program.RestartClusterFile + " "
		    + config.LabelFile);

		if (Program.ComparisonInputFileType == 1)
		{
			Program.CompareSolution = -1;
		}
		// initialize Experiment Number Assignments
		Program.ExperimentNumberAssigments = new int[DAVectorUtility.PointCount_Global];

		//  Read data for LC-MS Clustering Comparisons
		if (Program.CompareSolution > 0)
		{
			if (!Program.DoLCMS)
			{
				DAVectorUtility.printAndThrowRuntimeException("Invalid Compare Solution Option " + Program.CompareSolution);
			}

			int save1 = InputFileType;
			int save2 = SelectedInputLabel;
			InputFileType = ComparisonInputFileType;
			SelectedInputLabel = ComparisonSelectedInputLabel;
			Program.GoldenPeaks = new ArbitraryClustering(DAVectorUtility.PointCount_Global, "Golden Peaks");
			Program.MclustClusters = new ArbitraryClustering(DAVectorUtility.PointCount_Global, "Mclust");
			Program.MedeaClusters = new ArbitraryClustering(DAVectorUtility.PointCount_Global, "Medea");
			Program.OurClusters = new ArbitraryClustering(DAVectorUtility.PointCount_Global, "DAVectorSponge");
			GoldenExamination.GoldenID = new int[DAVectorUtility.PointCount_Global];
			GoldenExamination.GoldenLabel = new String[DAVectorUtility.PointCount_Global];
			GoldenExamination.PeakPosition = new double[DAVectorUtility.PointCount_Global][];
			for (int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
			{
				GoldenExamination.PeakPosition[GlobalPointIndex] = new double[Program.ParameterVectorDimension];
			}
			ReadDataTiming.startTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			DAVectorReadData.ReadLabelsFromFile(ComparisonClusterFile);
			ReadDataTiming.endTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			InputFileType = save1;
			SelectedInputLabel = save2;
		}else{
			//NOTE: This code is experiment dependent
			ReadDataTiming.startTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			DAVectorReadData.ReadExperimentNumbers(DistanceMatrixFile);
			ReadDataTiming.endTiming(ReadDataTiming.TimingTask.TOTAL_READ);
		}
		
		// Set up Decomposition of USED points
		DAVectorParallelism.SetParallelDecomposition();
		if (Program.Replicate > 1)
		{
			DAVectorUtility.SALSAPrint(1, "Replicate " + Program.Replicate);
			DAVectorUtility.PointCount_Global *= Program.Replicate;
			DAVectorUtility.PointCount_Largest *= Program.Replicate;
			DAVectorUtility.PointStart_Process *= Program.Replicate;
			DAVectorUtility.PointCount_Process *= Program.Replicate;
			for (int ProcessIndex = 0; ProcessIndex < DAVectorUtility.MPI_Size; ProcessIndex++)
			{
				DAVectorUtility.PointsperProcess[ProcessIndex] *= Program.Replicate;
				for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
				{
					DAVectorUtility.PointsperThreadperProcess[ProcessIndex][ThreadIndex] *= Program.Replicate;
				}
			}
			for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
			{
				DAVectorUtility.StartPointperThread[ThreadIndex] *= Program.Replicate;
				DAVectorUtility.PointsperThread[ThreadIndex] *= Program.Replicate;
			}
		}

		// Setup PointPosition and related arrays
		Program.PointOriginalIndex = new int[DAVectorUtility.PointCount_Process]; // Point Index on File
		Program.PointOriginalExperimentNumber = new int[DAVectorUtility.PointCount_Process]; // Point Index on File
		Program.PointLabel = new int[DAVectorUtility.PointCount_Process];
		Program.PointPosition = new double[DAVectorUtility.PointCount_Process][];
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    Program.PointPosition[alpha] = new double[Program.ParameterVectorDimension];
                }

            });// End loop initialing Point dependent quantities
        });

        Program.ClusterAssignments = new int[DAVectorUtility.PointCount_Global];

		Program.ClusterLimitforDistribution = Math.max(Program.ClusterLimitforDistribution, 4 * DAVectorUtility.MPI_Size);

		SetupTiming.endTiming(SetupTiming.TimingTask.SETUP);

		// Initial Processing Complete
        // Note - MPI Call - Barrier
        if (DAVectorUtility.MPI_Size > 1){
            DAVectorUtility.mpiOps.barrier(); // Make certain all processes have processed original data before writing updated
        }
		//  read data into memory
		// Timing Section 1
		SectionTiming.startTiming(SectionTiming.TimingTask.SC1);
		if (Program.DoKmeans)
		{
			GeneralTiming.startTiming(GeneralTiming.TimingTask.KMEANS);
			KmeansTriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_Kmeans, Program.MaxClusterLBsperPoint_Kmeans, Program.MaxCentersperCenter_Kmeans, Program.TriangleInequality_Delta1_old_Kmeans, Program.TriangleInequality_Delta1_current_Kmeans, Program.OldCenterOption_Kmeans, Program.DoBackwardFacingTests_Kmeans);
			Kmeans.InitializeKmeans(Program.PointPosition, Program.config.DistanceMatrixFile,
			    Program.ClusterIndexonInputLine, Program.FirstClusterValue, Program.StartPointPositiononInputLine, Program.InitialNcent, Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_Kmeans,
			    Program.ParameterVectorDimension, Program.KmeansCenterChangeStop, Program.KmeansIterationLimit);
			GeneralTiming.endTiming(GeneralTiming.TimingTask.KMEANS);
		}
		else if (Program.DoLCMS)
		{
			ReadDataTiming.startTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			GeneralTiming.startTiming(GeneralTiming.TimingTask.LCMS);
			DAVectorReadData.ReadDataFromFile(config.DistanceMatrixFile, 0);
			ReadDataTiming.endTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			GeneralTiming.endTiming(GeneralTiming.TimingTask.LCMS);
		}
		else
		{
			ReadDataTiming.startTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			DAVectorReadData.ReadDataFromFile(config.DistanceMatrixFile,
			                                  Program.ClusterIndexonInputLine, null,
			                                  Program.StartPointPositiononInputLine);
			ReadDataTiming.endTiming(ReadDataTiming.TimingTask.TOTAL_READ);
		}
		// Program.ClusterIndexonInputLine MUST be NULL

		//  Read 3D Data for Plotviz
		if (config.LabelFile.length() <= 0)
		{
			Program.RW3DData = -1;
		}
		if (Program.RW3DData > 0)
		{
			DAVectorUtility.SALSAPrint(0, "3D Data " + config.LabelFile);
			Program.FullPoint3DPosition = new double[DAVectorUtility.PointCount_Global][];
			for (int alpha = 0; alpha < DAVectorUtility.PointCount_Global; alpha++)
			{
				Program.FullPoint3DPosition[alpha] = new double[Program.RW3DData];
			}
			if (DAVectorUtility.MPI_Rank == 0)
			{
				ReadDataTiming.startTiming(ReadDataTiming.TimingTask.TOTAL_READ);
				DAVectorReadData.Read3DDataFromFile(config.LabelFile, Program.RW3DData, 1);
				ReadDataTiming.endTiming(ReadDataTiming.TimingTask.TOTAL_READ);
			}
		}
		SectionTiming.endTiming(SectionTiming.TimingTask.SC1);
		//  Set up Timing
		int nonMPITimings = 14;
		DAVectorUtility.InitializeTiming(13 + nonMPITimings);
		DAVectorUtility.SetUpMPISubTimers(nonMPITimings, "");
		DAVectorUtility.SetUpSubTimer(0, "Decide Splitting");
		DAVectorUtility.SetUpSubTimer(1, "Validity");
		DAVectorUtility.SetUpSubTimer(2, "EMIterate-1");
		DAVectorUtility.SetUpSubTimer(3, "EMIterate-2");
		DAVectorUtility.SetUpSubTimer(4, "ClusterAvgs");
		DAVectorUtility.SetUpSubTimer(5, "Do Splitting");
		DAVectorUtility.SetUpSubTimer(6, "Major Synch - 1");
		DAVectorUtility.SetUpSubTimer(7, "Minor Synch");
		DAVectorUtility.SetUpSubTimer(8, "Final Status");
		DAVectorUtility.SetUpSubTimer(9, "EMIterate-3");
		DAVectorUtility.SetUpSubTimer(10, "EMIterate-4");
		DAVectorUtility.SetUpSubTimer(11, "Major Synch - 2");
		DAVectorUtility.SetUpSubTimer(12, "Major Synch - 3");
		DAVectorUtility.SetUpSubTimer(13, "Redo Clusters per Point");
		DAVectorUtility.TimingOutputOrder[0] = 0;
		DAVectorUtility.TimingOutputOrder[1] = 5;
		DAVectorUtility.TimingOutputOrder[2] = 1;
		DAVectorUtility.TimingOutputOrder[3] = 13;
		DAVectorUtility.TimingOutputOrder[4] = 4;
		DAVectorUtility.TimingOutputOrder[5] = 8;
		DAVectorUtility.TimingOutputOrder[6] = 2;
		DAVectorUtility.TimingOutputOrder[7] = 3;
		DAVectorUtility.TimingOutputOrder[8] = 9;
		DAVectorUtility.TimingOutputOrder[9] = 10;
		DAVectorUtility.TimingOutputOrder[10] = 7;
		DAVectorUtility.TimingOutputOrder[11] = 6;
		DAVectorUtility.TimingOutputOrder[12] = 11;
		DAVectorUtility.TimingOutputOrder[13] = 12;
		if ((!Program.DoLCMS) && (!Program.DoKmeans))
		{
			DAVectorUtility.SetUpSubTimer(6, "Full Triangle Ineq");
			DAVectorUtility.SetUpSubTimer(11, "CenterFacing");
			DAVectorUtility.SetUpSubTimer(12, "Point LB");
		}
		if (Program.DoKmeans)
		{
			DAVectorUtility.SetUpSubTimer(0, "Pure Kmeans");
			DAVectorUtility.SetUpSubTimer(1, "End Iteration");
			DAVectorUtility.SetUpSubTimer(2, "Center Facing");
			DAVectorUtility.SetUpSubTimer(3, "Point LB");
			DAVectorUtility.SetUpSubTimer(4, "Cluster Centers");
			DAVectorUtility.SetUpSubTimer(5, "Individual Centers");
			DAVectorUtility.SetUpSubTimer(6, "Center Centers");
			DAVectorUtility.SetUpSubTimer(7, "Update Centers");
			for (int timingloop = 9; timingloop < nonMPITimings; timingloop++)
			{
				DAVectorUtility.SubTimingEnable[timingloop] = false;
			}
			for (int timingloop = 0; timingloop < nonMPITimings; timingloop++)
			{
				DAVectorUtility.TimingOutputOrder[timingloop] = timingloop;
			}
		}


		// Timing Section 2
		SectionTiming.startTiming(SectionTiming.TimingTask.SC2);
		//  Set up basic clusters
		ClusteringSolution.SetParameters(DAVectorUtility.PointCount_Process, Program.maxNcentCreated, Program.maxNcentTOTAL, Program.maxNcentperNode, Program.cachelinesize, Program.targetNcentperPoint, Program.targetMinimumNcentperPoint, Program.maxNcentperPoint, Program.ExpArgumentCut2);

		//  Set up Distributed Clusters
		DistributedClusteringSolution DistributedSetup;
		if (Program.DoLCMS)
		{
			GeneralTiming.startTiming(GeneralTiming.TimingTask.LCMS);
			DistributedSetup = new DistributedClusteringSolution(MaxMPITransportBuffer, MaxTransportedClusterStorage, MaxNumberAccumulationsperNode, MaxDoubleComponents, MaxIntegerComponents);
			GeneralTiming.endTiming(GeneralTiming.TimingTask.LCMS);
		}
		DAVectorUtility.SALSAPrint(0, "Setup Finished");

		//  Do Clustering
		VectorAnnealIterate RunVectorSpongeDA;
		ControlKmeans RunKmeansClustering;
		if (DoKmeans)
		{
			GeneralTiming.startTiming(GeneralTiming.TimingTask.KMEANS);
			RunKmeansClustering = new ControlKmeans();
			GeneralTiming.endTiming(GeneralTiming.TimingTask.KMEANS);
		}
		else
		{
			GeneralTiming.startTiming(GeneralTiming.TimingTask.DA);
			RunVectorSpongeDA = new VectorAnnealIterate();
			RunVectorSpongeDA.ControlVectorSpongeDA();
			GeneralTiming.endTiming(GeneralTiming.TimingTask.DA);

		}
		Program.ActualEndTemperatureafterconverging = ParallelClustering.runningSolution.Temperature;

		//  End Timing
		DAVectorUtility.EndTiming();

		// Calculate Cluster Statistics
		DAVectorUtility.StartSubTimer(8);

		// Calculate Occupation Counts and other statistics for LCMS
		//  This can alter "Sponge Confused Points"
		if (!Program.DoKmeans)
		{
            GeneralTiming.startTiming(GeneralTiming.TimingTask.LCMS);
            ParallelClustering.runningSolution.FindOccupationCounts();
            GeneralTiming.endTiming(GeneralTiming.TimingTask.LCMS);
        }

		if (Program.ClusterCountOutput >= 0 && (!Program.DoKmeans))
		{
			VectorAnnealIterate.OutputClusteringResults("Final");
		}

		if (Program.DoBasicDA && (Program.RW3DData > 0))
		{
			VectorAnnealIterate.Output3DClusterLabels("DA");
		}
		if (Program.DoKmeans)
		{
            GeneralTiming.startTiming(GeneralTiming.TimingTask.KMEANS);
            double[][] ClusterCenters = ControlKmeans.CaptureKmeans(Program.ClusterAssignments);
			if (Program.ClusterCountOutput >= 0)
			{
				VectorAnnealIterate.SimpleOutputClusteringResults("Kmeans", ClusterCenters);
				DAVectorUtility.SALSAPrint(0, "End Output of Full Clusters");
			}
			if (Program.RW3DData > 0)
			{
				VectorAnnealIterate.Output3DClusterLabels("Kmeans");
				DAVectorUtility.SALSAPrint(0, "End Output of 3D Clusters");
			}
            GeneralTiming.endTiming(GeneralTiming.TimingTask.KMEANS);
        }
		if (Program.DoLCMS)
		{
			// Set Cluster Quality Parameters

			Program.HistogramDistancesMax = config.HistogramDistancesMax;
			Program.HistogramOccupationMax = config.HistogramOccupationMax	;
			GeneralTiming.startTiming(GeneralTiming.TimingTask.LCMS);
            VectorAnnealIterate.CalculateClusterStatus();
            GeneralTiming.endTiming(GeneralTiming.TimingTask.LCMS);
        }


		SectionTiming.endTiming(SectionTiming.TimingTask.SC2);
		// Output Results
		//Section Timing 3
		SectionTiming.startTiming(SectionTiming.TimingTask.SC3);
		String nextline = "\nNode 0 Center Averages (Counts) [Distce] ";
		int Totnumber = ClusteringSolution.NumberLocalActiveClusters;
		if (Totnumber > Program.MaxNumberClustersToPrint)
		{
			nextline += Totnumber + " Truncated at " + Program.MaxNumberClustersToPrint + " ";
			Totnumber = Program.MaxNumberClustersToPrint;
		}
		for (int ActiveClusterIndex = 0; ActiveClusterIndex < Totnumber; ActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
			double meandist = ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex];
			String cees = "";
			if (!Program.DoKmeans)
			{
				cees = String.format("%1$3.2f", ParallelClustering.runningSolution.C_k_[RealClusterIndex]) + " - ";
			}
			nextline += RealClusterIndex + " " + cees + ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex] + " [Wid " + String.format("%1$2.1f", meandist) + "]";
			if (!Program.DoKmeans)
			{
				nextline += "[Frz " + String.format("%1$3.2E", ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex]) + "]";
			}
			nextline += " * ";
		}
		boolean save = DAVectorUtility.ConsoleDebugOutput;
		if (Totnumber > 100)
		{
			DAVectorUtility.ConsoleDebugOutput = false;
		}
		DAVectorUtility.SALSAPrint(0, nextline);
		DAVectorUtility.ConsoleDebugOutput = save;

		if (!Program.DoKmeans)
		{
			double tmp1 = 0.0;
			if (Program.NumberMdiffSums > 0)
			{
				Program.MdiffSum = Program.MdiffSum / Program.NumberMdiffSums;
			}
			if (Program.NumberIterationSteps > 0)
			{
				tmp1 = (double) Program.IterationsperStepSum / (double) Program.NumberIterationSteps;
			}

			DAVectorUtility.SALSAPrint(0, "\nT " + String.format("%1$6.5f", ParallelClustering.runningSolution.Temperature) + " Cluster " + ParallelClustering.runningSolution.Ncent_Global + " Iter " + VectorAnnealIterate.EMIterationCount + " Extra Iter " + VectorAnnealIterate.Extra_EMIterationCount + " Average Mdiff " + String.format("%1$6.5f", Program.MdiffSum) + " Number of Steps " + Program.NumberIterationSteps + " Iterations per Step " + String.format("%1$3.2f", tmp1));
		}

		DAVectorUtility.SALSAPrint(0, "\n" + DAVectorUtility.PatternLabel);
		DAVectorUtility.SALSAPrint(0, "Labels File: " + config.ClusterFile);
		DAVectorUtility.SALSAPrint(0, "Timing Output: " + config.TimingFile);
		String message = " Charge Value ";
		if (InputFileType == 1)
		{
			message = " Selected with Cluster ";
		}
		if (Program.DoKmeans)
		{
			message = " Selection ";
		}
		DAVectorUtility.SALSAPrint(0, "Data Points " + DAVectorUtility.PointCount_Global + message + Program.SelectedInputLabel);
		DAVectorUtility.SALSAPrint(0, "Vector Dimension: " + Program.ParameterVectorDimension);
		DAVectorUtility.SALSAPrint(0, "Continuous Clustering: " + Program.ContinuousClustering);

		String sigmamethodString = "Sigma Method: " + Program.SigmaMethod;
		if (Program.SigmaMethod > 0)
		{
			sigmamethodString += " Parameters:";
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				sigmamethodString += " " + String.format("%1$5.4E", Program.SigmaVectorParameters_i_[VectorIndex]);
			}
		}
		DAVectorUtility.SALSAPrint(0, sigmamethodString);
		if (Program.SigmaMethod == 3)
		{
			DAVectorUtility.SALSAPrint(0, "Initial Sigma[0] " + String.format("%1$5.4E", Program.InitialSigma0) + " Final Sigma[0] = " + String.format("%1$5.4E", Program.FinalTargetSigma0) + " *Center[0] at Temperature " + String.format("%1$5.4f", Program.FinalTargetTemperature));
		}

		DAVectorUtility.SALSAPrint(0, "Initial Number of Centers: " + Program.InitialNcent);
		DAVectorUtility.SALSAPrint(0, "Final Number of Centers: " + ParallelClustering.runningSolution.Ncent_Global);
		DAVectorUtility.SALSAPrint(0, "Maximum Number of Centers per node: " + Program.maxNcentperNode);
		DAVectorUtility.SALSAPrint(0, "Target Maximum number of cluster centers for each point (includes Sponge): " + Program.targetNcentperPoint);
		DAVectorUtility.SALSAPrint(0, "Actual Maximum number of cluster centers for each point (includes Sponge): " + Program.maxNcentperPoint);

		if (!Program.DoKmeans)
		{
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - int - max
//                ActualMaxMPITransportBuffer = DAVectorUtility.MPI_communicator.<Integer>Allreduce(ActualMaxMPITransportBuffer, Operation<Integer>.Max);
                ActualMaxMPITransportBuffer = DAVectorUtility.mpiOps.allReduce(ActualMaxMPITransportBuffer, MPI.MAX);

                // Note - MPI Call - Allreduce - int - max
//                ActualMaxNumberAccumulationsperNode = DAVectorUtility.MPI_communicator.<Integer>Allreduce(ActualMaxNumberAccumulationsperNode, Operation<Integer>.Max);
                ActualMaxNumberAccumulationsperNode = DAVectorUtility.mpiOps.allReduce(ActualMaxNumberAccumulationsperNode, MPI.MAX);

                // Note - MPI Call - Allreduce - int - max
//                ActualMaxTransportedClusterStorage = DAVectorUtility.MPI_communicator.<Integer>Allreduce(ActualMaxTransportedClusterStorage, Operation<Integer>.Max);
                ActualMaxTransportedClusterStorage = DAVectorUtility.mpiOps.allReduce(ActualMaxTransportedClusterStorage, MPI.MAX);
            }
			DAVectorUtility.SALSAPrint(0, "Needed Maximum: MaxMPITransportBuffer " + ActualMaxMPITransportBuffer + " MaxNumberAccumulationsperNode " + ActualMaxNumberAccumulationsperNode + " MaxTransportedClusterStorage " + ActualMaxTransportedClusterStorage);
			int MaxCreatedIndex = ClusteringSolution.CenterMaxforCreatedIndex;
            // Note - MPI Call - Allreduce - int - max
            if (DAVectorUtility.MPI_Size > 1){
//                MaxCreatedIndex = DAVectorUtility.MPI_communicator.<Integer>Allreduce(MaxCreatedIndex, Operation<Integer>.Max);
                MaxCreatedIndex = DAVectorUtility.mpiOps.allReduce(MaxCreatedIndex, MPI.MAX);
            }
			DAVectorUtility.SALSAPrint(0, "Maximum Created Index Space over all nodes " + MaxCreatedIndex + " times for each cluster " + ClusteringSolution.PACKINGMULTIPLIER);

			DAVectorUtility.SALSAPrint(0, "Iterations at the End: " + Program.Iterationatend);
			DAVectorUtility.SALSAPrint(0, "Limit on EM Convergence for each step: " + Program.ConvergenceLoopLimit);
			DAVectorUtility.SALSAPrint(0, "Change in Average M per point: " + String.format("%1$5.4E", Program.Malpha_MaxChange) + " In final Loop " + String.format("%1$5.4E", Program.Malpha_MaxChange1) + " Y Change (final only) " + String.format("%1$5.4E", YChangeSquared));
			DAVectorUtility.SALSAPrint(0, "Freezing Limit for convergence: " + String.format("%1$5.4E", Program.FreezingLimit));
			DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Iteration Code " + String.format("%1$5.4f", Program.ExpArgumentCut1));
			DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Associating Clusters with Points " + String.format("%1$5.4f", Program.ExpArgumentCut2) + " (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut");
			DAVectorUtility.SALSAPrint(0, "Mean Cluster Count per point " + String.format("%1$3.2f", VectorAnnealIterate.MeanClusterCount) + " Pts with Just 1 Cluster " + String.format("%1$5.4E", VectorAnnealIterate.PointswithClusterCount1));

			DAVectorUtility.SALSAPrint(0, "Tmininimum " + String.format("%1$5.4f", Program.Tminimum) + " If Positive this is minmum Temperature, If negative Divide maximum temperature by this to get target minimum");

			DAVectorUtility.SALSAPrint(0, "\nNumber of Clusters Split Simultaneously: " + Program.MaxNumberSplitClusters + " Initially " + SaveSplitNumber);

			DAVectorUtility.SALSAPrint(0, "Do not split Clusters smaller than this: " + String.format("%1$4.3f", Program.ToosmalltoSplit));
			DAVectorUtility.SALSAPrint(0, "Do not split Clusters with Scaled Squared Width Less than this " + String.format("%1$4.3f", Program.MinimumScaledWidthsquaredtosplit));
			DAVectorUtility.SALSAPrint(0, "Wait stages between splits: " + Program.Waititerations);
			DAVectorUtility.SALSAPrint(0, "Initial Cooling Factor in Annealing: " + String.format("%1$8.7f", Program.InitialCoolingFactor1) + " " + String.format("%1$8.7f", Program.InitialCoolingFactor2));
			DAVectorUtility.SALSAPrint(0, "Refined Cooling Factor in Annealing: " + String.format("%1$8.7f", Program.FineCoolingFactor1) + " " + String.format("%1$8.7f", Program.FineCoolingFactor2) + " Switching at " + String.format("%1$4.3f", Program.CoolingTemperatureSwitch));

			int cleanuplength = Program.MagicTemperatures.length;
			message = "Clean up Temperatures ";
			for (int loop = 0; loop < cleanuplength; loop++)
			{
				message += String.format("%1$3.2f", Program.MagicTemperatures[loop]) + " ";
			}
			DAVectorUtility.SALSAPrint(0, message);

			if (ParallelClustering.runningSolution.SpongeCluster >= 0)
			{
				DAVectorUtility.SALSAPrint(0, "\nSponge Cluster: " + ParallelClustering.runningSolution.SpongeCluster);
				if (Program.UseSponge)
				{
					DAVectorUtility.SALSAPrint(0, "Sponge Created Initially");
				}
				DAVectorUtility.SALSAPrint(0, "Final Sponge Factor: " + String.format("%1$3.2f", Program.SpongeFactor) + " Initial " + String.format("%1$3.2f", Program.SpongeFactor1) + " Started at " + String.format("%1$4.3f", Program.SpongeTemperature1) + " Aiming at " + String.format("%1$3.2f", Program.SpongeFactor2) + " at Temperature " + String.format("%1$4.3f", Program.SpongeTemperature2));
				DAVectorUtility.SALSAPrint(0, "Sponge p(k) Option: " + Program.SpongePoption + " Weighting " + String.format("%1$4.3E", Program.SpongePWeight));
				DAVectorUtility.SALSAPrint(0, "Nearby Sponge Point Limit: " + (new Double(Program.NearbySpongePointLimit)));
				if (ParallelClustering.runningSolution.SpongeCluster >= 0)
				{
					DAVectorUtility.SALSAPrint(0, " Points in Sponge Cluster " + String.format("%1$3.2f", ParallelClustering.runningSolution.C_k_[ParallelClustering.runningSolution.SpongeCluster]) + " Occupation Count " + ParallelClustering.runningSolution.OccupationCounts_k_[ParallelClustering.runningSolution.SpongeCluster]);
				}
				DAVectorUtility.SALSAPrint(0, "Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches " + String.format("%1$3.2f", Program.CreateSpongeScaledSquaredWidth));
				if (Program.ActualSpongeTemperature > 0.0)
				{
					double FractionTimeatSponge = Program.TimeatSponge / DAVectorUtility.HPDuration;
					DAVectorUtility.SALSAPrint(0, "Temperature where Sponge Introduced " + String.format("%1$5.4f", Program.ActualSpongeTemperature) + " Width Then " + String.format("%1$4.3f", Program.ActualWidth) + " Time fraction " + String.format("%1$5.4f", FractionTimeatSponge));
				}
			}
		}

		DAVectorUtility.SALSAPrint(0, "\nCluster Deletion Statistics: C_k_ Small " + Program.TotalClustersDeleted_CSmall + " Occ. Count Small " + Program.TotalClustersDeleted_OccCount + " Clusters Close " + Program.TotalClustersDeleted_Close);
		DAVectorUtility.SALSAPrint(0, "Minimum Cluster Count Average Point for Absorbing: " + String.format("%1$3.2f", Program.MinimumCountforCluster_C_k) + " With Sponge " + String.format("%1$3.2f", Program.MinimumCountforCluster_C_kwithSponge));
		DAVectorUtility.SALSAPrint(0, "Minimum Cluster Identified Point Count for Absorbing: " + Program.MinimumCountforCluster_Points);
		DAVectorUtility.SALSAPrint(0, "Minimum Cluster Average Count for considering as nonzero " + String.format("%1$6.5f", Program.CountforCluster_C_ktobezero));
		DAVectorUtility.SALSAPrint(0, "Scaled Squared Distance Cut for Closeness " + String.format("%1$3.2f", Program.ScaledSquaredDistanceatClosenessTest) + " Temperature " + String.format("%1$3.2f", Program.TemperatureforClosenessTest));

		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			double FractionTimeatDist = Program.TimeatDistribution / DAVectorUtility.HPDuration;
			DAVectorUtility.SALSAPrint(0, "\nDistributed Execution entered at Temperature " + String.format("%1$5.4f", Program.ActualTemperatureforDistribution) + " Time(fraction) " + String.format("%1$5.4f", FractionTimeatDist) + " Cut " + String.format("%1$5.4f", Program.TemperatureLimitforDistribution) + " Clusters " + Program.ActualClusterNumberforDistribution + " Cut " + Program.ClusterLimitforDistribution);
			DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Distributed Broadcast " + String.format("%1$5.4f", Program.ExpArgumentCut3) + " (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut -- 0th component only");
			DAVectorUtility.SALSAPrint(0, "Number of Minor Synchronizations " + Program.NumberMinorSynchs + " Number of Major Synchs while Global " + Program.NumberMajorSynchs1 + " Number of Major Synchs while Distributed " + Program.NumberMajorSynchs2);
			double ClustersperPipeline = 0.0;
			if (Program.NumberPipelineSteps > 0)
			{
				ClustersperPipeline = (double) Program.NumberofPipelineClusters / (double) Program.NumberPipelineSteps;
			}
			DAVectorUtility.SALSAPrint(0, " Number of Pipeline steps " + Program.NumberPipelineSteps + " with average cluster load " + String.format("%1$4.3f", ClustersperPipeline) + " In " + Program.NumberPipelineGroups + " sets");
			double AverageClusterCount = (double) Program.CountClusters2 / (double) Program.NumberMajorSynchs2;
			double[] AverageClusterCountvalues = new double[DAVectorUtility.MPI_Size];

            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allgather - double
//                AverageClusterCountvalues = DAVectorUtility.MPI_communicator.<Double>Allgather(AverageClusterCount);
                 DAVectorUtility.mpiOps.allGather(AverageClusterCount, AverageClusterCountvalues);
            } else {
                AverageClusterCountvalues[0] = AverageClusterCount;
            }
			double maxClusterCount = 0.0;
			double minClusterCount = 1.0E+08;
			double avgClusterCount = 0.0;
			for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
			{
				avgClusterCount += AverageClusterCountvalues[nodes];
				if (AverageClusterCountvalues[nodes] > maxClusterCount)
				{
					maxClusterCount = AverageClusterCountvalues[nodes];
				}
				if (AverageClusterCountvalues[nodes] < minClusterCount)
				{
					minClusterCount = AverageClusterCountvalues[nodes];
				}
			}
			avgClusterCount = avgClusterCount / DAVectorUtility.MPI_Size;
			AverageClusterCount = (double) ParallelClustering.runningSolution.Ncent_ThisNode;
            // Note - MPI Call - Allgather - double
            if (DAVectorUtility.MPI_Size > 1){
//                AverageClusterCountvalues = DAVectorUtility.MPI_communicator.<Double>Allgather(AverageClusterCount);
                DAVectorUtility.mpiOps.allGather(AverageClusterCount,AverageClusterCountvalues);
            } else {
                AverageClusterCountvalues[0] = AverageClusterCount;
            }
			double maxClusterCountFinal = 0.0;
			double minClusterCountFinal = 1.0E+08;
			double avgClusterCountFinal = 0.0;
			for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
			{
				avgClusterCountFinal += AverageClusterCountvalues[nodes];
				if (AverageClusterCountvalues[nodes] > maxClusterCountFinal)
				{
					maxClusterCountFinal = AverageClusterCountvalues[nodes];
				}
				if (AverageClusterCountvalues[nodes] < minClusterCountFinal)
				{
					minClusterCountFinal = AverageClusterCountvalues[nodes];
				}
			}
			avgClusterCountFinal = avgClusterCountFinal / DAVectorUtility.MPI_Size;
			DAVectorUtility.SALSAPrint(0, "Clusters per node Average " + String.format("%1$3.2f", avgClusterCount) + " Min " + String.format("%1$3.2f", minClusterCount) + " Max " + String.format("%1$3.2f", maxClusterCount));
			DAVectorUtility.SALSAPrint(0, "Clusters per node Final " + String.format("%1$3.2f", avgClusterCountFinal) + " Min " + String.format("%1$3.2f", minClusterCountFinal) + " Max " + String.format("%1$3.2f", maxClusterCountFinal));
		}

		DAVectorUtility.SALSAPrint(0, "\nActual Starting Temperature " + String.format("%1$5.4f", Program.ActualStartTemperature));
		DAVectorUtility.SALSAPrint(0, "Final Temperature after Converging " + String.format("%1$5.4f", Program.ActualEndTemperatureafterconverging));
		double FractionTimeatSplitStop = Program.TimeatSplittingStop / DAVectorUtility.HPDuration;
		DAVectorUtility.SALSAPrint(0, "Temperature where we decided to stop " + String.format("%1$5.4f", Program.ActualEndTemperature) + " Fractional Time " + String.format("%1$5.4f", FractionTimeatSplitStop));
		DAVectorUtility.SALSAPrint(0, "Temperature where we aimed to stop at " + String.format("%1$5.4f", Program.TargetEndTemperature));
		if (!Program.DoKmeans)
		{
			double TotalCalcs = Program.SumUsefulCalcs + Program.SumUselessCalcs + Program.SumIgnoredCalcs;
			DAVectorUtility.SALSAPrint(0, "\nNumber of distance calcs in EMIterate " + String.format("%1$5.4E", TotalCalcs) + " Useful " + String.format("%1$5.4f", Program.SumUsefulCalcs / TotalCalcs) + " Useless " + String.format("%1$5.4f", Program.SumUselessCalcs / TotalCalcs) + " Ignored " + String.format("%1$5.4f", Program.SumIgnoredCalcs / TotalCalcs));
			DAVectorUtility.SALSAPrint(0, "Scalar Products in Eigenvector part " + String.format("%1$5.4E", Program.SumEigenSPCalcs) + " Sum Useful + Useless " + String.format("%1$5.4E", Program.SumUsefulCalcs + Program.SumUselessCalcs));

			double MIterationAverage = Program.AccumulateMvalues / Program.NumberMsuccesses;
			DAVectorUtility.SALSAPrint(0, "Mchange Test Failures " + String.format("%1$5.4E", Program.NumberMfailures) + " Successes " + String.format("%1$5.4E", Program.NumberMsuccesses) + " Average Iterations " + String.format("%1$3.2f", MIterationAverage) + " Max " + Program.ConvergenceLoopLimit);
			double YIterationAverage = Program.AccumulateYvalues / Program.NumberYsuccesses;
			DAVectorUtility.SALSAPrint(0, "Ychange Test Failures " + String.format("%1$5.4E", Program.NumberYfailures) + " Successes " + String.format("%1$5.4E", Program.NumberYsuccesses) + " Average Iterations " + String.format("%1$3.2f", YIterationAverage) + " Max " + Program.ConvergenceLoopLimit);
		}
		if (Program.UseTriangleInequality_DA > 0)
		{
			double weighting = 1.0 / Program.ClustersperCenterDiagnostics[0];
			for (int DiagnosticLoop = 1; DiagnosticLoop < 8; DiagnosticLoop++)
			{
				Program.ClustersperCenterDiagnostics[DiagnosticLoop] *= weighting;
			}
			DAVectorUtility.SALSAPrint(0, "Dynamic Clusters in Triangle Inequality Failures " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[1]) + " Succeses " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[2]) + " Very Small " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[3]) + " Small " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[4]) + " Too Big " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[5]) + " Average reset M " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[6]) + " Average Changes in Clusters per Point " + String.format("%1$5.4f", Program.ClustersperCenterDiagnostics[7]));
		}

		DAVectorUtility.SALSAPrint(0, "\nPoints " + DAVectorUtility.PointCount_Global + " Selection " + Program.SelectedInputLabel + " Number of Clusters " + ParallelClustering.runningSolution.Ncent_Global + " Temperature Steps " + Program.NumberTemperatureSteps + " " + Program.FinalReason);

		//  Cluster Widths
		String widthmessage = "";
		if (!Program.CalculateIndividualWidths)
		{
			widthmessage = " Average Width " + String.format("%1$4.3E", ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
		}
		else
		{
			widthmessage = " Average Widths excluding Sponge Points (These are in Hamiltonian value) ";
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				widthmessage += String.format("%1$4.3E", ParallelClustering.runningSolution.AverageWidth[VectorIndex]) + " ";
			}
		}
		DAVectorUtility.SALSAPrint(0, widthmessage + " These are squared and divided by squares of appropriate sigmas with Hamiltonian " + String.format("%1$5.4E", ParallelClustering.runningSolution.PairwiseHammy));

		if (Program.DoBasicDA)
		{ //  General Statistics
			ClusterQuality.CaculateTemperatureClusterCountPlot();
			ClusteringSolution.TotalClusterSummary.HistogramClusterProperties();
			if (Program.UseTriangleInequality_DA > 0)
			{
				DATriangleInequality.PrintDiagnostics();
			}
		}

		/* Compute the duration between the initial and the end time ignoring print out. */
		DAVectorUtility.StopSubTimer(8);
		DAVectorUtility.SALSAPrint(0, "\nTotal Time excluding I/O  " + DAVectorUtility.formatElapsedMillis(DAVectorUtility.mainDuration) + " " + String.format("%1$5.4E", DAVectorUtility.HPDuration * .001));


        nextline = "Partial Times ";
		for (int jtimer = 0; jtimer < DAVectorUtility.NumberofSubTimings; jtimer++)
		{
			int itimer = DAVectorUtility.TimingOutputOrder[jtimer];
			if (itimer == nonMPITimings)
			{
				nextline += "\n\n";
			}
			if (DAVectorUtility.SubTimingEnable[itimer])
			{
				double tmp = DAVectorUtility.SubDurations[itimer] / DAVectorUtility.HPDuration;
				nextline += DAVectorUtility.SubTimingNames[itimer] + " Time=" + String.format("%1$5.4E", DAVectorUtility.SubDurations[itimer] * .001) + " #=" + DAVectorUtility.SubTimingCalls[itimer] + " " + String.format("%1$5.4f", tmp) + " ** ";
			}
		}
		DAVectorUtility.SALSAPrint(0, nextline);

		//  Add a Fast K means step using triangle inequality speed up
		if (Program.DoBasicDA && Program.FollowUpKmeans && (ParallelClustering.runningSolution.SpongeCluster < 0) && !ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			Program.DoKmeans = true;

			Program.TemperatureLimitforDistribution = -1.0;
			Program.ClusterLimitforDistribution = -1;

			Program.CalculateEigenvaluesfromMatrix = false;
			Program.UseSponge = false;
			Program.ContinuousClustering = false;
			Program.SigmaMethod = 0;
			Program.targetNcentperPoint = ParallelClustering.runningSolution.Ncent_Global;
			Program.maxNcentperPoint = ParallelClustering.runningSolution.Ncent_Global;
			Program.targetMinimumNcentperPoint = 2;

			Program.InitialNcent = ParallelClustering.runningSolution.Ncent_Global;
			Program.maxNcentperNode = ParallelClustering.runningSolution.Ncent_Global;
			Program.maxNcentTOTAL = ParallelClustering.runningSolution.Ncent_Global;
			Program.maxNcentTOTALforParallelism_Kmeans = 50;
			Program.maxNcentCreated = 200;
			Program.NumberDataPoints = 200000;

			Program.UseTriangleInequality_Kmeans = 0; // 0 is Pure K means
			Program.MaxClusterLBsperPoint_Kmeans = ParallelClustering.runningSolution.Ncent_Global;
			Program.MaxCentersperCenter_Kmeans = ParallelClustering.runningSolution.Ncent_Global;

			Program.TriangleInequality_Delta1_old_Kmeans = 0.2;
			Program.TriangleInequality_Delta1_current_Kmeans = 0.2;
			Program.KmeansCenterChangeStop = 0.00001;
			Program.KmeansIterationLimit = 1000;

			KmeansTriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_Kmeans, Program.MaxClusterLBsperPoint_Kmeans, Program.MaxCentersperCenter_Kmeans, Program.TriangleInequality_Delta1_old_Kmeans, Program.TriangleInequality_Delta1_current_Kmeans, Program.OldCenterOption_Kmeans, Program.DoBackwardFacingTests_Kmeans);
			Kmeans.InitializeKmeans(Program.PointPosition, "", Program.ClusterIndexonInputLine, Program.FirstClusterValue, Program.StartPointPositiononInputLine, Program.InitialNcent, Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_Kmeans, Program.ParameterVectorDimension, Program.KmeansCenterChangeStop, Program.KmeansIterationLimit);
			Kmeans.SetupKmeans(Program.ClusterAssignments);
			RunKmeansClustering = new ControlKmeans();

			double[][] ClusterCenters = ControlKmeans.CaptureKmeans(Program.ClusterAssignments);
			if (Program.ClusterCountOutput >= 0)
			{
				VectorAnnealIterate.SimpleOutputClusteringResults("Kmeans", ClusterCenters);
			}
			if (Program.RW3DData > 0)
			{
				VectorAnnealIterate.Output3DClusterLabels("Kmeans");
			}

		}

		SectionTiming.endTiming(SectionTiming.TimingTask.SC3);
		GeneralTiming.endTiming(GeneralTiming.TimingTask.TOTAL);

		if (DAVectorUtility.MPI_Rank == 0)
		{
			DAVectorUtility.writeClusterResults(config.SummaryFile, DAVectorUtility.CosmicOutput);
			WriteTimingFile(config.TimingFile, DAVectorUtility.mainDuration, DAVectorUtility.HPDuration,
                    DAVectorUtility.ThreadCount, DAVectorUtility.MPIperNodeCount, DAVectorUtility.NodeCount,
                    DAVectorUtility.PointCount_Process, DAVectorUtility.PointCount_Global, Program.maxNcentperNode,
                    DAVectorUtility.SubDurations[0] * 0.001, DAVectorUtility.SubDurations[0] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[1] * 0.001, DAVectorUtility.SubDurations[1] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[2] * 0.001, DAVectorUtility.SubDurations[2] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[3] * 0.001, DAVectorUtility.SubDurations[3] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[4] * 0.001, DAVectorUtility.SubDurations[4] / DAVectorUtility.HPDuration,
                    config.DistanceMatrixFile, new Date(), MPI.getProcessorName());
            WriteTimingFileNew(config.TimingFile.replace("timing.txt","timingNew.txt"));
		}

        try {
            // Finalize MPI and threads
            DAVectorParallelism.TearDownParallelism();
        } catch (MPIException e) {
            DAVectorUtility.printException(e);
        }
    }

	public static void ReadControlFile(CommandLine cmd)
	{
        config = ConfigurationMgr.LoadConfiguration(cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).daVectorSpongeSection;
        Program.ControlFileName = cmd.getOptionValue(Constants.CMD_OPTION_LONG_C);
        DAVectorUtility.NodeCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        DAVectorUtility.ThreadCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));

        // Override section's node and thread counts with command line values if different
        if (config.getNodeCount() != DAVectorUtility.NodeCount) {
            config.setNodeCount(DAVectorUtility.NodeCount);
        }
        if (config.getThreadCount() != DAVectorUtility.ThreadCount) {
            config.setThreadCount(DAVectorUtility.ThreadCount);
        }


        Program.ClusterFile = config.ClusterFile;
        Program.DistanceMatrixFile = config.DistanceMatrixFile;
        Program.LabelFile= config.LabelFile	;
        Program.ComparisonClusterFile = config.ComparisonClusterFile;
        Program.RestartClusterFile = config.RestartClusterFile;
        Program.TimingFile = config.TimingFile;
        Program.SummaryFile = config.SummaryFile;

        Program.UseSponge = config.UseSponge;
        Program.SpongeFactor1 = config.SpongeFactor1;
        Program.SpongeFactor2 = config.SpongeFactor2;
        Program.SpongePoption = config.SpongePOption;
        Program.SpongePWeight = config.SpongePWeight;
        Program.CreateSpongeScaledSquaredWidth = config.CreateSpongeScaledSquaredWidth;
        Program.ContinuousClustering = config.ContinuousClustering;
        Program.ParameterVectorDimension = config.ParameterVectorDimension;

        Program.SpongeTemperature1 = config.SpongeTemperature1;
        Program.SpongeTemperature2 = config.SpongeTemperature2;
        Program.RestartTemperature = config.RestartTemperature;

        Program.NumberDataPoints = config.NumberDataPoints;
        Program.SelectedInputLabel = config.SelectedInputLabel;
        Program.InputFileType = config.InputFileType;
        Program.Replicate = config.Replicate;

        Program.CompareSolution = config.CompareSolution;
        Program.ComparisonInputFileType = config.ComparisonInputFileType;
        Program.ComparisonSelectedInputLabel = config.ComparisonSelectedInputLabel;
        Program.RestartSelectedInputLabel = config.RestartSelectedInputLabel;
        Program.RestartInputFileType = config.RestartInputFileType;

        Program.SigmaMethod = config.SigmaMethod;
        Program.SigmaVectorParameters_i_ = config.SigmaVectorParameters_i;
        Program.FinalTargetTemperature = config.FinalTargetTemperature;
        Program.FinalTargetSigma0 = config.FinalTargetSigma0;
        Program.InitialSigma0 = config.InitialSigma0;

        Program.ClusterCountOutput = config.ClusterCountOutput;
        Program.NumberNearbyClusters = config.NumberNearbyClusters;
        Program.NearbySpongePointLimit = config.NearbySpongePointLimit;

        Program.ProcessingOption = config.ProcessingOption;

        Program.cachelinesize = config.CacheLineSize;
        Program.ClusterPrintNumber = config.ClusterPrintNumber;
        Program.PrintInterval = config.PrintInterval;
        Program.RemovalDiagnosticPrint = config.RemovalDiagnosticPrint;
        Program.MagicTemperatures = config.MagicTemperatures;
        Program.magicindex = config.MagicIndex;

        Program.maxNcentperNode = config.MaxNcentPerNode;
        Program.maxNcentTOTAL = config.MaxNcentTotal;
        Program.maxNcentCreated = config.maxNcentCreated;
        Program.targetNcentperPoint = config.TargetNcentPerPoint;
        Program.targetMinimumNcentperPoint = config.TargetMinimumNcentPerPoint;
        Program.maxNcentperPoint = config.MaxNcentPerPoint;

        Program.MaxIntegerComponents = config.MaxIntegerComponents;
        Program.MaxDoubleComponents = config.MaxDoubleComponents;
        Program.MaxMPITransportBuffer = config.MaxMPITransportBuffer;
        Program.MaxNumberAccumulationsperNode = config.MaxNumberAccumulationsPerNode;
        Program.MaxTransportedClusterStorage = config.MaxTransportedClusterStorage;

        Program.ExpArgumentCut1 = config.ExpArgumentCut1;
        Program.ExpArgumentCut2 = config.ExpArgumentCut2;
        Program.ExpArgumentCut3 = config.ExpArgumentCut3;
        Program.Tminimum = config.Tminimum;

        Program.InitialNcent = config.InitialNcent;
        Program.MinimumCountforCluster_C_k = config.MinimumCountForClusterCk;
        Program.MinimumCountforCluster_C_kwithSponge = config.MinimumCountForClusterCkWithSponge;
        Program.MinimumCountforCluster_Points = config.MinimuCountForClusterPoints;
        Program.CountforCluster_C_ktobezero = config.CountForClusterCkToBeZero;
        Program.AddSpongeScaledWidthSquared = config.AddSpongeScaledWidthSquared;

        Program.InitialCoolingFactor = config.InitialCoolingFactor;
        Program.InitialCoolingFactor1 = config.InitialCoolingFactor1;
        Program.InitialCoolingFactor2 = config.InitialCoolingFactor2;
        Program.FineCoolingFactor = config.FineCoolingFactor;
        Program.FineCoolingFactor1 = config.FineCoolingFactor1;
        Program.FineCoolingFactor2 = config.FineCoolingFactor2;
        Program.CoolingTemperatureSwitch = config.CoolingTemperatureSwitch;
        Program.Waititerations = config.WaitIterations;
        Program.Waititerations_Converge = config.WaitIterationsConverge;

        Program.Iterationatend = config.IterationAtEnd;
        Program.ConvergenceLoopLimit = config.ConvergenceLoopLimit;

        Program.FreezingLimit = config.FreezingLimit;
        Program.Malpha_MaxChange = config.MalphaMaxChange;
        Program.Malpha_MaxChange1 = config.MalphaMaxChange1;
        Program.MaxNumberSplitClusters = config.MaxNumberSplitClusters;
        Program.ConvergeIntermediateClusters = config.ConvergeIntermediateClusters;
        Program.ToosmalltoSplit = config.TooSmallToSplit;
        Program.MinimumScaledWidthsquaredtosplit = config.MinimumScaledWidthSquaredToSplit;
        Program.ScaledSquaredDistanceatClosenessTest = config.ScaledSquaredDistanceAtClosenessTest;

        Program.ClusterLimitforDistribution = config.ClusterLimitForDistribution;
        Program.TemperatureLimitforDistribution = config.TemperatureLimitForDistribution;
        Program.TemperatureforClosenessTest = config.TemperatureForClosenessTest;

        DAVectorUtility.DebugPrintOption = config.DebugPrintOption;
        DAVectorUtility.ConsoleDebugOutput = config.ConsoleDebugOutput;
	}

	public static void WriteControlFile()
	{

	    config.SpongeTemperature1 = Program.SpongeTemperature1;
	    config.SpongeTemperature2 = Program.SpongeTemperature2;
	    config.RestartTemperature = Program.RestartTemperature;

	    config.NumberDataPoints = Program.NumberDataPoints;
	    config.ProcessingOption = Program.ProcessingOption;
	    config.MaxNcentPerNode = Program.maxNcentperNode;
	    config.ThreadCount = DAVectorUtility.ThreadCount;
	    config.NodeCount = DAVectorUtility.NodeCount;
	    config.TooSmallToSplit = Program.ToosmalltoSplit;
	    config.WaitIterations = Program.Waititerations;
	    config.InitialCoolingFactor = Program.InitialCoolingFactor;
	    config.FineCoolingFactor = Program.FineCoolingFactor;
	    config.MalphaMaxChange = Program.Malpha_MaxChange;
	    config.IterationAtEnd = Program.Iterationatend;
	    config.ConvergenceLoopLimit = Program.ConvergenceLoopLimit;
	    config.FreezingLimit = Program.FreezingLimit;
	    config.ConvergeIntermediateClusters = Program.ConvergeIntermediateClusters;
	    config.DebugPrintOption = DAVectorUtility.DebugPrintOption;
	    config.ConsoleDebugOutput = DAVectorUtility.ConsoleDebugOutput;
	} 

	public static void WriteTimingFile(String fileName, long duration, double HPDuration, int ThreadCount, int MPIperNodeCount, int NodeCount, int PointCount_Process, int PointCount_Global, int maxNcent, double PTsplitting, double Ratio1, double MPIReduce, double Ratio2, double MPISRBasic, double Ratio3, double MPISREigen, double Ratio4, double MPIBdcast, double Ratio5, String DataFileName, java.util.Date CurrentTime, String ProcessorName)
	{
        String header =
                "Duration(ms)\tHPDuration(us)" +
                        "\tThread#\tMPIperNode\tNode#\tPattern\tParallelism\tPoint#/Process\tPoint#/Global\tmaxNcent" +
                        "\tPTsplitting\tRatio1\tMPIReduce\tRatio2\tMPISRBasic\tRatio3\tMPISREigen\tRatio4\tMPIBdcast" +
                        "\tRatio5\tDataFileName\tCurrentTime\tProcessorName";

        String format = "%1$s\t%2$s\t%3$s\t%4$s\t%5$s\t%6$s\t%7$s\t%8$s\t%9$s\t%10$s\t%11$s\t%12$s\t%13$s\t%14$s\t%15$s\t%16$s\t%17$s\t%18$s\t%19$s\t%20$s\t%21$s\t%22$s\t%23$s";
//        String format = "%1$\t%2$\t%3$\t%4$\t%5$\t%6$\t%7$\t%8$\t%9$\t%10$\t%11$\t%12$\t%13$\t%14$\t%15$\t%16$\t%17$\t%18$\t%19$\t%20$\t%21$\t%22$\t%23$";
        String line =
                String.format(format,
                              Math.round(duration),
                              Math.round(HPDuration),
                              ThreadCount,
                              MPIperNodeCount,
                              NodeCount,
                              String.format("%1$sx%2$sx%3$s", ThreadCount, MPIperNodeCount, NodeCount),
                              (ThreadCount * MPIperNodeCount * NodeCount),
                              PointCount_Process,
                              PointCount_Global,
                              maxNcent,
                              Math.round(PTsplitting),
                              new BigDecimal(Ratio1).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPIReduce),
                              new BigDecimal(Ratio2).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPISRBasic),
                              new BigDecimal(Ratio3).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPISREigen),
                              new BigDecimal(Ratio4).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPIBdcast),
                              new BigDecimal(Ratio5).setScale(4, RoundingMode.UP).doubleValue(),
                              DataFileName, CurrentTime,
                              ProcessorName); // Name of input data file -  MPIBdcast vs HPDuration -  MPI broadcast -  MPISREigen vs HPDuration -  MPI SR Eigen -  MPISRBasic vs HPDuration -  MPI SR Basic -  MPIReduce vs HPDuration -  MPI reduce -  PTsplitting vs HPDuration -  Partial times splitting -  cluster# -  Global points -  Local points -  Pattern -  Node# aquired -  Process# per node -  Thread# per process -  High performance timer -  Total time

        Path filePath = Paths.get(fileName);
        try {
            if (!Files.exists(filePath))
            {
                Files.write(filePath, header.getBytes(), StandardOpenOption.CREATE);
            }
            Files.write(filePath, (System.lineSeparator() + line).getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e) {
            System.err.format("Failed writing timing file due to I/O exception: %s%n", e);
        }
    }

    public static void WriteTimingFileNew(String filename){
        Path file = Paths.get(filename);
        StandardOpenOption option = (Files.exists(file)) ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;

        try(PrintWriter printWriter = new PrintWriter(Files.newBufferedWriter(file,option))){
            printWriter.println("General Timing K-Means " + DAVectorUtility.formatElapsedMillis(GeneralTiming.getTotalTime(GeneralTiming.TimingTask.KMEANS)));
            printWriter.println("General Timing LCMS " + DAVectorUtility.formatElapsedMillis(GeneralTiming.getTotalTime(GeneralTiming.TimingTask.LCMS)));
            printWriter.println("General Timing DA " + DAVectorUtility.formatElapsedMillis(GeneralTiming.getTotalTime(GeneralTiming.TimingTask.DA)));

            printWriter.println();

            printWriter.println("ReadData Timing " + DAVectorUtility.formatElapsedMillis(ReadDataTiming.getTotalTime(ReadDataTiming.TimingTask.TOTAL_READ)));
            printWriter.println("Setup Timing " + DAVectorUtility.formatElapsedMillis(SetupTiming.getTotalTime(SetupTiming.TimingTask.SETUP)));

            printWriter.println();

            printWriter.println("Section Timing SC1 " + DAVectorUtility.formatElapsedMillis(SectionTiming.getTotalTime(SectionTiming.TimingTask.SC1)));
            printWriter.println("Section Timing SC2 " + DAVectorUtility.formatElapsedMillis(SectionTiming.getTotalTime(SectionTiming.TimingTask.SC2)));
            printWriter.println("Section Timing SC3 " + DAVectorUtility.formatElapsedMillis(SectionTiming.getTotalTime(SectionTiming.TimingTask.SC3)));

			printWriter.println();
			printWriter.println("General Method Timing ClusterStatus.SetClusterStatistics " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.SET_CLUSTER_STATISTICS)));
			printWriter.println("General Method Timing ClusterStatus.SetPointStatistics " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.SET_POINT_STATISTICS)));
			printWriter.println("General Method Timing ClusterStatus.SetNearbyClusters " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.SET_NEARBY_CLUSTERS)));
			printWriter.println("General Method Timing ClusterStatus.OutputStatus " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.OUTPUT_STATUS)));
			printWriter.println("General Method Timing ClusterStatus.ExperimentAnalysis " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.EXPERIMENT_ANALYSIS)));
			printWriter.println("General Method Timing OurClusters.setup " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.SETUP)));
			printWriter.println("General Method Timing LCMSAnalyze.ClusterComparison " + DAVectorUtility.formatElapsedMillis(GeneralMethodTiming.getTotalTime(GeneralMethodTiming.TimingTask.CLUSTER_COMPARISON)));

			printWriter.println();
			printWriter.println("Total Run time " + DAVectorUtility.formatElapsedMillis(GeneralTiming.getTotalTime(GeneralTiming.TimingTask.TOTAL)));


			printWriter.flush();
            printWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

	public static void SetupDA()
	{
		if (!Program.DoBasicDA)
		{
			return;
		}

		//  Avoid Distributed Execution
		Program.TemperatureLimitforDistribution = -1.0;
		Program.ClusterLimitforDistribution = -1;

		Program.RW3DData = 3; // Read and Write Plotviz

		Program.CalculateEigenvaluesfromMatrix = false;
		Program.UseSponge = false;
		Program.ContinuousClustering = true;
		Program.SigmaMethod = 0;

		Program.maxNcentperNode = 140;
		Program.maxNcentperNode = 10;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 8;
		Program.InitialNcent = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
		Program.maxNcentCreated = 500 * Program.maxNcentperNode;

		Program.NumberDataPoints = 200000;
		Program.StartPointPositiononInputLine = 0;
		Program.ParameterVectorDimension = 74;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;

		Program.MagicTemperatures[0] = -1.0;
		Program.MagicTemperatures[1] = -1.0;
		Program.MagicTemperatures[2] = -1.0;
		Program.MagicTemperatures[3] = -1.0;
		Program.MagicTemperatures[4] = -1.0;

		Program.Iterationatend = 800;
		Program.ConvergenceLoopLimit = 400;
		Program.Malpha_MaxChange = 0.005;
		Program.Malpha_MaxChange1 = 0.0005;
		Program.YChangeSquared = 0.00000001;
		Program.ExpArgumentCut1 = 20.0;
		Program.ExpArgumentCut2 = 40.0;
		Program.ExpArgumentCut2 = 20.0;
		Program.Waititerations = 1;
		Program.Waititerations_Converge = 4;
		Program.MinimumCountforCluster_Points = -1;
		Program.InitialCoolingFactor1 = 0.95;
		Program.FineCoolingFactor1 = 0.995;
		Program.InitialCoolingFactor2 = 0.95;
		Program.FineCoolingFactor2 = 0.995;
		Program.MinimumScaledWidthsquaredtosplit = -1.0;
		Program.ScaledSquaredDistanceatClosenessTest = 0.01;
		Program.TemperatureforClosenessTest = 0.01;
		Program.MinimumCountforCluster_C_k = 8.0; // Remove clusters with fewer points than this (based on C_k_) (This only used for converged clusters)
		Program.MinimumCountforCluster_Points = 10; // Remove clusters with fewer points than this (based on assigned points)

		Program.MaxNumberSplitClusters = 16;
		Program.ToosmalltoSplit = 200.0;
		Program.eigenvaluechange = 0.001; // Limit 1 on Eigenvalue Changes
		Program.eigenvectorchange = 0.001; // Limit 2 on Eigenvalue Changes

		Program.UseTriangleInequality_DA = 1;
		Program.OldCenterOption_DA = -1;
		Program.maxNcentTOTALforParallelism_DA = 20; // Use Center parallelism when this limit hit
		Program.TriangleInequality_Delta1_old_DA = 0.1; // Test for center change and old lower bounds (normalized by radius)
		Program.TriangleInequality_Delta1_current_DA = 0.1; // Test for Center change and current lower bounds (normalized by radius)
		Program.MaxClusterLBsperPoint_DA = 140; // Maximum number of Lower Bound values
		Program.MaxCentersperCenter_DA = 200; // Maximum number of Centers in Center Difference Array

		Program.Tminimum = -20000.0;

		Program.PrintInterval = 5;
		Program.ClusterPrintNumber = 10;
		DAVectorUtility.DebugPrintOption = 2;

		/* 85399 54D
		Program.maxNcentperNode = 100;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 4;
		Program.InitialNcent = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
		Program.maxNcentCreated = 500 * Program.maxNcentperNode;

		Program.NumberDataPoints = 85399;
		Program.StartPointPositiononInputLine = 0;
		Program.ParameterVectorDimension = 54;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;
		Program.MaxNumberSplitClusters = 4;
		85399 54D */


		/*  Cmeans FLAME Test
		Program.MaxNumberSplitClusters = 1;
		Program.UseTriangleInequality_DA = 1;

		Program.maxNcentperNode = 5;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 5;
		Program.InitialNcent = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
		Program.maxNcentCreated = 500 * Program.maxNcentperNode;

		Program.NumberDataPoints = 22014;
		Program.StartPointPositiononInputLine = 1;
		Program.ParameterVectorDimension = 4;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;

		Program.UseTriangleInequality_DA = 0;
		*/
		/* Lung
		Program.MaxNumberSplitClusters = 1;
		Program.UseTriangleInequality_DA = 1;

		Program.maxNcentperNode = 20;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 5;
		Program.InitialNcent = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
		Program.maxNcentCreated = 500 * Program.maxNcentperNode;

		Program.NumberDataPoints = 20054;
		Program.StartPointPositiononInputLine = 0;
		Program.ParameterVectorDimension = 4;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;

		Program.UseTriangleInequality_DA = 1;
		Lung */ 

		//  1000 point Dating Run
		Program.MaxNumberSplitClusters = 4;
		Program.Waititerations = 4;
		Program.InitialCoolingFactor1 = 0.995;
		Program.FineCoolingFactor1 = 0.9995;
		Program.InitialCoolingFactor2 = 0.995;
		Program.FineCoolingFactor2 = 0.9995;
		Program.ToosmalltoSplit = 15.0;

		Program.maxNcentperNode = 30;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 30;
		Program.InitialNcent = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
		Program.maxNcentCreated = 500 * Program.maxNcentperNode;

		Program.NumberDataPoints = 1000;
		Program.StartPointPositiononInputLine = 1;
		Program.ParameterVectorDimension = 3;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;

		Program.UseTriangleInequality_DA = 0;

	} // End SetupDA

	public static void SetupKmeans()
	{
		if (!Program.DoKmeans)
		{
			return;
		}

		Program.RW3DData = 3; // Read and Write Plotviz

		//  Avoid Distributed Execution
		Program.TemperatureLimitforDistribution = -1.0;
		Program.ClusterLimitforDistribution = -1;
		Program.SelectedInputLabel = 0;

		Program.CalculateEigenvaluesfromMatrix = false;
		Program.UseSponge = false;
		Program.ContinuousClustering = false;
		Program.SigmaMethod = 0;
		Program.targetNcentperPoint = 124;

		Program.maxNcentperNode = 140;
		Program.maxNcentperNode = 10;
		Program.InitialNcent = Program.maxNcentperNode;
		Program.targetNcentperPoint = Program.maxNcentperNode;
		Program.maxNcentperPoint = Program.maxNcentperNode;
		Program.targetMinimumNcentperPoint = 8;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentCreated = 200;
		Program.maxNcentTOTALforParallelism_Kmeans = 20;

		Program.UseTriangleInequality_Kmeans = 0; // 0 is Pure K means
		Program.UseTriangleInequality_Kmeans = 1; // 0 is Pure K means

		Program.MaxClusterLBsperPoint_Kmeans = Program.maxNcentperNode;
		Program.MaxCentersperCenter_Kmeans = Program.maxNcentperNode;

		Program.NumberDataPoints = 200000;
		Program.ClusterIndexonInputLine = 75;
		Program.ClusterIndexonInputLine = -1;
		Program.ParameterVectorDimension = 74;
		Program.StartPointPositiononInputLine = 0;
		Program.FirstClusterValue = 0;

		Program.KmeansCenterChangeStop = 0.00001;
		Program.KmeansIterationLimit = 1000;
		Program.TriangleInequality_Delta1_old_Kmeans = 0.2;
		Program.TriangleInequality_Delta1_current_Kmeans = 0.2;

        // Dating-1000-1
        Program.ContinuousClustering = true;
        Program.InitialNcent = 1;
        Program.maxNcentperNode = 3;
        Program.targetNcentperPoint = Program.maxNcentperNode;
        Program.maxNcentperPoint = Program.maxNcentperNode;
        Program.UseTriangleInequality_Kmeans = 0; // 0 is Pure K means
        Program.NumberDataPoints = 1000;
        Program.ParameterVectorDimension = 3;

		// Small Crandall Data
		/*Program.RW3DData = -1;
		Program.ParameterVectorDimension = 2048;
		Program.StartPointPositiononInputLine = 4;
		Program.ClusterCountOutput = -1; // No output

		Program.SelectedInputLabel = 0;
		Program.InputFileType = 0;
		Program.Replicate = 1;
		Program.FirstClusterValue = -1;

		Program.NumberDataPoints = 76800;

		Program.maxNcentperNode = 3200;
		Program.InitialNcent = Program.maxNcentperNode;
		Program.targetNcentperPoint = 1;
		Program.maxNcentperPoint = 1;
		Program.targetMinimumNcentperPoint = 1;
		Program.maxNcentTOTAL = Program.maxNcentperNode;
		Program.maxNcentCreated = Program.maxNcentperNode;
		Program.maxNcentTOTALforParallelism_Kmeans = 20;

		Program.UseTriangleInequality_Kmeans = 0; // 0 is Pure K means
		Program.ClusterIndexonInputLine = -3;
		Program.OldCenterOption_Kmeans = 100;
		Program.DoBackwardFacingTests_Kmeans = false;

		int CutBounds1 = 8;
		int CutBounds2 = 8;
		Program.MaxClusterLBsperPoint_Kmeans = Program.maxNcentperNode / CutBounds1;
		Program.MaxCentersperCenter_Kmeans = Program.maxNcentperNode / CutBounds2;*/

		/*  Cmeans/Flame
		Program.NumberDataPoints = 22014;
		Program.StartPointPositiononInputLine = 1;
		Program.ParameterVectorDimension = 4;
		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = 5;
		Program.InputFileType = 0;
		Program.Replicate = 1;
		Program.FirstClusterValue = 1;


		Program.InitialNcent = 5;
		Program.targetNcentperPoint = Program.InitialNcent;
		Program.maxNcentperNode = Program.InitialNcent;
		Program.maxNcentTOTAL = Program.InitialNcent;
		Program.maxNcentTOTALforParallelism_Kmeans = Program.InitialNcent + 1;
		Program.maxNcentCreated = Program.InitialNcent;

		Program.UseTriangleInequality_Kmeans = 0;
		*/

		/* Lung Data

		Program.InitialNcent = 100;
		Program.targetNcentperPoint = Program.InitialNcent;
		Program.maxNcentperNode = Program.InitialNcent;
		Program.maxNcentTOTAL = Program.InitialNcent;
		Program.maxNcentTOTALforParallelism_Kmeans = Program.InitialNcent + 1;
		Program.maxNcentCreated = Program.InitialNcent;

		Program.NumberDataPoints = 85399;
		Program.StartPointPositiononInputLine = 0;
		Program.ParameterVectorDimension = 54;

		Program.SelectedInputLabel = 0;
		Program.ClusterIndexonInputLine = -1;
		Program.InputFileType = 0;
		Program.Replicate = 1;
		Program.MaxNumberSplitClusters = 4;
		Program.RW3DData = 3;   // Read and Write Plotviz
		*/

		/* Read Plot file
		Program.ParameterVectorDimension = 3;
		Program.StartPointPositiononInputLine = 1;
		*/

	} // End SetupKmeans()

	public static void SetupLCMS()
	{ // Set up LC MS 2D analysis

		if (!Program.DoLCMS)
		{
			return;
		}

		int LCMSmode = config.LCMSmode;

		//  Generic Set up
		LCMSAnalyze.InitializeLCMS();

		//  Do a Harvard Format Full Analysis
		if (LCMSmode == 1)
		{
			LCMSAnalyze.FreshFullAnalysis();
			return;
		}

		//  Analyse Sponge Data from a Previous Run
		if (LCMSmode == 2)
		{
			LCMSAnalyze.SpongePointAnalysis();
			return;
		}

		//  Join and Equilibriate 2 Files together (original Full Analysis and Analaysis of its Sponge)
		//  Or just equilibriate a single file
		if (LCMSmode == 3)
		{
			LCMSAnalyze.JoinEquilibriate();
			return;
		}

		//  Change Sponge Factor
		if (LCMSmode == 4)
		{
			LCMSAnalyze.ChangeSpongeFactor();
			return;
		}

		//  Just Analyze a file
		if (LCMSmode == 5)
		{
			LCMSAnalyze.JustAnalyze();
        }


	}
}