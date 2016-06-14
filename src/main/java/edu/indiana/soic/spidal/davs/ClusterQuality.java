package edu.indiana.soic.spidal.davs;

import com.google.common.io.Files;
import edu.indiana.soic.spidal.general.Box;
import mpi.MPIException;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

//  Data structure to record quality of clusters
//  SpongePointsinCut Number of Sponge Points within Distance Cut
//  NumberClustersinCut Number of Clusters within Distance Cut
//  NearbyClusterIndices Labels of Nearby Clusters sorted by Distance
//  NearbyClusterDistances Distance of Nearby Clusters
public class ClusterQuality
{
	private int NumberClusters = 0;
	private int SpongeOption = -1;
	private int NumberNearbyClusters = 0;
	private double NearbyCut = 0.0;
	private double TotalConfusedClusterPoints = 0.0;
	private double TotalConfusedSpongePoints1 = 0.0;
	private double TotalConfusedSpongePoints2 = 0.0;
	private double AverageClustersinCut = 0.0;
	private int TotalIllegalPoints = 0;
	private int HistogramOccupationMax = 202;
	private int HistogramDistancesMax = 200;

	private int[] SpongePointsinCut;
	private int[] NonSpongePointsinCut;
	private int[] NumberClustersinCut;
	private double[] ClusterMaxWidths;
	private int[] ClusterIllegalPoints;
	private int[][] NearbyClusterIndices;
	private double[][] NearbyClusterDistances;
	private int[] HistogramOccupationValues;
	private int[] HistogramDistancesValues;

	//  InputNumberClusters -- Number of Clusters
	//  InputSpongeOption -- Index of Sponge or if < 0 No sponge
	//  InputNumberNearbyClusters -- List this number of nearby clusters
	//  InputNearbyCut -- Count Sponge Points and Clusters within this cut times sigma
	public ClusterQuality(int InputNumberClusters, int InputSpongeOption, int InputNumberNearbyClusters, double InputNearbyCut)
	{
		NumberClusters = InputNumberClusters;
		SpongeOption = InputSpongeOption;
		NumberNearbyClusters = Math.min(InputNumberNearbyClusters, NumberClusters - 1);
		NearbyCut = InputNearbyCut * InputNearbyCut;

		HistogramOccupationValues = new int[HistogramOccupationMax];
		HistogramDistancesValues = new int[HistogramDistancesMax];

		SpongePointsinCut = new int[NumberClusters];
		NonSpongePointsinCut = new int[NumberClusters];
		NumberClustersinCut = new int[NumberClusters];
		NearbyClusterIndices = new int[NumberClusters][];
		ClusterIllegalPoints = new int[NumberClusters];
		ClusterMaxWidths = new double[NumberClusters];
		NearbyClusterDistances = new double[NumberClusters][];

		if (InputNumberNearbyClusters > 0)
		{
			for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
			{
				NearbyClusterIndices[ClusterIndex] = new int[NumberNearbyClusters];
				NearbyClusterDistances[ClusterIndex] = new double[NumberNearbyClusters];
			}
			if (SpongeOption >= 0)
			{
				for (int spongespecial = 0; spongespecial < NumberNearbyClusters; spongespecial++)
				{
					NearbyClusterIndices[SpongeOption][spongespecial] = 0;
					NearbyClusterDistances[SpongeOption][spongespecial] = 0.0;
				}
				NumberClustersinCut[SpongeOption] = 0;
			}
		}
	} // End Initializer ClusterQuality(int InputNumberClusters, int InputSpongeOption, int InputNumberNearbyClusters, double InputNearbyCut)

    public final void SetClusterStatistics() throws MPIException {
        GlobalReductions.FindDoubleArraySum PointsperCluster =
                new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, HistogramOccupationMax);
        GlobalReductions.FindDoubleArraySum DistancesfromCluster =
                new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, HistogramDistancesMax);
        GlobalReductions.FindIntSum TotalClustersinCut = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);
        double DistanceCut = (double) HistogramDistancesMax;

        Range[] ParallelClusterRange = RangePartitioner.Partition(NumberClusters, DAVectorUtility.ThreadCount);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                PointsperCluster.startthread(threadIndex);
                DistancesfromCluster.startthread(threadIndex);
                if (DAVectorUtility.MPI_Rank == 0) {
                    int beginindex = ParallelClusterRange[threadIndex].getStartIndex();
                    int indexlength = ParallelClusterRange[threadIndex].getLength();
                    double[] ClusterSigma = new double[Program.ParameterVectorDimension];
                    for (int ClusterIndex1 = beginindex; ClusterIndex1 < beginindex + indexlength; ClusterIndex1++) {
                        int count = 0;
                        if (ClusterIndex1 == SpongeOption) {
                            continue;
                        }
                        int histcount = ClusteringSolution.TotalClusterSummary.OccupationCount[ClusterIndex1];
                        if (histcount > (HistogramOccupationMax - 1)) {
                            histcount = HistogramOccupationMax - 1;
                        }
                        PointsperCluster.addapoint(threadIndex, histcount);
                        for (int ClusterIndex2 = 0; ClusterIndex2 < NumberClusters; ClusterIndex2++) {
                            if (ClusterIndex2 == SpongeOption) {
                                continue;
                            }
                            if (ClusterIndex1 == ClusterIndex2) {
                                continue;
                            }
                            Box<double[]> tempRef_ClusterSigma =
                                    new Box<>(ClusterSigma);
                            Program
                                    .CalculateSigma(
                                            ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                                            tempRef_ClusterSigma);
                            ClusterSigma = tempRef_ClusterSigma.content;
                            double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(
                                    ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                                    ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex2], ClusterSigma);
                            if (tmp <= DistanceCut) {
                                int itmp = (int) Math.floor(tmp + 0.5);
                                itmp = Math.min(itmp, HistogramDistancesMax - 1);
                                DistancesfromCluster.addapoint(threadIndex, itmp);
                            }
                            if (tmp < NearbyCut) {
                                ++count;
                            }
                        }
                        NumberClustersinCut[ClusterIndex1] = count;
                        TotalClustersinCut.addapoint(threadIndex, count);
                    }
                }

            });
        });
        TotalClustersinCut.sumoverthreadsandmpi();

        double ClusterSum = NumberClusters;
        if (SpongeOption >= 0) {
            ClusterSum -= 1.0;
        }
        AverageClustersinCut = (double) TotalClustersinCut.TotalInt / (2.0 * ClusterSum);

        PointsperCluster.sumoverthreadsandmpi();
        DistancesfromCluster.sumoverthreadsandmpi();
        for (int histloop = 0; histloop < HistogramOccupationMax; histloop++) {
            HistogramOccupationValues[histloop] = (int) (PointsperCluster.TotalSum[histloop] + 0.001);
        }
        for (int histloop = 0; histloop < HistogramDistancesMax; histloop++) {
            HistogramDistancesValues[histloop] = (int) (DistancesfromCluster.TotalSum[histloop] + 0.001);
        }

    } // End SetClusterStatistics()

	public final void SetPointStatistics() throws MPIException {
		TotalIllegalPoints = 0;

		GlobalReductions.FindDoubleArraySum AccumulateSpongePoints = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, NumberClusters);
		GlobalReductions.FindDoubleArraySum AccumulateNonSpongePoints = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, NumberClusters);
		GlobalReductions.FindDoubleSum ConfusedSpongePoints1 = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
		GlobalReductions.FindDoubleSum ConfusedSpongePoints2 = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
		GlobalReductions.FindDoubleSum ConfusedClusterPoints = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
		GlobalReductions.FindDoubleMax[] FindMaxExtension = new GlobalReductions.FindDoubleMax[NumberClusters];
		GlobalReductions.FindIntSum[] FindIllegalPoints = new GlobalReductions.FindIntSum[NumberClusters];
		for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
		{
			FindMaxExtension[ClusterIndex] = new GlobalReductions.FindDoubleMax(DAVectorUtility.ThreadCount);
			FindIllegalPoints[ClusterIndex] = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);
		}

		//  Parallel Section over points
        //  Remove Sponge Confused Points
        //  Not in Sponge
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                double[] ClusterSigma = new double[Program.ParameterVectorDimension];
                AccumulateSpongePoints.startthread(threadIndex);
                AccumulateNonSpongePoints.startthread(threadIndex);
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                    int AssignedCluster = Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process];
                    boolean isSponge = false;
                    if (AssignedCluster == SpongeOption) {
                        isSponge = true;
                    }
                    if (isSponge) {
                        int nearbyrealcluster = -1;
                        double NearesetDistce = -1.0;
                        for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
                            if (ClusterIndex == SpongeOption) {
                                continue;
                            }
                            if (ClusterIndex == AssignedCluster) {
                                continue;
                            }
                            Box<double[]> tempRef_ClusterSigma = new Box<>(ClusterSigma);
                            Program.CalculateSigma(
                                    ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex],
                                    tempRef_ClusterSigma);
                            ClusterSigma = tempRef_ClusterSigma.content;
                            double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(
                                    Program.PointPosition[index],
                                    ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ClusterSigma);
                            if (tmp < NearbyCut) {
                                if ((nearbyrealcluster < 0) || (tmp < NearesetDistce)) {
                                    nearbyrealcluster = ClusterIndex;
                                    NearesetDistce = tmp;
                                }
                            }
                        }
                        if (nearbyrealcluster >= 0) {
                            Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process] = nearbyrealcluster;
                            AssignedCluster = nearbyrealcluster;
                            isSponge = false;
                            ConfusedSpongePoints1.addapoint(threadIndex, 1.0);
                        }
                    }
                    if (!isSponge) {
                        if ((AssignedCluster < 0) || (AssignedCluster >= NumberClusters)) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    "Point " + (index + DAVectorUtility.PointStart_Process) + " Assigned " + AssignedCluster + " Max " + NumberClusters);

                        }
                        Box<double[]> tempRef_ClusterSigma2 = new Box<>(ClusterSigma);
                        Program.CalculateSigma(
                                ClusteringSolution.TotalClusterSummary.CenterPosition[AssignedCluster],
                                tempRef_ClusterSigma2);
                        ClusterSigma = tempRef_ClusterSigma2.content;
                        double distce = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(
                                Program.PointPosition[index],
                                ClusteringSolution.TotalClusterSummary.CenterPosition[AssignedCluster], ClusterSigma);
                        FindMaxExtension[AssignedCluster].addapoint(threadIndex, distce);
                        if (distce > NearbyCut) {
                            FindIllegalPoints[AssignedCluster].addapoint(threadIndex, 1);
                            if (SpongeOption >= 0) {
                                AssignedCluster = SpongeOption;
                                Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process] = AssignedCluster;
                            }
                        }
                    }
                    boolean SpongeConfused = false;
                    boolean ClusterConfused = false;
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
                        if (ClusterIndex == SpongeOption) {
                            continue;
                        }
                        if (ClusterIndex == AssignedCluster) {
                            continue;
                        }
                        Box<double[]> tempRef_ClusterSigma3 = new Box<>(ClusterSigma);
                        Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex],
                                tempRef_ClusterSigma3);
                        ClusterSigma = tempRef_ClusterSigma3.content;
                        double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(
                                Program.PointPosition[index],
                                ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ClusterSigma);
                        if (tmp < NearbyCut) {
                            if (isSponge) {
                                AccumulateSpongePoints.addapoint(threadIndex, ClusterIndex);
                                SpongeConfused = true;
                            } else {
                                AccumulateNonSpongePoints.addapoint(threadIndex, ClusterIndex);
                                ClusterConfused = true;
                            }
                        }
                    }
                    if (SpongeConfused) {
                        ConfusedSpongePoints2.addapoint(threadIndex, 1.0);
                    }
                    if (ClusterConfused) {
                        ConfusedClusterPoints.addapoint(threadIndex, 1.0);
                    }
                }

            });
        });

        AccumulateSpongePoints.sumoverthreadsandmpi();
		AccumulateNonSpongePoints.sumoverthreadsandmpi();
		ConfusedClusterPoints.sumoverthreadsandmpi();
		ConfusedSpongePoints1.sumoverthreadsandmpi();
		ConfusedSpongePoints2.sumoverthreadsandmpi();

		for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
		{
			FindIllegalPoints[ClusterIndex].sumoverthreadsandmpi();
			FindMaxExtension[ClusterIndex].sumoverthreadsandmpi();
			SpongePointsinCut[ClusterIndex] = (int)Math.floor(AccumulateSpongePoints.TotalSum[ClusterIndex] + 0.5);
			NonSpongePointsinCut[ClusterIndex] = (int)Math.floor(AccumulateNonSpongePoints.TotalSum[ClusterIndex] + 0.5);
			ClusterIllegalPoints[ClusterIndex] = FindIllegalPoints[ClusterIndex].TotalInt;
			ClusterMaxWidths[ClusterIndex] = FindMaxExtension[ClusterIndex].TotalMax;
			TotalIllegalPoints += ClusterIllegalPoints[ClusterIndex];
		}
		TotalConfusedClusterPoints = ConfusedClusterPoints.Total;
		TotalConfusedSpongePoints1 = ConfusedSpongePoints1.Total;
		TotalConfusedSpongePoints2 = ConfusedSpongePoints2.Total;

	} // End SetPointStatistics()

	public final void SetNearbyClusters()
	{
		Range[] ParallelClusterRange = RangePartitioner.Partition(NumberClusters, DAVectorUtility.ThreadCount);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) -> {
                int beginindex = ParallelClusterRange[threadIndex].getStartIndex();
                int indexlength = ParallelClusterRange[threadIndex].getLength();
                double[] ClusterSigma = new double[Program.ParameterVectorDimension];
                int[] localindices = new int[NumberNearbyClusters];
                double[] localdistces = new double[NumberNearbyClusters];
                for (int ClusterIndex1 = beginindex; ClusterIndex1 < beginindex + indexlength; ClusterIndex1++) {
                    if (ClusterIndex1 == SpongeOption) {
                        continue;
                    }
                    for (int nearbyloop = 0; nearbyloop < NumberNearbyClusters; nearbyloop++) {
                        localindices[nearbyloop] = -1;
                        localdistces[nearbyloop] = 0.0;
                    }
                    int worstposition = -1;
                    for (int ClusterIndex2 = 0; ClusterIndex2 < NumberClusters; ClusterIndex2++) {
                        if (ClusterIndex2 == SpongeOption) {
                            continue;
                        }
                        if (ClusterIndex1 == ClusterIndex2) {
                            continue;
                        }
                        Box<double[]> tempRef_ClusterSigma = new Box<>(ClusterSigma);
                        Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                                tempRef_ClusterSigma);
                        ClusterSigma = tempRef_ClusterSigma.content;
                        double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(
                                ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                                ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex2], ClusterSigma);
                        Box<Integer> tempRef_worstposition = new Box<>(worstposition);
                        GlobalReductions.FindMinimumSet(tmp, ClusterIndex2, tempRef_worstposition, localdistces,
                                localindices, NumberNearbyClusters);
                        worstposition = tempRef_worstposition.content;
                    }
                    for (int countpos = 0; countpos < NumberNearbyClusters; countpos++) {
                        int nearbyloop = -1;
                        for (int nearbyloop1 = 0; nearbyloop1 < NumberNearbyClusters; nearbyloop1++) {
                            if (localindices[nearbyloop1] < 0) {
                                continue;
                            }
                            if ((nearbyloop < 0) || (localdistces[nearbyloop1] <= localdistces[nearbyloop])) {
                                nearbyloop = nearbyloop1;
                            }
                        }
                        NearbyClusterIndices[ClusterIndex1][countpos] = localindices[nearbyloop];
                        NearbyClusterDistances[ClusterIndex1][countpos] = localdistces[nearbyloop];
                        localindices[nearbyloop] = -1;
                        localdistces[nearbyloop] = 0.0;
                    }
                }

            });
        });

    } // End SetNearbyClusters()

	public final void OutputStatus()
	{
		if (DAVectorUtility.MPI_Rank != 0)
		{
			return;
		}
		String nextline = "\n---------------- Cluster Status with Cut " + String.format("%1$5.4f", NearbyCut) + " Nearby Cluster Limit " + NumberNearbyClusters + " WRITTEN to status file " + " Sponge " + SpongeOption + "\nIllegal Points outside cut of center MOVED to Sponge " + TotalIllegalPoints + " Confused Sponge Points MOVED to Real Cluster " + String.format("%1$5.4E", TotalConfusedSpongePoints1) + " After Correction " + String.format("%1$5.4E", TotalConfusedSpongePoints2) + " Confused Cluster Points near >1 cluster " + String.format("%1$5.4E", TotalConfusedClusterPoints) + " Confused Cluster Fraction (Clusters within Cut)" + String.format("%1$3.2f", AverageClustersinCut);
		DAVectorUtility.SALSAPrint(0, nextline);

		nextline = "Cluster Occupation Count Histogram\n";
		for (int histloop = 0; histloop < HistogramOccupationMax; histloop++)
		{
			nextline += HistogramOccupationValues[histloop] + ", ";
		}
		nextline += "\n\nInter Cluster Distance Histogram\n";
		for (int histloop = 0; histloop < HistogramDistancesMax; histloop++)
		{
			nextline += HistogramDistancesValues[histloop] + ", ";
		}
		DAVectorUtility.SALSAPrint(0, nextline);

		//  Output Global Counts versus Temperature
		ClusterQuality.CaculateTemperatureClusterCountPlot();

		String directory = (new File(Program.config.ClusterFile)).getParent();
		String file = Files.getNameWithoutExtension(Program.config.ClusterFile) + "Status1" + "-M" + Program.maxNcentperNode + "-C" + ParallelClustering.runningSolution.Ncent_Global + "." + Files.getFileExtension(
                Program.config.ClusterFile);
		String StatusFileName = Paths.get(directory, file).toString();
		try (PrintWriter writer = new PrintWriter(java.nio.file.Files.newBufferedWriter(Paths.get(StatusFileName),
                Charset.defaultCharset(), StandardOpenOption.CREATE),true))
		{
            nextline = "\nLabel\tY-0\tY-1\tCount\tNearby Clusters\tNearby Sponge Points\tIllegal Points\tMax Extension\tNearby Cluster Points";
            writer.println(nextline);
            for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
            {
                if (ClusterIndex == SpongeOption)
                {
                    continue;
                }
                nextline = DAVectorUtility.PrintFixedInteger(ClusterIndex, 5);
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    double tmp = ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex][VectorIndex];
                    if ((Program.SigmaMethod >= 2) && (VectorIndex == 0))
                    {
                        tmp = Math.log(tmp);
                    }
                    tmp = tmp / Program.SigmaVectorParameters_i_[VectorIndex];
                    nextline += "\t" + DAVectorUtility.PadString(String.format("%1$5.4f", tmp), 10);
                }
                nextline += "\t" + DAVectorUtility.PrintFixedInteger(ClusteringSolution.TotalClusterSummary.OccupationCount[ClusterIndex], 4) + "\t" + DAVectorUtility.PrintFixedInteger(NumberClustersinCut[ClusterIndex], 4) + "\t" + DAVectorUtility.PrintFixedInteger(SpongePointsinCut[ClusterIndex], 4) + "\t" + DAVectorUtility.PrintFixedInteger(ClusterIllegalPoints[ClusterIndex], 4) + "\t" + DAVectorUtility.PadString(String.format("%1$4.3f", ClusterMaxWidths[ClusterIndex]), 8) + "\t" + DAVectorUtility.PrintFixedInteger(NonSpongePointsinCut[ClusterIndex], 4);

                String inmessage = "";
                for (int nearbyloop = 0; nearbyloop < NumberNearbyClusters; nearbyloop++)
                {
                    int ClusterNeighbor = NearbyClusterIndices[ClusterIndex][nearbyloop];
                    if (ClusterIndex != ClusterNeighbor)
                    {
                        inmessage += "\t" + ClusterNeighbor + "(" + String.format("%1$3.2f", NearbyClusterDistances[ClusterIndex][nearbyloop]) + ") ";
                    }
                }
                if (inmessage.length() > 0)
                {
                    nextline += inmessage;
                }
                writer.println(nextline);
            }
            writer.close();
		} catch (IOException e) {
            System.err.format("Failed writing status due to I/O exception: %s%n", e);
        }
	} // End OutputStatus()

    public final void ExperimentAnalysis(){
        int NumberofExperiments = 37;

        if(NumberofExperiments <= 1) return;

        int ExperimentPoints[]  = new int[NumberofExperiments];
        int SpongeExptPoints[]  = new int[NumberofExperiments];
        double xshift[] = new double[NumberofExperiments];
        double yshift[] = new double[NumberofExperiments];
        double xwidth[] = new double[NumberofExperiments];
        double ywidth[] = new double[NumberofExperiments];

        int SpongeClusterNumber=ClusteringSolution.TotalClusterSummary.SpongeCluster;
        for (int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
        {
            int Clusterforpoint = Program.ClusterAssignments[GlobalPointIndex];


            if(Clusterforpoint < 0 ) continue;
            //if(Clusterforpoint == SpongeClusterNumber ) continue;

            int expt =  Program.ExperimentNumberAssigments[GlobalPointIndex] - 1;
            if(expt < 0 ) DAVectorUtility.SALSAPrint(0, "Error: experiment number cannot be less than 0");
            if(expt >= NumberofExperiments) DAVectorUtility.SALSAPrint(0, "Error: experiment number cannot be greater than " + NumberofExperiments);

            if(Clusterforpoint == SpongeClusterNumber ) {
                ++SpongeExptPoints[expt];
                continue;
            }

            ExperimentPoints[expt]++;

            double xpoint = GoldenExamination.PeakPosition[GlobalPointIndex][0];
            double ypoint = GoldenExamination.PeakPosition[GlobalPointIndex][1];
            double xcenter = ClusteringSolution.TotalClusterSummary.CenterPosition[Clusterforpoint][0];
            double ycenter = ClusteringSolution.TotalClusterSummary.CenterPosition[Clusterforpoint][1];
            double Sigmax = Program.SigmaVectorParameters_i_[0] * xcenter;
            double Sigmay = Program.SigmaVectorParameters_i_[1];
            double tmp = (xpoint-xcenter)/Sigmax;
            xshift[expt] += tmp;
            xwidth[expt] += tmp*tmp;
            tmp = (ypoint-ycenter)/Sigmay;
            yshift[expt] += tmp;
            ywidth[expt] += tmp*tmp;
        }

        // output  deviations by experiment number
        for(int expt = 0; expt < NumberofExperiments; expt++) {
            int numberofpointsinexpt = ExperimentPoints[expt];
            if (numberofpointsinexpt == 0) continue;
            double MZshift = xshift[expt] / numberofpointsinexpt;
            double MZSD = Math.sqrt((xwidth[expt] / numberofpointsinexpt) - MZshift * MZshift);
            double RTshift = yshift[expt] / numberofpointsinexpt;
            double RTSD = Math.sqrt((ywidth[expt] / numberofpointsinexpt) - RTshift * RTshift);

            DAVectorUtility.SALSAPrint(0,expt + " Number of Points " +  numberofpointsinexpt + " in Sponge " + SpongeExptPoints[expt]);
            DAVectorUtility.SALSAPrint(0," m/z Normalized shift " + String.format("%1$5.5f", MZshift) + " width " + String.format("%1$5.5f", MZSD));
            DAVectorUtility.SALSAPrint(0," RT Normalized shift " + String.format("%1$5.5f", RTshift) +  " width " + String.format("%1$5.5f", RTSD));

            String file = Program.config.SummaryFile.replace("summary.txt","experimentalShifts.cvs");
            DAVectorUtility.writeExperimentalShifts(file,expt + 1, MZshift,MZSD,RTshift,RTSD);
        }

    }
	public static void CaculateTemperatureClusterCountPlot()
	{
		//  Output Global Counts versus Temperature
		int TemperatureCount = DAVectorUtility.TemperatureValues.size();
		int step = TemperatureCount / 500;
		step = Math.max(1, step);
		int loopnumber = 1 + (TemperatureCount - 1) / step;
		int Tindex = 0;
		String OutputLine = "\nCluster Count v Temperature\n";
		for (int loop = 0; loop < loopnumber; loop++)
		{
			if (Tindex >= TemperatureCount)
			{
				Tindex = TemperatureCount - 1;
			}
			double value = (Double)DAVectorUtility.TemperatureValues.get(Tindex);
			OutputLine += String.format("%1$5.4E", value) + ", ";
			Tindex += step;
		}
		OutputLine += "\n";
		Tindex = 0;
		for (int loop = 0; loop < loopnumber; loop++)
		{
			if (Tindex >= TemperatureCount)
			{
				Tindex = TemperatureCount - 1;
			}
			int ivalue = (Integer)DAVectorUtility.ClusterCountValues.get(Tindex);
			OutputLine += ivalue + ", ";
			Tindex += step;
		}
		DAVectorUtility.SALSAPrint(0, OutputLine);
	}

} // End Class ClusterQuality