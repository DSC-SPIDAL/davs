package edu.indiana.soic.spidal.davs;

import com.google.common.io.Files;
import edu.indiana.soic.spidal.general.Box;
import edu.indiana.soic.spidal.general.IntArray;
import edu.indiana.soic.spidal.mpi.MPIPacket;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

// Calculate Vector Deterministic Annealing Algorithm
public class VectorAnnealIterate
{

	public static int[] ListofClusterstoSplit; // List of Clusters to Split
	public static double[] EigsofClusterstoSplit; // List of Eigenvalues of Clusters to Split
	public static int[] PrioritiesofClusterstoSplit; // Priorities
	public static int Numberthatcanbesplit = 0; // Length of entries in ListofClusterstoSplit which is at most Program.MaxNumberSplitClusters

	public static double Tmin; //minimal temperature
	public static int ActualMaxNcent; // Maximum reduced if small clusters
	public static int CountValidityFailures = 0; // Count failures in validity checks of final solution
	public static boolean OnLastleg = false; // If True we are on final convergence step with full number of clusters
	public static boolean ArithmeticError = false; // Set True if overflow
	public static double MeanClusterCount = 0.0; // Set to Mean Cluster Count
	public static double LocalUselessCalcs = 0.0; // Useless Calcs this iteration
	public static double LocalUsefulCalcs = 0.0; // Useful Calcs this iteration
	public static double LocalIgnoredCalcs = 0.0; // Ignored Calcs this iteration
	public static double PointswithClusterCount1 = 0.0; // Number of Clusters with Cluster Count 1
	public static double C_k_Sum = 0.0; // Set to Sum of C(k)'s for debugging
	public static double MalphaDiffAvg = 0.0; // M value change for test
	public static double AverageY_k_SquaredChange = 0.0; // Average over Clusters of Squared Changes in Y
	public static double NumberCountChanges = 0.0; // Average over points of # Cluster Changes
	public static double TemperatureatLastCloseTest = -1.0;

	public static int EMIterationCount = 0; // iteration number in EM method
	public static int EMIterationStepCount = 0;
	public static int EMIterationStepCount1 = -1; // M convergence
	public static int EMIterationStepCount2 = -1; // Y convergence
	public static int Extra_EMIterationCount = 0; // Extra iterations in EM method
	public static int ActualWaititerations = 0;
	public static int countAfterFixingClusterCount = 0; // Count loops after end on number of centers
	public static int IterationNumberPrinted = -1; // Iteration Number Printed
	public static int convergence = 0; // If zero NOT converged; = 1 Converged but more to go; =2 Converged at end
	public static int CountBetweenSplits = 0; // Count actual iterations between splits. Limit WaitIterations

	public static boolean HitConvergenceLoopLimit = false;
	public static boolean FinalLoop = false; // If true just converging -- No splits
	public static int SplitFailures = 0; // Counts splitting failures
	public static boolean DistributeNextTime = false;
	public static boolean initialized = false;
	public static boolean HammyNotSet = true; // If True Hamiltonian Set
	public static double ChangeinHammy = 0.0; // Change in Hamiltonian

	public static boolean StopReason1 = false;
	public static boolean StopReason2 = false;
	public static boolean StopReason3 = false;
	public static boolean StopReason4 = false;
	public static boolean StopReason5 = false;
	public static boolean StopReason6 = false;

	/*
	 * The basic function in the program. It implements all steps of the Deterministic Annealing Vector Clustering algorithm.
	 *  
	 */
	public final void ControlVectorSpongeDA() throws MPIException {

		//allocate memory on first and indeed only call
		if (!initialized)
		{
			if (!Program.ContinuousClustering)
			{
				DAVectorUtility.printAndThrowRuntimeException(" Invalid Continuous Clustering");

			}

			initialized = true;
			int cachelinesize = Program.cachelinesize;

			ParallelClustering.runningSolution = new ClusteringSolution(Program.UseSponge);
			ParallelClustering.savedSolution = new ClusteringSolution(Program.UseSponge);
			ParallelClustering.bestSolution = new ClusteringSolution(Program.UseSponge);
			DAVectorUtility.SALSAPrint(0, "Clustering Solutions Created");

			Program.InitialNcent = 1; //  We only support starting with one center (plus sponge if needed)
			if (Program.UseSponge)
			{
				++Program.InitialNcent;
			}

			Program.MaxNumberSplitClusters = Math.max(Program.MaxNumberSplitClusters, 1);
			VectorAnnealIterate.ListofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
			VectorAnnealIterate.PrioritiesofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
			VectorAnnealIterate.EigsofClusterstoSplit = new double[Program.MaxNumberSplitClusters];

		} //end Initialization

		//  Do EM calculation
		VectorAnnealIterate.ActualMaxNcent = Program.maxNcentperNode;
		VectorAnnealIterate.OnLastleg = false;

		//  Set up Triangle Inequality
		if (Program.UseTriangleInequality_DA > 0)
		{
			DAVectorUtility.SALSAPrint(0, "Triangle Inequality Initialized Option " + Program.UseTriangleInequality_DA);
			DATriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_DA, Program.MaxClusterLBsperPoint_DA, Program.MaxCentersperCenter_DA, Program.TriangleInequality_Delta1_old_DA, Program.TriangleInequality_Delta1_current_DA, Program.OldCenterOption_DA);

			DATriangleInequality.InitializeTriangleInequality(Program.PointPosition, ParallelClustering.runningSolution.Y_k_i_, ParallelClustering.runningSolution.Sigma_k_i_, ParallelClustering.runningSolution.LocalStatus, ClusteringSolution.RealClusterIndices, ParallelClustering.runningSolution, Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_DA, Program.ParameterVectorDimension);

			DATriangleInequality.SetExternalFunctions(ParallelClustering.runningSolution::SetClusterRadius);
		}

		// Initialize Clusters
		if (Program.RestartTemperature > 0.0)
		{
			// Restart from Previous Runs
			VectorAnnealIterate.Restart();
			ClusteringSolution.CopySolution(ParallelClustering.runningSolution, ParallelClustering.bestSolution);

			//  Set up initial cluster and Sponge if it exists
			if (Program.UseTriangleInequality_DA > 0)
			{
				for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
				{
					DATriangleInequality.AddCenter(ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]);
				}
			}
			DistributedClusteringSolution.ManageMajorSynchronization(true);
			VectorAnnealIterate.ActualWaititerations = 10 * Program.Waititerations;
			ParallelClustering.runningSolution.ActualCoolingFactor = Program.FineCoolingFactor1;
			if (ParallelClustering.runningSolution.Temperature < Program.CoolingTemperatureSwitch)
			{
				ParallelClustering.runningSolution.ActualCoolingFactor = Program.FineCoolingFactor2;
			}
		}
		else
		{
			// Set Initial Temperature Cluster Center P_k_ and occupations Malpha_k_
			VectorAnnealIterate.InitializeSolution(ParallelClustering.runningSolution);

			//  Set up initial cluster and Sponge if it exists
			if (Program.UseTriangleInequality_DA > 0)
			{
				for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
				{
					DATriangleInequality.AddCenter(ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex]);
				}
			}

			DistributedClusteringSolution.ManageMajorSynchronization(true);
			VectorAnnealIterate.ActualWaititerations = Program.Waititerations; // Wait this number of Temperature iterations before splitting
			// This is changed in special circumstances where more convergence needed
		}


		VectorAnnealIterate.EMIterationCount = 0; // Initialize iteration count -- there is no limit
		VectorAnnealIterate.EMIterationStepCount = 0; // This is number of counts in current step converging EM at given temperature
		VectorAnnealIterate.EMIterationStepCount1 = -1;
		VectorAnnealIterate.EMIterationStepCount2 = -1;
		VectorAnnealIterate.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached (this decreases temperature)
		int HammyViolations = 0; // Counts number of (illegal) consecutive INCREASES in Hamiltonian
		VectorAnnealIterate.CountBetweenSplits = 0; // Count iterations between splits. Limit WaitIterations
		boolean CurrentJobFinished = false;

		//  Variables tracking convergence
		VectorAnnealIterate.HitConvergenceLoopLimit = false; // If true, EMIterationStepCount has hit limit
		VectorAnnealIterate.SplitFailures = 0;
		VectorAnnealIterate.FinalLoop = false; // If true we are in a special loop freezing current case

		VectorAnnealIterate.HammyNotSet = true; // To note that OldHammy not set

		ClusteringSolution.ClustersDeleted = 0;
		ClusteringSolution.ClustersMoved = 0;
		ClusteringSolution.ClustersSplit = 0;

		//  Proper Deterministic Annealing
		//	Loop over EM calculations
		while (Program.Refinement)
		{
			VectorAnnealIterate.EMIterationCount++; // Increment EM Loop Count
			VectorAnnealIterate.EMIterationStepCount++; // Iterate within a fixed temperature converging M and p

			//  Set Cooling Factors
			Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
			Program.FineCoolingFactor = Program.FineCoolingFactor1;
			if (ParallelClustering.runningSolution.Temperature < Program.CoolingTemperatureSwitch)
			{
				Program.InitialCoolingFactor = Program.InitialCoolingFactor2;
				Program.FineCoolingFactor = Program.FineCoolingFactor2;
			}

			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - int - min
//                ParallelClustering.runningSolution.DiffMalpha_k_Set = DAVectorUtility.MPI_communicator.<Integer>Allreduce(ParallelClustering.runningSolution.DiffMalpha_k_Set, Operation<Integer>.Min);
                ParallelClustering.runningSolution.DiffMalpha_k_Set = DAVectorUtility.mpiOps.allReduce(ParallelClustering.runningSolution.DiffMalpha_k_Set, MPI.MIN);

            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

			DAVectorEMIterate(); // Basic Update loop
			ParallelClustering.runningSolution.SolutionSet = true;
			ParallelClustering.runningSolution.IterationSetAt = VectorAnnealIterate.EMIterationCount;

			//	Now see if we are done
			String ReasonforConvergence = null;
			double MvalueChange = Program.Malpha_MaxChange;
			if (VectorAnnealIterate.FinalLoop)
			{
				MvalueChange = Program.Malpha_MaxChange1;
			}
			Box<String> tempRef_ReasonforConvergence = new Box<>(ReasonforConvergence);
			convergence = VectorAnnealIterate.convergenceTest(MvalueChange, tempRef_ReasonforConvergence);
			ReasonforConvergence = tempRef_ReasonforConvergence.content;

			convergence = DAVectorUtility.SynchronizeMPIInteger(convergence);

			//  If not converged at this temperature and cluster number just proceed with while(true) iteration
			if (convergence == 0)
			{
				continue;
			}

			if (VectorAnnealIterate.EMIterationStepCount1 == -2)
			{
				++Program.NumberMfailures;
			}
			if (VectorAnnealIterate.EMIterationStepCount1 >= 0)
			{
				Program.AccumulateMvalues += VectorAnnealIterate.EMIterationStepCount1;
				++Program.NumberMsuccesses;
			}
			VectorAnnealIterate.EMIterationStepCount1 = -1;
			if (VectorAnnealIterate.EMIterationStepCount2 == -2)
			{
				++Program.NumberYfailures;
			}
			if (VectorAnnealIterate.EMIterationStepCount2 >= 0)
			{
				Program.AccumulateYvalues += VectorAnnealIterate.EMIterationStepCount2;
				++Program.NumberYsuccesses;
			}
			VectorAnnealIterate.EMIterationStepCount2 = -1;

			// Case we are converged = 1 in middle or = 2 at end of refinement
			Program.NumberIterationSteps++; // Accumulate diagnostic statistics
			Program.IterationsperStepSum += VectorAnnealIterate.EMIterationStepCount;

			//  Update Solution
			DistributedClusteringSolution.ManageMajorSynchronization(true);

			//  Calculate ClusterSquaredWidth_k_
			//  Used in Print and Calculation of Objective Function and Split decisions
			ParallelClustering.runningSolution.SetClusterWidths();

			VectorAnnealIterate.EMIterationStepCount = 0;
			VectorAnnealIterate.EMIterationStepCount1 = -1;
			VectorAnnealIterate.EMIterationStepCount2 = -1;

			//  Update Hamiltonian and see if decreasing
			boolean decreasing = UpdateHamiltonian();

			//  Case when at end of stage (e.g. given number of clusters or of whole job). This still needs to be iterated to low Temperatures
			if (convergence == 2 || VectorAnnealIterate.FinalLoop)
			{
				if (VectorAnnealIterate.countAfterFixingClusterCount < Program.Iterationatend) //    do Iterationatend iterations after reaching the maximum cluster#
				{
					if (VectorAnnealIterate.countAfterFixingClusterCount == 0)
					{ // First step of "justconverging stage"
						HammyViolations = 0;
						Program.ActualEndTemperature = ParallelClustering.runningSolution.Temperature;
						DAVectorUtility.InterimTiming();
						Program.TimeatSplittingStop = DAVectorUtility.HPDuration;
						Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
						ParallelClustering.runningSolution.ActualCoolingFactor = Program.InitialCoolingFactor;
					}

					//  Check Freezing measure -- only use to stop if Y_k_ and P(k) converged
					int toobigfreezing01 = 0;
					int toobigfreezing2 = 0;
					double freezemax01 = 0.0;
					double freezemax2 = 0.0;
					for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
					{
						if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0 || ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 2)
						{
							continue;
						}
						if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
						{
							freezemax2 = Math.max(freezemax2, ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex]);
							if (ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] > Program.FreezingLimit)
							{
								++toobigfreezing2;
							}
						}
						else
						{
							freezemax01 = Math.max(freezemax01, ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex]);
							if (ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] > Program.FreezingLimit)
							{
								++toobigfreezing01;
							}
						}
					}

					if (ParallelClustering.runningSolution.DistributedExecutionMode)
					{
                        freezemax2 = toobigfreezing2;
						DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                        if (DAVectorUtility.MPI_Size > 1){
                            // Note - MPI Call - Allreduce - int - max
//                            freezemax2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(toobigfreezing2, Operation<Double>.Max);
                            freezemax2 = DAVectorUtility.mpiOps.allReduce(toobigfreezing2, MPI.MAX);
                            // Note - MPI Call - Allreduce - int - sum
//                            toobigfreezing2 = DAVectorUtility.MPI_communicator.<Integer>Allreduce(toobigfreezing2, Operation<Integer>.Add);
                            toobigfreezing2 = DAVectorUtility.mpiOps.allReduce(toobigfreezing2, MPI.SUM);
                        }
						DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
					}
					double freezemax = Math.max(freezemax01, freezemax2);

					toobigfreezing01 = DAVectorUtility.SynchronizeMPIInteger(toobigfreezing01);

					int toobigfreezing = toobigfreezing01 + toobigfreezing2;
					if (toobigfreezing == 0)
					{
						DAVectorUtility.SALSAPrint(1, " ****** Stop as Freezing Measures small " + ParallelClustering.runningSolution.FreezingMeasure_k_[0]);
					}

					boolean TemperatureLimit = (ParallelClustering.runningSolution.Temperature < VectorAnnealIterate.Tmin);
					TemperatureLimit = DAVectorUtility.SynchronizeMPIBoolean(TemperatureLimit);
					VectorAnnealIterate.countAfterFixingClusterCount++;
					if (decreasing && (!TemperatureLimit))
					{
						VectorAnnealIterate.countAfterFixingClusterCount = Math.min(Program.Iterationatend - 10, VectorAnnealIterate.countAfterFixingClusterCount);
					}
					String looptype = "Final Clean Up";
					int printtype = -1;
					if ((VectorAnnealIterate.countAfterFixingClusterCount == (Program.Iterationatend - 1)) || (toobigfreezing == 0) || (!decreasing))
					{
						printtype = ParallelClustering.runningSolution.IterationSetAt;
					}
					looptype += " Frz Viol# " + toobigfreezing + " Max " + String.format("%1$4.3E", freezemax) + " " + ReasonforConvergence;
					PrintIteration(looptype, printtype);
					if (decreasing)
					{
						HammyViolations = 0;
					}
					else
					{
						++HammyViolations;
						if (HammyViolations < 5)
						{
							decreasing = true;
						}
					}
					decreasing = DAVectorUtility.SynchronizeMPIBoolean(decreasing);
					if (DAVectorUtility.MPI_Rank == 0)
					{
						CurrentJobFinished = !((toobigfreezing > 0) && decreasing);
						if (VectorAnnealIterate.ArithmeticError)
						{
							CurrentJobFinished = true;
						}
					}
					CurrentJobFinished = DAVectorUtility.SynchronizeMPIBoolean(CurrentJobFinished);

					// Iteration at end still going Reduce Temperature and proceed
					if (!CurrentJobFinished)
					{
						ReduceTemperature();
						continue; // Continue while over EM Loop
					}

					// We have finished this job
					if (VectorAnnealIterate.ArithmeticError)
					{
						StopReason1 = true;
					}
					if (!decreasing)
					{
						StopReason2 = true;
					}
					if (toobigfreezing > 0)
					{
						StopReason3 = true;
					}

				} // End Convergence=2  or justconverging doing Final Program.Iterationatend iterations counted by countAfterFixingClusterCount

				else
				{ // Case when Program.Iterationatend iterations counted by countAfterFixingClusterCount exceeded
					StopReason6 = true;
					CurrentJobFinished = true;
				} // End Convergence=2 or justconverging case where extra iteration count completed


				// This completes processing for convergence=2 and justconverging=true case
				if (CurrentJobFinished)
				{
					if (ParallelClustering.runningSolution.Temperature < VectorAnnealIterate.Tmin)
					{
						StopReason4 = true;
					}
					if (VectorAnnealIterate.SplitFailures > 1)
					{
						StopReason5 = true;
					}


					// Real end of job except check validity
					DAVectorUtility.StartSubTimer(1);
					boolean validity = VectorAnnealIterate.CheckValidSolution(true);
					DAVectorUtility.StopSubTimer(1);
					if ((!validity) || (!VectorAnnealIterate.FinalLoop))
					{ // Restart with small clusters removed or do a final improved precision run
						if (!validity)
						{
							DAVectorUtility.SALSAPrint(1, "Restart after small clusters removed");
						}
						if (!VectorAnnealIterate.FinalLoop)
						{
							DAVectorUtility.SALSAPrint(1, "Add in Y convergence Test");
						}
						VectorAnnealIterate.ActualMaxNcent = ParallelClustering.runningSolution.Ncent_ThisNode;
						ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;
						ParallelClustering.runningSolution.YPreviousSet = -1;
						VectorAnnealIterate.OnLastleg = true;
						VectorAnnealIterate.countAfterFixingClusterCount = Math.max(1, Program.Iterationatend - 100);
						VectorAnnealIterate.CountBetweenSplits = 0;
						VectorAnnealIterate.HammyNotSet = true;
						HammyViolations = 0;
						VectorAnnealIterate.FinalLoop = true;
						CurrentJobFinished = false;
						continue;
					}

					if (StopReason1)
					{
						Program.FinalReason += "Arithmetic Error ";
					}
					if (StopReason2)
					{
						Program.FinalReason += "Stop as Hamiltonian Increasing with Change " + String.format("%1$5.4E", ChangeinHammy) + " ";
					}
					if (StopReason3)
					{
						Program.FinalReason += "Stop as Freezing Measures smaller than " + (Program
                                .FreezingLimit) + " ";
					}
					if (StopReason4)
					{
						Program.FinalReason += "Tmin Reached " + String.format("%1$5.4E", VectorAnnealIterate.Tmin) + " ";
					}
					if (StopReason5)
					{
						Program.FinalReason += "Consecutive Split Failures ";
					}
					if (StopReason6)
					{
						Program.FinalReason += " Stop as Final Iteration Count larger than " + Program.Iterationatend + " ";
					}
					DAVectorUtility.SALSAPrint(1, Program.FinalReason);
					break;
				} // end case when task finished

			} // End justconverging = true or convergence =2 case

			DAVectorUtility.StartSubTimer(1);
			boolean convergedvalidity = VectorAnnealIterate.CheckValidSolution(false);
			DAVectorUtility.StopSubTimer(1);
			if (!convergedvalidity)
			{ // Restart with small clusters removed
				PrintIteration(" Restart After Cluster Removal ", -1);
				CountBetweenSplits = 0;
				VectorAnnealIterate.ActualWaititerations = Program.Waititerations_Converge;
				HammyNotSet = true;
				continue;
			}
			//  This section results either in EM Loop Continue (for ongoing refinement for a nonsplittable case) or EM Loop Break (to end) 

			//  Convergence = 1 Case
			//	Converged so test for split -  -- take results from Rank 0 but rest do arithmetic
			//  For restarted task that decision has already been made

			// Need to decide if split to occur
			boolean ResultofSplittingTest = false;
			++CountBetweenSplits;
			if (CountBetweenSplits > VectorAnnealIterate.ActualWaititerations)
			{
				CountBetweenSplits = 0;
				VectorAnnealIterate.ActualWaititerations = Program.Waititerations;
			}
			if ((CountBetweenSplits > 0) && (VectorAnnealIterate.ActualWaititerations > 0) && (ParallelClustering.runningSolution.Ncent_Global > 1))
			{
				PrintIteration("Ongoing Annealing", -1);
				ReduceTemperature();
				continue; // EM Loop for compulsory iterations between splits
			}

			boolean ChangedClusters = VectorAnnealIterate.CheckNumberofClustersTooBig();
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - boolean - logical OR
//                ChangedClusters = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(ChangedClusters, Operation<Boolean>.LogicalOr);
                ChangedClusters = DAVectorUtility.mpiOps.allReduce(ChangedClusters, MPI.LOR);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
			if (ChangedClusters)
			{
				DistributedClusteringSolution.ManageMajorSynchronization(true);
			}

			//  Decide if to split
			DAVectorUtility.StartSubTimer(0);
			ResultofSplittingTest = shouldweSplit();
			if (ParallelClustering.runningSolution.DistributedExecutionMode)
			{
				DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                if (DAVectorUtility.MPI_Size > 1){
                    // Note - MPI Call - Allreduce - boolean - logical OR
//                    DAVectorUtility.MPI_communicator.<Boolean>Allreduce(ResultofSplittingTest, Operation<Boolean>.LogicalOr);
                    ResultofSplittingTest = DAVectorUtility.mpiOps.allReduce(ResultofSplittingTest, MPI.LOR);
                }
				DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
			}
			else
			{
				ResultofSplittingTest = DAVectorUtility.SynchronizeMPIBoolean(ResultofSplittingTest);
                ParallelClustering.runningSolution.ClustertoSplit = DAVectorUtility.SynchronizeMPIInteger(ParallelClustering.runningSolution.ClustertoSplit);
			}
			DAVectorUtility.StopSubTimer(0);

			//  Diagnostic Output for splitting

			if (ResultofSplittingTest && (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0))
			{
				this.diagnosticsplitprint(ResultofSplittingTest);
			}


			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - boolean - logical OR
//                ResultofSplittingTest = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(ResultofSplittingTest, Operation<Boolean>.LogicalOr);
                ResultofSplittingTest = DAVectorUtility.mpiOps.allReduce(ResultofSplittingTest, MPI.LOR);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

			//	If split indicated perform this                          
			if (ResultofSplittingTest)
			{
				if (ParallelClustering.runningSolution.Ncent_Global > 1)
				{
					if (Program.ClusterCountOutput > 0)
					{
						VectorAnnealIterate.OutputClusteringResults("Inter2");
					}
				}

				//  Need to perform already determined split (maybe zero in one node if distributed)
				DAVectorUtility.StartSubTimer(5);

				ClusteringSolution.ClustersSplit = 0;
				String SplitString = "";
				if (VectorAnnealIterate.Numberthatcanbesplit > 0)
				{
					for (int splitlistloop = 0; splitlistloop < VectorAnnealIterate.Numberthatcanbesplit; splitlistloop++)
					{
						ParallelClustering.runningSolution.ClustertoSplit = VectorAnnealIterate.ListofClusterstoSplit[splitlistloop];
						SplitString += ParallelClustering.runningSolution.ClustertoSplit + "(" + (new Double(ParallelClustering.runningSolution.C_k_[ParallelClustering.runningSolution.ClustertoSplit])) + ")" + " " + String.format("%1$5.4E", VectorAnnealIterate.EigsofClusterstoSplit[splitlistloop]) + " * ";
						this.dothesplit();

					}
					VectorAnnealIterate.Numberthatcanbesplit = 0;
				}
				if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
				{
					int TotalNumberSplit = ClusteringSolution.ClustersSplit;
					if (ParallelClustering.runningSolution.DistributedExecutionMode)
					{
						DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                        if (DAVectorUtility.MPI_Size > 1){
                            // Note - MPI Call - Allreduce - int - sum
//                            TotalNumberSplit = DAVectorUtility.MPI_communicator.<Integer>Allreduce(TotalNumberSplit, Operation<Integer>.Add);
                            TotalNumberSplit = DAVectorUtility.mpiOps.allReduce(TotalNumberSplit, MPI.SUM);
                        }
						DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
						if (TotalNumberSplit > 0)
						{
							DAVectorUtility.SALSASyncPrint(1, "Clusters Split ", SplitString);
						}
					}
					else
					{
						if (TotalNumberSplit > 0)
						{
							DAVectorUtility.SALSAPrint(1, "Clusters Split " + SplitString);
						}
					}
				}
				DAVectorUtility.StopSubTimer(5);

				DistributedClusteringSolution.ManageMajorSynchronization(true);
				VectorAnnealIterate.SplitFailures = 0;
				VectorAnnealIterate.ActualWaititerations = Program.Waititerations;
			} // End ResultofSplittingTest == true (either changed number of clusters or spun off a convergence task)

			//  Final portion of loop changing Temperature if needed
			if (!ResultofSplittingTest)
			{ //    Reduce T and continue iterating if converged and no split
				ReduceTemperature();
			}
			else
			{ // Don't reduce T as clusters split
				ParallelClustering.runningSolution.ActualCoolingFactor = Program.FineCoolingFactor;
				CountBetweenSplits = 0;
			}

		} // End while EM Loop

		//  Check if solution best
		if (!ParallelClustering.bestSolution.SolutionSet)
		{
			return;
		}
		if (ParallelClustering.bestSolution.Ncent_Global != ParallelClustering.runningSolution.Ncent_Global)
		{
			DAVectorUtility.SALSAPrint(1, " Best Solution Not Used as Number of Centers " + ParallelClustering.bestSolution.Ncent_Global + " Different from Running Solution with " + ParallelClustering.runningSolution.Ncent_Global);
			return;
		}
		boolean changesolution = ParallelClustering.bestSolution.PairwiseHammy < ParallelClustering.runningSolution.PairwiseHammy;
		changesolution = DAVectorUtility.SynchronizeMPIBoolean(changesolution);
		if (changesolution)
		{
			DAVectorUtility.SALSAPrint(1, " Solution at Iteration " + ParallelClustering.bestSolution.IterationSetAt + " Chisq " + String.format("%1$5.4E", ParallelClustering.bestSolution.PairwiseHammy) + " Taken rather than Iteration " + ParallelClustering.runningSolution.IterationSetAt + " Chisq " + String.format("%1$5.4E", ParallelClustering.runningSolution.PairwiseHammy));
			ClusteringSolution.CopySolution(ParallelClustering.bestSolution, ParallelClustering.runningSolution);
		}
		else if (VectorAnnealIterate.ArithmeticError)
		{
			DAVectorUtility.SALSAPrint(1, " Solution at Iteration " + ParallelClustering.bestSolution.IterationSetAt + " Taken rather than Iteration " + ParallelClustering.runningSolution.IterationSetAt + " Due to Arithmetic Error ");
			ClusteringSolution.CopySolution(ParallelClustering.bestSolution, ParallelClustering.runningSolution);
		}
		int save = Program.ClusterPrintNumber;
		Program.ClusterPrintNumber = ParallelClustering.runningSolution.Ncent_ThisNode;
		PrintIteration(" Final Solution ", ParallelClustering.runningSolution.IterationSetAt);
		Program.ClusterPrintNumber = save;

	} // End of ControlVectorSpongeDA()

	public final void ReduceTemperature() throws MPIException {
		// Add in a Sponge if needed
		if (VectorAnnealIterate.AddSpongeCluster())
		{
			ParallelClustering.bestSolution.SolutionSet = false;
			CountBetweenSplits = 0;
			HammyNotSet = true;
			return;
		}
		else
		{
			if (Program.ChangeSpongeFactor(ParallelClustering.runningSolution.ActualCoolingFactor,
                    ParallelClustering.runningSolution))
			{
				HammyNotSet = true;
			}
		}
		//  Reduce Cluster Sigmas for those being Annealed
		if (Program.ChangeClusterSigmas(ParallelClustering.runningSolution.ActualCoolingFactor,
                ParallelClustering.runningSolution))
		{
			ParallelClustering.runningSolution.SetClusterSizes();
			ParallelClustering.runningSolution.SetClusterWidths();
			HammyNotSet = true;
		}

		//  Switch into Distributed Mode if necessary
		if (!ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			boolean godistributed = (ParallelClustering.runningSolution.Temperature <= Program.TemperatureLimitforDistribution) || ((Program.ClusterLimitforDistribution > 0) && (ParallelClustering.runningSolution.Ncent_Global >= Program.ClusterLimitforDistribution));
			godistributed = DAVectorUtility.SynchronizeMPIBoolean(godistributed);
			if (godistributed)
			{
				if (VectorAnnealIterate.DistributeNextTime)
				{
					ParallelClustering.runningSolution.DistributedExecutionMode = true;

					int newsplitnumber = Program.MaxNumberSplitClusters / DAVectorUtility.MPI_Size;
					if ((newsplitnumber * Program.MaxNumberSplitClusters) < Program.MaxNumberSplitClusters)
					{
						++newsplitnumber;
					}
					Program.MaxNumberSplitClusters = newsplitnumber;
					++DAVectorUtility.MPIREDUCETiming1;
					Program.ActualTemperatureforDistribution = ParallelClustering.runningSolution.Temperature;
					DAVectorUtility.InterimTiming();
					Program.TimeatDistribution = DAVectorUtility.HPDuration;
					Program.ActualClusterNumberforDistribution = ParallelClustering.runningSolution.Ncent_Global;
					DAVectorUtility.SALSAPrint(1, "Start Distributed Execution Mode " + String.format("%1$3.2f", ParallelClustering.runningSolution.Temperature) + " Iteration Count " + ParallelClustering.runningSolution.IterationSetAt + " " + EMIterationCount + " Clusters " + ParallelClustering.runningSolution.Ncent_Global + " Time " + String.format("%1$5.4E", DAVectorUtility.HPDuration));
					DistributedClusteringSolution.ManageMajorSynchronization(false);
					ParallelClustering.bestSolution.SolutionSet = false;
					CountBetweenSplits = 0;
					HammyNotSet = true;
					return;
				}
				else
				{
					VectorAnnealIterate.DistributeNextTime = true;
					VectorAnnealIterate.CompleteCleanUp();
					return;
				}
			}
		}

		//  Generate Magic Temperature Actions
		if (Program.magicindex < Program.MagicTemperatures.length)
		{
			boolean abracadabra = ParallelClustering.runningSolution.Temperature <= Program.MagicTemperatures[Program.magicindex];
			abracadabra = DAVectorUtility.SynchronizeMPIBoolean(abracadabra);
			if (abracadabra)
			{
//                    Dist.OutputClusterLabels("Temp" + ParallelClustering.RunningDAVectorSolt.Temperature.ToString("F2"));
				VectorAnnealIterate.CompleteCleanUp();
				if (ParallelClustering.runningSolution.DistributedExecutionMode)
				{
					DistributedClusteringSolution.ManageMajorSynchronization(true);
				}
				CountBetweenSplits = 0;
				++Program.magicindex;
				return;
			}
		}
		ParallelClustering.runningSolution.Temperature = ParallelClustering.runningSolution.ActualCoolingFactor * ParallelClustering.runningSolution.Temperature;
		++Program.NumberTemperatureSteps;
		DAVectorUtility.TemperatureValues.add(ParallelClustering.runningSolution.Temperature);
		DAVectorUtility.ClusterCountValues.add(ParallelClustering.runningSolution.Ncent_Global);

	} // End ReduceTemperature

	public final void PrintIteration(String looptype, int linecheck)
	{
		if (linecheck < 0)
		{
			if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval != 0)
			{
				return;
			}
		}
		else if (linecheck == VectorAnnealIterate.IterationNumberPrinted)
		{
			return;
		}
		VectorAnnealIterate.IterationNumberPrinted = ParallelClustering.runningSolution.IterationSetAt;

		DAVectorUtility.InterimTiming();
		String endinfo = "";
		if (!Program.CalculateIndividualWidths)
		{
			endinfo = " Average Width " + String.format("%1$4.3E", ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
		}
		else
		{
			endinfo = " Average Widths ";
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				endinfo += String.format("%1$4.3E", ParallelClustering.runningSolution.AverageWidth[VectorIndex]) + " ";
			}
		}
		if (Program.SigmaMethod == 3)
		{
			endinfo += " Sigma[0] Coeff " + String.format("%1$5.4E", Program.SigmaVectorParameters_i_[0]);
		}

		DAVectorUtility.SALSAPrint(1, "B) Clusters " + ParallelClustering.runningSolution.Ncent_Global + " " + looptype + " Cluster # Chg " + String.format("%1$4.3E", VectorAnnealIterate.NumberCountChanges) + " M-Change " + String.format("%1$4.3E", VectorAnnealIterate.MalphaDiffAvg) + " Y Chg " + String.format("%1$4.3E", VectorAnnealIterate.AverageY_k_SquaredChange / ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth) + " Cnvg " + VectorAnnealIterate.convergence + " Iter " + VectorAnnealIterate.EMIterationCount + " Major " + Program.NumberMajorSynchs1 + " T " + String.format("%1$5.4E", ParallelClustering.runningSolution.Temperature) + " PWHammy " + String.format("%1$5.4E", ParallelClustering.runningSolution.PairwiseHammy) + " Useful Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalUsefulCalcs) + " Useless Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalUselessCalcs) + " Ignored Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalIgnoredCalcs) + " Arithmetic Error " + (Boolean.valueOf(
                VectorAnnealIterate.ArithmeticError)) + " Mean Cluster Count per point " + String.format("%1$3.2f", VectorAnnealIterate.MeanClusterCount) + " Pts with Just 1 Cluster " + String.format("%1$5.4E", VectorAnnealIterate.PointswithClusterCount1) + " Sum of C(k) " + String.format("%1$3.2f", VectorAnnealIterate.C_k_Sum) + endinfo + " Time " + String.format("%1$5.4E", DAVectorUtility.HPDuration));

		clusterprint();

	} // End Print Iteration from getDist

	public final void diagnosticsplitprint(boolean ResultofSplittingTest)
	{
		String endinfo = "";
		if (!Program.CalculateIndividualWidths)
		{
			endinfo = " Average Width " + String.format("%1$4.3E", ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
		}
		else
		{
			endinfo = " Average Widths ";
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				endinfo += String.format("%1$4.3E", ParallelClustering.runningSolution.AverageWidth[VectorIndex]) + " ";
			}
		}
		if (Program.SigmaMethod == 3)
		{
			endinfo += " Sigma[0] Coeff " + String.format("%1$5.4E", Program.SigmaVectorParameters_i_[0]);
		}
		DAVectorUtility.InterimTiming();
		String nextline1 = "A) Clusters " + ParallelClustering.runningSolution.Ncent_Global + " Iter " + VectorAnnealIterate.EMIterationCount + " Major " + Program.NumberMajorSynchs1 + " T " + String.format("%1$5.4E", ParallelClustering.runningSolution.Temperature) + " PWHammy " + String.format("%1$5.4E", ParallelClustering.runningSolution.PairwiseHammy) + " Useful Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalUsefulCalcs) + " Useless Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalUselessCalcs) + " Ignored Calcs " + String.format("%1$5.4E", VectorAnnealIterate.LocalIgnoredCalcs) + " Arithmetic Error " + (Boolean.valueOf(
                VectorAnnealIterate.ArithmeticError)) + " Mean Cluster Count per point " + String.format("%1$3.2f", VectorAnnealIterate.MeanClusterCount) + " Pts with Just 1 Cluster " + String.format("%1$5.4E", VectorAnnealIterate.PointswithClusterCount1) + " Sum of C(k) " + String.format("%1$3.2f", VectorAnnealIterate.C_k_Sum) + endinfo + " Time " + String.format("%1$5.4E", DAVectorUtility.HPDuration);
		if (ParallelClustering.runningSolution.ClustertoSplit >= 0)
		{
			if (ResultofSplittingTest)
			{
				nextline1 += " C# To Split " + ParallelClustering.runningSolution.ClustertoSplit;
			}
			else
			{
				nextline1 += " No Split";
			}
			nextline1 += " " + (Boolean.valueOf(ResultofSplittingTest)) + " Status " + ParallelClustering.runningSolution.Splittable_k_[ParallelClustering.runningSolution.ClustertoSplit] + " Eig " + String.format("%1$5.4E", ParallelClustering.runningSolution.Eigenvalue_k[ParallelClustering.runningSolution.ClustertoSplit]);
		}
		DAVectorUtility.SALSAPrint(1, nextline1);
		clusterprint();
		DAVectorUtility.SALSAPrint(1, " ");

    } // End diagnosticsplitprint

	public final void clusterprint()
	{
		String nextline = "";
		int ClusterIndex = ParallelClustering.runningSolution.SpongeCluster;
		if (ClusterIndex < 0)
		{
			ClusterIndex = 0;
		}
		int count = 0;
		while (true)
		{
			String spongelabel = "";
			String Pformat = "F4";
			if (ParallelClustering.runningSolution.LocalStatus[ClusterIndex] >= 0)
			{
				if (ClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
				{
					if (count > 0)
					{
						ClusterIndex++;
						continue;
					}
					spongelabel = "Sponge " + String.format("%1$3.2f", Program.SpongeFactor) + " ";
					Pformat = "E3";
				}
				double tmp = ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[ClusterIndex];
				if (count != 0)
				{
					nextline += "* ";
				}
				nextline += spongelabel + ClusterIndex + "(" + ParallelClustering.runningSolution.LocalCreatedIndex[ClusterIndex] + ") C " + String.format("%1$2.1f", ParallelClustering.runningSolution.C_k_[ClusterIndex]) + " (Frz " + String.format("%1$7.6f", ParallelClustering.runningSolution.FreezingMeasure_k_[ClusterIndex]) + ") " + "[Wdth " + String.format("%1$5.4f", tmp) + "] ";
				if (Program.ContinuousClustering)
				{
					nextline += "P " + String.format(Pformat, ParallelClustering.runningSolution.P_k_[ClusterIndex]) + " ";
				}
			}
			if ((count >= Program.ClusterPrintNumber) || (ClusterIndex >= ParallelClustering.runningSolution.Ncent_ThisNode-1))
			{
				break;
			}
			count++;
			if ((count == 1) && (ParallelClustering.runningSolution.SpongeCluster >= 0))
			{
				ClusterIndex = 0;
			}
			else
			{
				ClusterIndex++;
			}
		}
		int ClusterLimit = Math.max(ClusterIndex + 1, ParallelClustering.runningSolution.Ncent_ThisNode - Program.ClusterPrintNumber);
		for (int ClusterEndIndex = ClusterLimit; ClusterEndIndex < ParallelClustering.runningSolution.Ncent_ThisNode; ClusterEndIndex++)
		{
			if (ClusterEndIndex == ParallelClustering.runningSolution.SpongeCluster)
			{
				continue;
			}
			if (ParallelClustering.runningSolution.LocalStatus[ClusterEndIndex] < 0)
			{
				continue;
			}
			double tmp = ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[ClusterEndIndex];
			if (ClusterEndIndex != 0)
			{
				nextline += "* ";
			}
			nextline += ClusterEndIndex + "(" + ParallelClustering.runningSolution.LocalCreatedIndex[ClusterEndIndex] + ") C " + String.format("%1$2.1f", ParallelClustering.runningSolution.C_k_[ClusterEndIndex]) + " (Frz " + String.format("%1$7.6f", ParallelClustering.runningSolution.FreezingMeasure_k_[ClusterEndIndex]) + ") " + "[Wdth " + String.format("%1$5.4f", tmp) + "] ";
			if (Program.ContinuousClustering)
			{
				nextline += "P " + String.format("%1$5.4f", ParallelClustering.runningSolution.P_k_[ClusterEndIndex]) + " ";
			}
		}
		DAVectorUtility.SALSAPrint(1, nextline);

	} // End clusterprint()

	//  Update Hamiltonian and return boolean decreasing to indicate if decreasing
	public final boolean UpdateHamiltonian() throws MPIException {
		double hammy2 = 0.0;
		double hammy01 = 0.0;
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0 || ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 2)
			{
				continue;
			}
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
			{
				hammy2 += 0.5 * ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.runningSolution.C_k_[RealClusterIndex];
			}
			else
			{
				hammy01 += 0.5 * ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.runningSolution.C_k_[RealClusterIndex];
			}
		}
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - double - sum
//                hammy2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(hammy2, Operation<Double>.Add);
                hammy2 = DAVectorUtility.mpiOps.allReduce(hammy2, MPI.SUM);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		}
		ParallelClustering.runningSolution.PairwiseHammy = hammy2 + hammy01;

		//  Set Best Solution and progress flag decreasing
		boolean decreasing;
		if (HammyNotSet)
		{
			HammyNotSet = false;
			decreasing = true;
		}
		else
		{
			decreasing = ParallelClustering.runningSolution.PairwiseHammy <= ParallelClustering.runningSolution.OldHammy;
		}
		decreasing = DAVectorUtility.SynchronizeMPIBoolean(decreasing);
		if (decreasing)
		{
			ClusteringSolution.CopySolution(ParallelClustering.runningSolution, ParallelClustering.bestSolution);
		}
		VectorAnnealIterate.ChangeinHammy = ParallelClustering.runningSolution.PairwiseHammy - ParallelClustering.runningSolution.OldHammy;
		ParallelClustering.runningSolution.OldHammy = ParallelClustering.runningSolution.PairwiseHammy;
		return decreasing;

	} // End UpdateHamiltonian()

	//  Process CreatedIndex for a Point Cluster pointer
	// Note ActiveCluster is defined even for Remote Clusters
	public static void ClusterPointersforaPoint(int alpha, int IndirectClusterIndex, Box<Integer> RealClusterIndex, Box<Integer> ActiveClusterIndex, Box<Integer> RemoteIndex)
	{
		int CreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
		int MappedClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
		int ClusterIterationNo = ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet;
		if (ClusterIterationNo < ClusteringSolution.CurrentIteration)
		{
			int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
			String errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + IndirectClusterIndex + " Created Index " + CreatedIndex + " Actual Iteration " + ClusteringSolution.CurrentIteration + " Cluster Iteration " + ClusterIterationNo + " Full Set Created Indices ";
			for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
			{
				errormessage += ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop] + " ";
			}
			DAVectorUtility.printAndThrowRuntimeException(errormessage);

		}
		if (MappedClusterIndex > 0)
		{
			RealClusterIndex.content = MappedClusterIndex - 1;
			if (RealClusterIndex.content >= ParallelClustering.runningSolution.Ncent_ThisNode)
			{
				int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
				String errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + IndirectClusterIndex + " Bad Cluster Number " + RealClusterIndex.content + " Created Index " + CreatedIndex + " Actual Iteration " + ClusteringSolution.CurrentIteration + " Cluster Iteration " + ClusterIterationNo + " Full Set Created Indices ";
				for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
				{
					errormessage += ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop] + " ";
				}
				DAVectorUtility.printAndThrowRuntimeException(errormessage);

			}
			ActiveClusterIndex.content = ClusteringSolution.ActiveClusterIndices[RealClusterIndex.content];
        }
		else if (MappedClusterIndex == 0)
		{
			int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
			String errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + IndirectClusterIndex + " Zero Mapped Index " + " Created Index " + CreatedIndex + " Actual Iteration " + ClusteringSolution.CurrentIteration + " Cluster Iteration " + ClusterIterationNo + " Full Set Created Indices ";
			for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
			{
				errormessage += ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop] + " ";
			}
			DAVectorUtility.printAndThrowRuntimeException(errormessage);

		}
		else
		{
			if (!ParallelClustering.runningSolution.DistributedExecutionMode)
			{
				DAVectorUtility.printAndThrowRuntimeException(" DAVectorEMIterate() Illegal Created Index " + CreatedIndex + " Point " + (alpha + DAVectorUtility.PointStart_Process));

			}
			RemoteIndex.content = -MappedClusterIndex - 1;
			if (RemoteIndex.content >= DistributedClusteringSolution.StorageforTransportedClusters.SizeOfTransportedArray)
			{
				int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
				String errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process) + " Cluster Index " + IndirectClusterIndex + " Remote Index " + RemoteIndex.content + " Created Index " + CreatedIndex + " Actual Iteration " + ClusteringSolution.CurrentIteration + " Cluster Iteration " + ClusterIterationNo + " Full Set Created Indices ";
				for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
				{
					errormessage += ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop] + " ";
				}
				DAVectorUtility.printAndThrowRuntimeException(errormessage);

			}
			ActiveClusterIndex.content = RemoteIndex.content + ClusteringSolution.NumberLocalActiveClusters;
        }

	} // End ClusterPointersforaPoint( int alpha, int IndirectClusterIndex, int IndirectSize, ref int RealClusterIndex, ref int ActiveClusterIndex, ref int RemoteIndex)

    //  Update in EM order Malpha_k_ setting Previous_Malpha_k_, P_k_ Y_k_ and Differences in Malpha for convergence test
    //  Set new values of C_k_, P_k_ Y_k_ and Freezing Factors
    public final void DAVectorEMIterate() throws MPIException {
        DAVectorUtility.StartSubTimer(2);
        final double SpongeTerm;
        if (ParallelClustering.runningSolution.SpongeCluster >= 0)
        {
            SpongeTerm = Program.SpongeFactor * Program.SpongeFactor;
        } else {
            SpongeTerm = 0.0;
        }
        double MinP_k_ = 0.0000001;

        final GlobalReductions.FindIndirectVectorDoubleSum FindDiffMalpha_k_;
        final GlobalReductions.FindIndirectVectorDoubleSum FindC_k_;
        final GlobalReductions.FindIndirectVectorDoubleSum FindFreezingMeasure_k_;
        final DistributedReductions.FindIndirectMultiVectorDoubleSum Find3Components;
        final int BeginFindDiffMalpha_k_;
        final int BeginFindC_k_;
        final int BeginFindFreezingMeasure_k_;

        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        {
            Find3Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginFindDiffMalpha_k_ = Find3Components.AddComponents(1);
            BeginFindC_k_ = Find3Components.AddComponents(1);
            BeginFindFreezingMeasure_k_ = Find3Components.AddComponents(1);
            Find3Components.NodeInitialize();
            FindDiffMalpha_k_ = null;
            FindC_k_ = null;
            FindFreezingMeasure_k_ = null;
        }
        else
        {
            FindDiffMalpha_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            FindC_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            FindFreezingMeasure_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            Find3Components = null;
            BeginFindDiffMalpha_k_ = -1;
            BeginFindC_k_ = -1;
            BeginFindFreezingMeasure_k_ = -1;
        }

        GlobalReductions.FindBoolOr FindArithmeticError = new GlobalReductions.FindBoolOr(DAVectorUtility.ThreadCount);
        GlobalReductions.FindDoubleMean FindMeanClusterCount = new GlobalReductions.FindDoubleMean(DAVectorUtility.ThreadCount);
        GlobalReductions.FindDoubleSum FindSinglePoints = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

        VectorAnnealIterate.LocalUsefulCalcs = 0.0;
        VectorAnnealIterate.LocalUselessCalcs = 0.0;
        VectorAnnealIterate.LocalIgnoredCalcs = 0.0;
        GlobalReductions.FindDoubleSum FindUsefulCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
        GlobalReductions.FindDoubleSum FindUselessCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
        GlobalReductions.FindDoubleSum FindIgnoredCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                    Find3Components.ThreadInitialize(threadIndex);
                } else {
                    FindDiffMalpha_k_.startthread(threadIndex);
                    FindC_k_.startthread(threadIndex);
                    FindFreezingMeasure_k_.startthread(threadIndex);
                }

                //  Arrays used in each thread and re-used by points in thread
                int ArraySize = Math.min(ClusteringSolution.NumberAvailableActiveClusters,
                        ClusteringSolution.MaximumClustersperPoint);
                double[] AbsChangeinMalpha_k_ = new double[ArraySize];
                double[] Malphax1minusMalpha = new double[ArraySize];
                double[] Save_CurrentP_k_TimesExponential = new double[ArraySize];
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] Save_Term_NegativeExponential = new double[ArraySize];
                int[] ThreadStorePosition = new int[ArraySize];

                //  Loop over Points
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int MinimumTerm_IndirectClusterIndex = -1;
                    double MinimumTerm_P_k_ = 0.0;
                    double MinimumTerm_NegativeExponential = 0.0;

                    //  Number of Clusters used for this point
                    int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];

                    //  Accumulate mean number of clusters per point
                    FindMeanClusterCount.addapoint(threadIndex, (double) IndirectSize);

                    //  Set count of points with just one cluster
                    double value = 0.0;
                    if (IndirectSize == 1) {
                        value = 1.0;
                    }
                    FindSinglePoints.addapoint(threadIndex, value);
                    double NumSponge = 0.0;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        // First Pass loop sets up a safer piece of arithmetic by subtracting smallest term in exponential
                        // It will cancel in top and bottom when calculating M
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        Box<Integer> tempRef_RealClusterIndex = new Box<>(RealClusterIndex);
                        Box<Integer> tempRef_ActiveClusterIndex = new Box<>(ActiveClusterIndex);
                        Box<Integer> tempRef_RemoteIndex = new Box<>(RemoteIndex);
                        ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex,
                                tempRef_ActiveClusterIndex, tempRef_RemoteIndex);
                        RealClusterIndex = tempRef_RealClusterIndex.content;
                        ActiveClusterIndex = tempRef_ActiveClusterIndex.content;
                        RemoteIndex = tempRef_RemoteIndex.content;
                        if (RemoteIndex < 0) {
                            ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                            if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                                ThreadStorePosition[IndirectClusterIndex] = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                            }
                        } else {
                            ThreadStorePosition[IndirectClusterIndex] = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                        }
                        double Term_NegativeExponential = 0.0;
                        if ((ParallelClustering.runningSolution.SpongeCluster >= 0) && (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)) {
                            Term_NegativeExponential = SpongeTerm;
                            NumSponge = 1.0;
                        } else {
                            Term_NegativeExponential = DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(
                                    alpha,
                                    ActiveClusterIndex, ParallelClustering.runningSolution);
                        }
                        Term_NegativeExponential = Term_NegativeExponential / ParallelClustering.runningSolution.Temperature;
                        Save_Term_NegativeExponential[IndirectClusterIndex] = Term_NegativeExponential;

                        //  Set P_k_ for this cluster
                        double CurrentP_k_;
                        if (RemoteIndex < 0) {
                            CurrentP_k_ = ParallelClustering.runningSolution.P_k_[RealClusterIndex];
                        } else {
                            CurrentP_k_ = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
                        }
                        Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = CurrentP_k_;

                        //  Set Selected (aka Best) Index to use for subtraction
                        if ((MinimumTerm_IndirectClusterIndex == -1) || (Term_NegativeExponential < MinimumTerm_NegativeExponential)) {
                            MinimumTerm_NegativeExponential = Term_NegativeExponential;
                            MinimumTerm_IndirectClusterIndex = IndirectClusterIndex;
                            MinimumTerm_P_k_ = CurrentP_k_;
                        }
                    } // End First Pass Loop over Clusters for this point

                    boolean LocalArithmeticError = false;
                    double MalphaDenominator = 0.0;
                    boolean ThereisanOKTerm_NegativeExponential = false;

                    //  Second Pass: Now calculate quantities subtracting smallest term at top and bottom
                    double NumUseful = 0.0;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        double ExponentiatedTerm;
                        double SubtractedTerm_NegativeExponential = Save_Term_NegativeExponential[IndirectClusterIndex] - MinimumTerm_NegativeExponential;
                        if (SubtractedTerm_NegativeExponential < (2.0 * Program.ExpArgumentCut1)) {
                            NumUseful += 1.0;
                            ExponentiatedTerm = Math.exp(-0.5 * SubtractedTerm_NegativeExponential);
                            if (Double.isNaN(ExponentiatedTerm) || Double.isInfinite(ExponentiatedTerm)) {
                                ExponentiatedTerm = 0.0;
                                LocalArithmeticError = true;
                            } else {
                                ThereisanOKTerm_NegativeExponential = true;
                            }
                        } else {
                            ExponentiatedTerm = 0.0;
                        }
                        Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = Math.max(MinP_k_,
                                Save_CurrentP_k_TimesExponential[IndirectClusterIndex]) * ExponentiatedTerm;
                        MalphaDenominator += Save_CurrentP_k_TimesExponential[IndirectClusterIndex];
                    } // End Second Pass Loop over Clusters for this point

                    //  Sponge Term is never rejected
                    FindUsefulCalcs.addapoint(threadIndex, NumUseful - NumSponge);
                    FindUselessCalcs.addapoint(threadIndex, (double) IndirectSize - NumUseful);
                    FindIgnoredCalcs.addapoint(threadIndex,
                            (double) ParallelClustering.runningSolution.Ncent_ThisNode - (double) IndirectSize);

                    //  Cope with ill defined arithmetic
                    if (ThereisanOKTerm_NegativeExponential) {
                        MalphaDenominator = 1.0 / MalphaDenominator;
                        if (Double.isNaN(MalphaDenominator) || Double.isInfinite(MalphaDenominator)) {
                            ThereisanOKTerm_NegativeExponential = false;
                            LocalArithmeticError = true;
                        }
                    }
                    if (!ThereisanOKTerm_NegativeExponential) {
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                            Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = 0.0;
                        }
                        if (MinimumTerm_IndirectClusterIndex < 0) {
                            DAVectorUtility.printAndThrowRuntimeException(
                                    " Error in Choosing Default Index Point " + (alpha + DAVectorUtility.PointStart_Process) + " Number of Clusters " + IndirectSize);

                        }
                        Save_CurrentP_k_TimesExponential[MinimumTerm_IndirectClusterIndex] = Math.max(MinP_k_,
                                MinimumTerm_P_k_);
                        MalphaDenominator = 1.0 / Save_CurrentP_k_TimesExponential[MinimumTerm_IndirectClusterIndex];
                    }
                    //  End Cope with ill defined arithmetic

                    //  Third Pass -- Finally set M and its change
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        double CandidateMvalue = Save_CurrentP_k_TimesExponential[IndirectClusterIndex] * MalphaDenominator;
                        if (Double.isNaN(CandidateMvalue) || Double.isInfinite(CandidateMvalue)) {
                            String message = "";
                            for (int loop = 0; loop < IndirectSize; loop++) {
                                int CreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][loop];
                                double Cvalue = 0.0;
                                double Pvalue = 0.0;
                                int ClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                                if (ClusterIndex > 0) {
                                    Cvalue = ParallelClustering.runningSolution.C_k_[ClusterIndex - 1];
                                    Pvalue = ParallelClustering.runningSolution.P_k_[ClusterIndex - 1];
                                } else {
                                    Pvalue = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedP_t[-ClusterIndex - 1];
                                }
                                message += CreatedIndex + " P " + String.format("%1$6.5f",
                                        Pvalue) + " C " + String.format(
                                        "%1$6.5f", Cvalue) + " Temp " + String.format("%1$5.4E",
                                        Save_CurrentP_k_TimesExponential[loop]) + " Term " + String.format("%1$7.6E",
                                        Save_Term_NegativeExponential[loop]) + " Old Malpha " + String.format("%1$6.5f",
                                        ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][loop]) + " * ";
                            }
                            DAVectorUtility.printAndThrowRuntimeException(
                                    "Arithmetic Error Point " + (alpha + DAVectorUtility.PointStart_Process) + " Number of Clusters " + IndirectSize + " Indirect Index " + IndirectClusterIndex + " Created Index " + ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] + " Selected Indirect Index " + MinimumTerm_IndirectClusterIndex + " ZeeX " + String.format(
                                            "%1$5.4E", MalphaDenominator) + "\n" + message
                            );

                        }

                        //  Finally set basic accumulation values
                        if (ParallelClustering.runningSolution.DiffMalpha_k_Set >= 0) {
                            AbsChangeinMalpha_k_[IndirectClusterIndex] = Math.abs(
                                    CandidateMvalue - ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                        }
                        ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = CandidateMvalue;
                        Malphax1minusMalpha[IndirectClusterIndex] = CandidateMvalue * (1.0 - CandidateMvalue);
                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        ParallelClustering.runningSolution.LegalCluster(alpha, IndirectClusterIndex);
                    }

                    //  Finally Load Contributions -- Distributed Cluster Execution
                    if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                            int LocationInThreadArray = ThreadStorePosition[IndirectClusterIndex];
                            if ((LocationInThreadArray < 0) || (LocationInThreadArray >= DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[threadIndex])) {
                                DAVectorUtility.printAndThrowRuntimeException(
                                        "Bad Thread Location " + LocationInThreadArray + " Thread " + threadIndex + " Number Clusters " + IndirectSize + " Cluster Index " + ActiveClustersperPoint[IndirectClusterIndex] + " Number Local " + ClusteringSolution.NumberLocalActiveClusters);

                            }
                            if (ParallelClustering.runningSolution.DiffMalpha_k_Set >= 0) {
                                Find3Components.addapoint(threadIndex, LocationInThreadArray, BeginFindDiffMalpha_k_,
                                        AbsChangeinMalpha_k_[IndirectClusterIndex]);
                            }
                            Find3Components.addapoint(threadIndex, LocationInThreadArray, BeginFindC_k_,
                                    ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                            Find3Components.addapoint(threadIndex, LocationInThreadArray, BeginFindFreezingMeasure_k_,
                                    Malphax1minusMalpha[IndirectClusterIndex]);
                        }
                    } else {
                        //  Finally Load Contributions -- Global Cluster Execution
                        if (ParallelClustering.runningSolution.DiffMalpha_k_Set >= 0) {
                            FindDiffMalpha_k_.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                    AbsChangeinMalpha_k_);
                        }
                        FindC_k_.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                ParallelClustering.runningSolution.M_alpha_kpointer_[alpha]);
                        FindFreezingMeasure_k_.addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                Malphax1minusMalpha);
                    }

                    //  Error is always Global over nodes
                    FindArithmeticError.addapoint(threadIndex, LocalArithmeticError);
                } // End loop over points

            }); // End parallel Thread loop
        });
        DAVectorUtility.StopSubTimer(2);

        //  Parallel Accumulations. First those using all nodes
        DAVectorUtility.StartSubTimer(3);
        FindArithmeticError.sumoverthreadsandmpi();
        VectorAnnealIterate.ArithmeticError = FindArithmeticError.TotalOr;
        FindMeanClusterCount.sumoverthreadsandmpi();
        VectorAnnealIterate.MeanClusterCount = FindMeanClusterCount.Totalmean;
        FindSinglePoints.sumoverthreadsandmpi();
        VectorAnnealIterate.PointswithClusterCount1 = FindSinglePoints.Total;

        FindUsefulCalcs.sumoverthreadsandmpi();
        FindUselessCalcs.sumoverthreadsandmpi();
        FindIgnoredCalcs.sumoverthreadsandmpi();
        VectorAnnealIterate.LocalUsefulCalcs = FindUsefulCalcs.Total;
        VectorAnnealIterate.LocalUselessCalcs = FindUselessCalcs.Total;
        VectorAnnealIterate.LocalIgnoredCalcs = FindIgnoredCalcs.Total;
        Program.SumUsefulCalcs += VectorAnnealIterate.LocalUsefulCalcs;
        Program.SumUselessCalcs += VectorAnnealIterate.LocalUselessCalcs;
        Program.SumIgnoredCalcs += VectorAnnealIterate.LocalIgnoredCalcs;

        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        { // Distributed Mode Cluster Accumulations
            Find3Components.sumoverthreadsandmpi();

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int AccumulationPosition = ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex];
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                ParallelClustering.runningSolution.C_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindC_k_];
                ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindFreezingMeasure_k_];
                if (ParallelClustering.runningSolution.DiffMalpha_k_Set >= 0) // Set DiffMalpha_k_
                {
                    ParallelClustering.runningSolution.DiffMsummed_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindDiffMalpha_k_];
                }
            }
        }
        else
        { // Non distributed Mode
            if (ParallelClustering.runningSolution.DiffMalpha_k_Set >= 0)
            { // Set DiffMalpha_k_
                FindDiffMalpha_k_.sumoverthreadsandmpi();
                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    ParallelClustering.runningSolution.DiffMsummed_k_[RealClusterIndex] = FindDiffMalpha_k_.TotalVectorSum[LocalActiveClusterIndex];
                }
            }

            FindC_k_.sumoverthreadsandmpi();
            FindFreezingMeasure_k_.sumoverthreadsandmpi();
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                ParallelClustering.runningSolution.C_k_[RealClusterIndex] = FindC_k_.TotalVectorSum[LocalActiveClusterIndex];
                ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] = FindFreezingMeasure_k_.TotalVectorSum[LocalActiveClusterIndex];

            } // End Code setting  C_k_ FreezingMeasure_k_
        }

        // Broadcast zero size information for global clusters
        boolean[] ZeroSizeClusters = new boolean[ClusteringSolution.NumberGlobalClusters];
        int CountGlobalClusters = 0;
        for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
        {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
            if ((ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 1))
            {
                continue;
            }
            boolean zerosizecluster = false;
            if (RealClusterIndex != ParallelClustering.runningSolution.SpongeCluster)
            {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }
            ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
            ++CountGlobalClusters;
        }
        DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
        // Note - MPI Call - Broadcast - boolean []
        if (DAVectorUtility.MPI_Size > 1){
            //            DAVectorUtility.MPI_communicator.<Boolean>Broadcast(tempRef_ZeroSizeClusters, 0);
            DAVectorUtility.mpiOps.broadcast(ZeroSizeClusters, 0);
        }
        DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        CountGlobalClusters = 0;
        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        { // Distributed Mode Cluster Accumulations

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                boolean largeenoughcluster = (ParallelClustering.runningSolution.C_k_[RealClusterIndex] > Program.CountforCluster_C_ktobezero);
                if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)
                {
                    largeenoughcluster = !ZeroSizeClusters[CountGlobalClusters];
                    ++CountGlobalClusters;
                }
                if (largeenoughcluster)
                {
                    ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] = ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] / ParallelClustering.runningSolution.C_k_[RealClusterIndex];
                }
            }
        }
        else
        { // Non distributed Mode
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                boolean zerosizecluster = false;
                if ((ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 0) || (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 1))
                {
                    zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                    ++CountGlobalClusters;
                }
                if (!zerosizecluster)
                {
                    ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] = ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] / ParallelClustering.runningSolution.C_k_[RealClusterIndex];
                }

            } // End Code normalizing FreezingMeasure_k_
        }

        DAVectorUtility.StopSubTimer(3);

        DAVectorUtility.StartSubTimer(9);
        //  Indicate status of setting Malpha differences
        ++ParallelClustering.runningSolution.DiffMalpha_k_Set;

        //  Set new value of P_k_ including case where value for Sponge Cluster is special
        double wgt = 1.0 / DAVectorUtility.PointCount_Global;
        double C_k_Sum01 = 0.0;
        double C_k_Sum2 = 0.0;
        for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
        {
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0)
            {
                continue;
            }
            ParallelClustering.runningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.runningSolution.C_k_[RealClusterIndex];
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
            {
                C_k_Sum2 += ParallelClustering.runningSolution.C_k_[RealClusterIndex];
            }
            else
            {
                C_k_Sum01 += ParallelClustering.runningSolution.C_k_[RealClusterIndex];
            }
        }
        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        {
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - double - sum
                //                C_k_Sum2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(C_k_Sum2, Operation<Double>.Add);
                C_k_Sum2 = DAVectorUtility.mpiOps.allReduce(C_k_Sum2, MPI.SUM);
            }
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
        }
        VectorAnnealIterate.C_k_Sum = C_k_Sum01 + C_k_Sum2;

        // Correctness check
        double CDiff = C_k_Sum - (double)DAVectorUtility.PointCount_Global;
        if (Math.abs(CDiff) > 0.5)
        {
            DAVectorUtility.printAndThrowRuntimeException("C_k_Sum bad " + String.format("%1$3.2f", VectorAnnealIterate.C_k_Sum) + " Should be " + DAVectorUtility.PointCount_Global);

        }

        if ((ParallelClustering.runningSolution.SpongeCluster != -1) && (Program.SpongePoption == 1))
        {
            wgt = 1.0 - ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster];
            ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster] = Program.SpongePWeight / ParallelClustering.runningSolution.Ncent_Global;
            wgt = (1.0 - ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster]) / wgt;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
                {
                    continue;
                }
                if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 0)
                {
                    ParallelClustering.runningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.runningSolution.P_k_[RealClusterIndex];
                }
            }
        }
        DAVectorUtility.StopSubTimer(9);
        DAVectorUtility.StartSubTimer(10);

        //  Set Y_k_
        GlobalReductions.FindIndirectVectorDoubleSum[] FindY_k_;
        FindY_k_ = new GlobalReductions.FindIndirectVectorDoubleSum[Program.ParameterVectorDimension];
        final DistributedReductions.FindIndirectMultiVectorDoubleSum FindY_k_Component;
        final int BeginFindY_k_;
        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        {
            FindY_k_Component = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
            BeginFindY_k_ = FindY_k_Component.AddComponents(Program.ParameterVectorDimension);
            FindY_k_Component.NodeInitialize();
        }
        else
        {
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                FindY_k_[VectorIndex] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            }
            FindY_k_Component = null;
            BeginFindY_k_ = -1;
        }

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                    FindY_k_Component.ThreadInitialize(threadIndex);
                } else {
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                        FindY_k_[VectorIndex].startthread(threadIndex);
                    }
                }
                int ArraySize = Math.min(ClusteringSolution.NumberAvailableActiveClusters,
                        ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[][] CenterPositionsperCluster = new double[Program.ParameterVectorDimension][];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                    CenterPositionsperCluster[VectorIndex] = new double[ArraySize];
                }
                double[] TerminY_k_ = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                    int ThreadStorePosition = -1;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        int RealClusterIndex = -2;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        Box<Integer> tempRef_RealClusterIndex2 = new Box<>(RealClusterIndex);
                        Box<Integer> tempRef_ActiveClusterIndex2 = new Box<>(ActiveClusterIndex);
                        Box<Integer> tempRef_RemoteIndex2 = new Box<>(RemoteIndex);
                        ClusterPointersforaPoint(alpha, IndirectClusterIndex, tempRef_RealClusterIndex2,
                                tempRef_ActiveClusterIndex2, tempRef_RemoteIndex2);
                        RealClusterIndex = tempRef_RealClusterIndex2.content;
                        ActiveClusterIndex = tempRef_ActiveClusterIndex2.content;
                        RemoteIndex = tempRef_RemoteIndex2.content;
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        if (RemoteIndex < 0) {
                            if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][threadIndex];
                            }
                        } else {
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][threadIndex];
                        }
                        if ((ParallelClustering.runningSolution.SpongeCluster >= 0) && (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)) {
                            continue;
                        }
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                            TerminY_k_[VectorIndex] = Program.PointPosition[alpha][VectorIndex] * ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            CenterPositionsperCluster[VectorIndex][IndirectClusterIndex] = TerminY_k_[VectorIndex];
                        }
                        if (ParallelClustering.runningSolution.DistributedExecutionMode) {
                            FindY_k_Component.addapoint(threadIndex, ThreadStorePosition, BeginFindY_k_,
                                    Program.ParameterVectorDimension, TerminY_k_);
                            int CreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        }
                    }
                    if (!ParallelClustering.runningSolution.DistributedExecutionMode) {
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                            FindY_k_[VectorIndex].addapoint(threadIndex, IndirectSize, ActiveClustersperPoint,
                                    CenterPositionsperCluster[VectorIndex]);
                        }
                    }
                }

            }); // End Parallel Section calculating center positions
        });

        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        {
            FindY_k_Component.sumoverthreadsandmpi();
        }
        else
        {
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                FindY_k_[VectorIndex].sumoverthreadsandmpi();
            }
        }
        // End section calculating centers

        //  Normalize Centers correctly
        CountGlobalClusters = 0;
        for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
        {
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 0)
            {
                continue;
            }
            boolean zerosizecluster = false;
            if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)
            {
                zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                ++CountGlobalClusters;
            }
            else if (ParallelClustering.runningSolution.DistributedExecutionMode)
            {
                zerosizecluster = ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
            }

            if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
            {
                continue;
            }

            if (!zerosizecluster)
            { // If zero size cluster, leave Y unchanged
                if (ParallelClustering.runningSolution.DistributedExecutionMode)
                {
                    int NodeAccumulationIndex = ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex];
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        if (ParallelClustering.runningSolution.YPreviousSet > -1)
                        {
                            ParallelClustering.runningSolution.YPrevious_k_i_[RealClusterIndex][VectorIndex] = ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                            ParallelClustering.runningSolution.YPreviousActuallySet[RealClusterIndex] = true;
                        }
                        ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][VectorIndex] = FindY_k_Component.TotalVectorSum[NodeAccumulationIndex][VectorIndex] / ParallelClustering.runningSolution.C_k_[RealClusterIndex];
                    }
                }
                else
                {
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        if (ParallelClustering.runningSolution.YPreviousSet > -1)
                        {
                            ParallelClustering.runningSolution.YPrevious_k_i_[RealClusterIndex][VectorIndex] = ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                            ParallelClustering.runningSolution.YPreviousActuallySet[RealClusterIndex] = true;
                        }
                        ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][VectorIndex] = FindY_k_[VectorIndex].TotalVectorSum[ActiveClusterIndex] / ParallelClustering.runningSolution.C_k_[RealClusterIndex];
                    }
                }

                if (Program.SigmaMethod > 1)
                {
                    Box<double[]> tempRef_Object = new Box<>(ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex]);
                    Program.CalculateSigma(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], tempRef_Object);
                    ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex] = tempRef_Object.content;
                }
            }
        } // End section setting Y_k_

        //  Section calculating change in Y
        VectorAnnealIterate.AverageY_k_SquaredChange = 0.0;
        if (ParallelClustering.runningSolution.YPreviousSet > -1)
        {
            double SumY_k_SquaredChange01 = 0.0;
            double SumY_k_SquaredChange2 = 0.0;
            double NumGlobalSet = 0.0;
            double NumDistributedSet = 0.0;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] <= 0)
                {
                    continue;
                }
                if (!ParallelClustering.runningSolution.YPreviousActuallySet[RealClusterIndex])
                {
                    continue;
                }
                double SquaredDifference_k_ = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], ParallelClustering.runningSolution.YPrevious_k_i_[RealClusterIndex], ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex]);
                if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
                {
                    SumY_k_SquaredChange2 += SquaredDifference_k_;
                    NumDistributedSet += 1.0;
                }
                else
                {
                    SumY_k_SquaredChange01 += SquaredDifference_k_;
                    NumGlobalSet += 1.0;
                }
            }
            if (ParallelClustering.runningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                if (DAVectorUtility.MPI_Size > 1){
                    // Note - MPI Call - Allreduce - double - sum
                    //                    SumY_k_SquaredChange2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(SumY_k_SquaredChange2, Operation<Double>.Add);
                    SumY_k_SquaredChange2 = DAVectorUtility.mpiOps.allReduce(SumY_k_SquaredChange2, MPI.SUM);
                    // Note - MPI Call - Allreduce - double - sum
                    //                    NumDistributedSet = DAVectorUtility.MPI_communicator.<Double>Allreduce(NumDistributedSet, Operation<Double>.Add);
                    NumDistributedSet = DAVectorUtility.mpiOps.allReduce(NumDistributedSet, MPI.SUM);
                }
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            if ((NumDistributedSet + NumGlobalSet) > 0.5)
            {
                VectorAnnealIterate.AverageY_k_SquaredChange = (SumY_k_SquaredChange01 + SumY_k_SquaredChange2) / (NumDistributedSet + NumGlobalSet);
            }
        } // End computation of change of vector difference squared
        ++ParallelClustering.runningSolution.YPreviousSet;

        ParallelClustering.runningSolution.CorrelationsSet = false;
        DAVectorUtility.StopSubTimer(10);
        if (ParallelClustering.runningSolution.DistributedExecutionMode)
        {
            DistributedClusteringSolution.MinorSynchronizationTransportDistributedClusterCenters();
        }

    } // End DAVectorEMIterate()

	// The change of sum(points) Delta(Malpha)/# Points should be less than argument ChangeLimit summed over points/clusters
	//  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
	//  Note divided by number of points here not in calculation
	// Note hitting Program.ConvergenceLoopLimit limit is NOT fatal. Just continues
	public static int convergenceTest(double ChangeLimit, Box<String> ReasonforConvergence) throws MPIException {
		ReasonforConvergence.content = "Not Converged";
		if (VectorAnnealIterate.ArithmeticError)
		{
			VectorAnnealIterate.OnLastleg = true;
			ReasonforConvergence.content = "Arithmetic Error";
			return 2;
		}
		VectorAnnealIterate.HitConvergenceLoopLimit = false;
		if (ParallelClustering.runningSolution.DiffMalpha_k_Set < 1)
		{
			return 0;
		}
		double MalphaDiffAvg01 = 0.0;
		double MalphaDiffAvg2 = 0.0;
		for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
			{
				MalphaDiffAvg2 += ParallelClustering.runningSolution.DiffMsummed_k_[RealClusterIndex];
			}
			else
			{
				MalphaDiffAvg01 += ParallelClustering.runningSolution.DiffMsummed_k_[RealClusterIndex];
			}
		}
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - double - sum
//                MalphaDiffAvg2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(MalphaDiffAvg2, Operation<Double>.Add);
                MalphaDiffAvg2 = DAVectorUtility.mpiOps.allReduce(MalphaDiffAvg2, MPI.SUM);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		}
		VectorAnnealIterate.MalphaDiffAvg = (MalphaDiffAvg01 + MalphaDiffAvg2) / DAVectorUtility.PointCount_Global;

		++Program.NumberMdiffSums; // Accumulate diagnostic statistics
		Program.MdiffSum += VectorAnnealIterate.MalphaDiffAvg;

		boolean notconverged1 = VectorAnnealIterate.MalphaDiffAvg > ChangeLimit;
		if (VectorAnnealIterate.EMIterationStepCount1 < 0)
		{
			VectorAnnealIterate.EMIterationStepCount1 = -2;
		}
		if ((VectorAnnealIterate.EMIterationStepCount1 < 0) && (!notconverged1))
		{
			VectorAnnealIterate.EMIterationStepCount1 = VectorAnnealIterate.EMIterationStepCount;
		}

		boolean notconverged2 = false;
		double YchangeTest = 0.0;
		ReasonforConvergence.content = "M cnvgd";
		if (VectorAnnealIterate.FinalLoop)
		{
			YchangeTest = ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth;
			notconverged2 = VectorAnnealIterate.AverageY_k_SquaredChange > Program.YChangeSquared * YchangeTest;
			if (VectorAnnealIterate.EMIterationStepCount2 < 0)
			{
				VectorAnnealIterate.EMIterationStepCount2 = -2;
			}
			if ((VectorAnnealIterate.EMIterationStepCount2 < 0) && (!notconverged2))
			{
				VectorAnnealIterate.EMIterationStepCount2 = VectorAnnealIterate.EMIterationStepCount;
			}
			ReasonforConvergence.content = "Y and M cnvgd";
		}
		boolean notconverged = notconverged1 || notconverged2;
		notconverged = DAVectorUtility.SynchronizeMPIBoolean(notconverged);
		if (notconverged)
		{
			if (notconverged1)
			{
				ReasonforConvergence.content = "Converging Malpha " + String.format("%1$5.4E", VectorAnnealIterate.MalphaDiffAvg) + " " + String.format("%1$5.4E", ChangeLimit) + " ";
			}
			if (notconverged2)
			{
				ReasonforConvergence.content += "Converging Center Change " + String.format("%1$5.4E", VectorAnnealIterate.AverageY_k_SquaredChange) + " " + String.format("%1$5.4E", YchangeTest);
			}
			if (VectorAnnealIterate.EMIterationStepCount > Program.ConvergenceLoopLimit)
			{
				if (notconverged)
				{
					DAVectorUtility.SALSAPrint(1, " Too many Iterations " + VectorAnnealIterate.EMIterationStepCount + " Warning Message " + ReasonforConvergence.content);
				}
				VectorAnnealIterate.HitConvergenceLoopLimit = true;
			}
			else
			{
				return 0;
			}
		}

		//  "Converged" either due StepCount Limit or Malpha (plus Y) change
		if (VectorAnnealIterate.OnLastleg || VectorAnnealIterate.FinalLoop)
		{
			return 2;
		}

		boolean toomanyclusters = (ParallelClustering.runningSolution.Ncent_ThisNode >= VectorAnnealIterate.ActualMaxNcent) || (ParallelClustering.runningSolution.Ncent_Global >= Program.maxNcentTOTAL);
		DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
        if (DAVectorUtility.MPI_Size > 1){
            // Note - MPI Call - Allreduce - boolean - logical OR
//            toomanyclusters = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(toomanyclusters, Operation<Boolean>.LogicalOr);
            toomanyclusters = DAVectorUtility.mpiOps.allReduce(toomanyclusters,MPI.LOR);
        }

		DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

		boolean toolowtemperature = ParallelClustering.runningSolution.Temperature < VectorAnnealIterate.Tmin;
		toolowtemperature = DAVectorUtility.SynchronizeMPIBoolean(toolowtemperature);
		if (toomanyclusters || toolowtemperature || (VectorAnnealIterate.SplitFailures > 1))
		{
			String reason = "";
			if (toomanyclusters)
			{
				reason = "Too many Clusters";
			}
			if (toolowtemperature)
			{
				reason += " Temperature Limit Reached";
			}
			if (VectorAnnealIterate.SplitFailures > 1)
			{
				reason += " Too Many Split Failures";
			}
			DAVectorUtility.SALSAPrint(1, " On last Leg due to " + reason);
			VectorAnnealIterate.OnLastleg = true;
			ReasonforConvergence.content = " On last Leg due to " + reason;
			return 2;
		}
		return 1;

	} // End convergenceTest

	public final void SaveCurrentTask()
	{ // Save current task that should be split but we need to converge first

		ClusteringSolution.CopySolution(ParallelClustering.runningSolution, ParallelClustering.savedSolution);

    } // End saving current task that should be split but we need to converge first

	public final void RestorePreviousTask() throws MPIException { // Restore previous task that should be split but we needed to converge current cluster configuration first

		ClusteringSolution.CopySolution(ParallelClustering.savedSolution, ParallelClustering.runningSolution);
		DistributedClusteringSolution.ManageMajorSynchronization(true);

    } // End Restore previous task that should be split but we needed to converge current cluster configuration first

	// Decide if to split based on negative eigenvalue of second derivative matrix
	// MinimumEigenvalue is MINIMUM eigenvalue
	// ClustertoSplit is cluster number of this
	// In distributed mode ONLY distributed clusters split and this is done asynchronously  in each node 
	//  In non distributed mode, clusters are split synchronously
	public final boolean shouldweSplit() throws MPIException {
		//	Calculate Cluster with minimum eigenvalue -- find cluster number and eigenvalue (which could be positive)
		VectorAnnealIterate.Numberthatcanbesplit = 0;
		int LimitonSplits = Math.min(Program.MaxNumberSplitClusters, VectorAnnealIterate.ActualMaxNcent - ParallelClustering.runningSolution.Ncent_ThisNode);
		ParallelClustering.runningSolution.ClustertoSplit = -1;
		boolean eigenvaluesexist = false;

		VectorClass vc = new VectorClass();
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			if (ParallelClustering.runningSolution.DistributedExecutionMode)
			{
				if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] != 2)
				{
					ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 0;
					continue;
				}
			}
			else
			{
				if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] != 1)
				{
					ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 0;
					continue;
				}
			}
			ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 1;
			if (ParallelClustering.runningSolution.SplitPriority_k_[RealClusterIndex] == 0)
			{
				ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 0;
			}

			if ((ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] <= Program.MinimumScaledWidthsquaredtosplit) || (ParallelClustering.runningSolution.C_k_[RealClusterIndex] <= Program.ToosmalltoSplit) || (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster))
			{
				ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 0;
			}
		}
		if (!ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int[]
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_Splittable_k_, 0);
                DAVectorUtility.mpiOps.broadcast(ParallelClustering.runningSolution.Splittable_k_, 0);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
		}

		//  Set up eigenvalues
		if (!Program.CalculateEigenvaluesfromMatrix)
		{
			vc.SetAllEigenvaluesIteratively(ParallelClustering.runningSolution);
		}

		// Only implement Continuous Clustering 
		//  Unlike pairwise call geteigenvalues separately for each cluster
		//   Set Correlation and Change Direction
		ParallelClustering.runningSolution.SetClusterCorrelations();
		int[] eigenvaluenegative = new int[ParallelClustering.runningSolution.Ncent_ThisNode];
		int[] eigenvaluestatus = new int[ParallelClustering.runningSolution.Ncent_ThisNode];
		for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.runningSolution.Ncent_ThisNode; ClusterToRefine++)
		{
			eigenvaluenegative[ClusterToRefine] = -1;
			eigenvaluestatus[ClusterToRefine] = -1;
			if (ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] == 0)
			{
				ParallelClustering.runningSolution.Eigenvalue_k[ClusterToRefine] = 0.0;
				continue;
			}
			if (Program.CalculateEigenvaluesfromMatrix)
			{
				double[][] secondderivmatrix = new double[Program.ParameterVectorDimension][Program.ParameterVectorDimension];
				for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
				{
					for (int VectorIndex2 = 0; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
					{
						secondderivmatrix[VectorIndex1][VectorIndex2] = -ParallelClustering.runningSolution.Correlation_k_i_j[ClusterToRefine][VectorIndex1][VectorIndex2];
						if (VectorIndex1 == VectorIndex2)
						{
							secondderivmatrix[VectorIndex1][VectorIndex2] += ParallelClustering.runningSolution.C_k_[ClusterToRefine] * 0.5;
						}
					}
				}
				vc.getEigenvaluefromMatrix(secondderivmatrix);
			}
			else
			{
				vc.getEigenvaluefromIteration(ClusterToRefine);
			}
			eigenvaluestatus[ClusterToRefine] = vc.EigenStatus;
			if (vc.EigenStatus > 0)
			{
				eigenvaluesexist = true;
				if (Program.CalculateEigenvaluesfromMatrix)
				{ // In iterative case, data already stored
                    System.arraycopy(vc.Eigenvector, 0,
                            ParallelClustering.runningSolution.Eigenvector_k_i[ClusterToRefine], 0,
                            Program.ParameterVectorDimension);
					ParallelClustering.runningSolution.Eigenvalue_k[ClusterToRefine] = vc.Eigenvalue;
				}

				boolean EigenvalueNegative = vc.Eigenvalue < 0.0;
				if (EigenvalueNegative)
				{
					eigenvaluenegative[ClusterToRefine] = 0;
				}
				else
				{
					eigenvaluenegative[ClusterToRefine] = 1;
				}
			}
			if (ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine] >= 1)
			{
				++ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine];
				if (ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine] >= 9)
				{
					ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine] = -1;
				}

				if (ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine] < 6)
				{
					ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] = 0;
                }

			}
		}
		if (!ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int[]
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_eigenvaluenegative, 0);
                DAVectorUtility.mpiOps.broadcast(eigenvaluenegative, 0);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
		}

		for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.runningSolution.Ncent_ThisNode; ClusterToRefine++)
		{
			if (ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] == 0)
			{
				continue;
			}
			if (eigenvaluenegative[ClusterToRefine] == -1)
			{
				continue;
			}
			double CurrentClusterMinimumEigenvalue = ParallelClustering.runningSolution.Eigenvalue_k[ClusterToRefine];
			boolean EigenvalueNegative = false;
			if (eigenvaluenegative[ClusterToRefine] == 0)
			{
				EigenvalueNegative = true;
			}
			ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] = 2;
			if (EigenvalueNegative)
			{
				ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] = 3;
			}
		}
		if ((DAVectorUtility.MPI_Rank == 0) || ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.runningSolution.Ncent_ThisNode; ClusterToRefine++)
			{
				if (ParallelClustering.runningSolution.Splittable_k_[ClusterToRefine] == 0)
				{
					continue;
				}
				if (eigenvaluenegative[ClusterToRefine] == -1)
				{
					continue;
				}
				int CurrentClusterPriority = ParallelClustering.runningSolution.SplitPriority_k_[ClusterToRefine];
				if (CurrentClusterPriority > 0)
				{
					CurrentClusterPriority += -10;
				}
				double CurrentClusterMinimumEigenvalue = ParallelClustering.runningSolution.Eigenvalue_k[ClusterToRefine];

				if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
				{
					String correlationmessage = "";
					if (Program.CalculateCorrelationMatrix && Program.ParameterVectorDimension == 2)
					{
						correlationmessage = " Correls " + String.format("%1$5.4E", ParallelClustering.runningSolution.Correlation_k_i_j[ClusterToRefine][0][0]) + " " + String.format("%1$5.4E", ParallelClustering.runningSolution.Correlation_k_i_j[ClusterToRefine][1][1]) + " " + String.format("%1$5.4E", ParallelClustering.runningSolution.Correlation_k_i_j[ClusterToRefine][0][1]);
					}
					String eigenvectormessageY = "";
					String eigenvectormessageDeltaY = "";
					if (Program.Printeigenvectors)
					{
						eigenvectormessageY = " Y ";
						eigenvectormessageDeltaY = "Delta Y ";
						for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
						{
							eigenvectormessageY += String.format("%1$5.4E", ParallelClustering.runningSolution.Y_k_i_[ClusterToRefine][VectorIndex]) + " ";
							eigenvectormessageDeltaY += String.format("%1$5.4E", ParallelClustering.runningSolution.DC_k_DY_k_i_[ClusterToRefine][VectorIndex]) + " ";
						}
					}

					if ((ClusterToRefine < Program.ClusterPrintNumber) || (ClusterToRefine >= (ParallelClustering.runningSolution.Ncent_ThisNode - Program.ClusterPrintNumber)))
					{
						DAVectorUtility.SALSAPrint(1, "       Cluster " + ClusterToRefine + "(" + ParallelClustering.runningSolution.LocalCreatedIndex[ClusterToRefine] + ")" + " Status " + eigenvaluestatus[ClusterToRefine] + " Priority " + CurrentClusterPriority + " Eigen " + String.format("%1$5.4E", CurrentClusterMinimumEigenvalue) + " Size " + String.format("%1$2.1f", ParallelClustering.runningSolution.C_k_[ClusterToRefine]) + eigenvectormessageY + "  " + eigenvectormessageDeltaY + correlationmessage);
					}
				}

				boolean EigenvalueNegative = false;
				if (eigenvaluenegative[ClusterToRefine] == 0)
				{
					EigenvalueNegative = true;
				}
				if (EigenvalueNegative)
				{ // Candidate for Split List

					double test1 = 2.0 * ((double) DAVectorUtility.PointCount_Global) / ((double) ParallelClustering.runningSolution.Ncent_ThisNode);
					if (ParallelClustering.runningSolution.C_k_[ClusterToRefine] >= test1)
					{
						CurrentClusterPriority += 20;
					}
					if (ParallelClustering.runningSolution.C_k_[ClusterToRefine] >= (3.0 * test1))
					{
						CurrentClusterPriority += 20;
					}

					if (VectorAnnealIterate.Numberthatcanbesplit == 0) // Initialize Split List
					{
						VectorAnnealIterate.Numberthatcanbesplit = 1;
						VectorAnnealIterate.EigsofClusterstoSplit[0] = CurrentClusterMinimumEigenvalue;
						VectorAnnealIterate.ListofClusterstoSplit[0] = ClusterToRefine;
						VectorAnnealIterate.PrioritiesofClusterstoSplit[0] = CurrentClusterPriority;
					}
					else // Add to Split List
					{
						int position = VectorAnnealIterate.Numberthatcanbesplit;
						for (int positionloop = 0; positionloop < VectorAnnealIterate.Numberthatcanbesplit; positionloop++)
						{
							if (CurrentClusterPriority < VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop])
							{
								continue;
							}
							if (CurrentClusterPriority == VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop])
							{
								boolean eigtest = (CurrentClusterMinimumEigenvalue >= VectorAnnealIterate.EigsofClusterstoSplit[positionloop]);
								if (eigtest)
								{
									continue;
								}
							}
							position = positionloop;
							break;
						}
						if (position >= LimitonSplits)
						{
							continue;
						}
						for (int positionloop = VectorAnnealIterate.Numberthatcanbesplit - 1; positionloop >= position; positionloop--)
						{
							if (positionloop == (LimitonSplits - 1))
							{
								continue;
							}
							VectorAnnealIterate.EigsofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.EigsofClusterstoSplit[positionloop];
							VectorAnnealIterate.ListofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.ListofClusterstoSplit[positionloop];
							VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop];
						}
						VectorAnnealIterate.Numberthatcanbesplit = Math.min(VectorAnnealIterate.Numberthatcanbesplit + 1, LimitonSplits);
						VectorAnnealIterate.EigsofClusterstoSplit[position] = CurrentClusterMinimumEigenvalue;
						VectorAnnealIterate.ListofClusterstoSplit[position] = ClusterToRefine;
						VectorAnnealIterate.PrioritiesofClusterstoSplit[position] = CurrentClusterPriority;
					}
				}
			}
		}
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - boolean - logical OR
//                eigenvaluesexist = DAVectorUtility.MPI_communicator.<Boolean>Allreduce(eigenvaluesexist, Operation<Boolean>.LogicalOr);
                eigenvaluesexist = DAVectorUtility.mpiOps.allReduce(eigenvaluesexist, MPI.LOR);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		}
		else
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_Numberthatcanbesplit, 0);
                VectorAnnealIterate.Numberthatcanbesplit = DAVectorUtility.mpiOps.broadcast(VectorAnnealIterate.Numberthatcanbesplit,0);
            }
            // Note - MPI Call - Broadcast - double []
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Double>Broadcast(tempRef_EigsofClusterstoSplit, 0);
                DAVectorUtility.mpiOps.broadcast(VectorAnnealIterate.EigsofClusterstoSplit, 0);
            }
            // Note - MPI Call - Broadcast - int []
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_ListofClusterstoSplit, 0);
                DAVectorUtility.mpiOps.broadcast(VectorAnnealIterate.ListofClusterstoSplit, 0);
            }
            // Note - MPI Call - Broadcast - int []
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_PrioritiesofClusterstoSplit, 0);
                DAVectorUtility.mpiOps.broadcast(VectorAnnealIterate.PrioritiesofClusterstoSplit, 0);

            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
		}
		if (!eigenvaluesexist)
		{
			VectorAnnealIterate.SplitFailures++;
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_SplitFailures, 0);
                VectorAnnealIterate.SplitFailures = DAVectorUtility.mpiOps.broadcast(VectorAnnealIterate.SplitFailures,0);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
			return false;
		}
		if (VectorAnnealIterate.Numberthatcanbesplit == 0)
		{
			return false;
		}
		ParallelClustering.runningSolution.ClustertoSplit = VectorAnnealIterate.ListofClusterstoSplit[0];
		return true;

	} // End shouldweSplit

	//  Do a split on the identified cluster
	//  New clusters stored in last cluster and old position
	//  Cluster count is incremented by one
	//  ClustertoSplit is cluster number to split
	public final void dothesplit() throws MPIException { // This just gets Malpha_k_ P_k_ and Y_k_

		if (ParallelClustering.runningSolution.ClustertoSplit < 0)
		{
			return;
		}
		int NewCenterIndex = ParallelClustering.runningSolution.Ncent_ThisNode;
		int OldCenterIndex = ParallelClustering.runningSolution.ClustertoSplit;
		++ClusteringSolution.ClustersSplit;

		//  Calculate Perturbed center positions
		// First Estimate Scale as minumum of one that produces .05 change in population of new centers and one that shifts 0.1 times cluster width
		double ambitiousscale = 0.0;
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			ambitiousscale += ParallelClustering.runningSolution.Eigenvector_k_i[OldCenterIndex][VectorIndex] * ParallelClustering.runningSolution.DC_k_DY_k_i_[OldCenterIndex][VectorIndex];
		}
		ambitiousscale = Math.abs(0.025 * ParallelClustering.runningSolution.C_k_[OldCenterIndex] / ambitiousscale);
		double actualscale = Math.min(ambitiousscale, 0.2 * Math.sqrt(ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[OldCenterIndex]));
		actualscale = Math.max(actualscale, 0.05 * Math.sqrt(ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[OldCenterIndex]));

		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			double tmp1 = ParallelClustering.runningSolution.Y_k_i_[OldCenterIndex][VectorIndex];
			double tmp2 = ParallelClustering.runningSolution.Eigenvector_k_i[OldCenterIndex][VectorIndex] * actualscale * Math.sqrt(ParallelClustering.runningSolution.Sigma_k_i_[OldCenterIndex][VectorIndex]);
			ParallelClustering.runningSolution.Y_k_i_[OldCenterIndex][VectorIndex] = tmp1 + tmp2;
			ParallelClustering.runningSolution.Y_k_i_[NewCenterIndex][VectorIndex] = tmp1 - tmp2;
		}
		ParallelClustering.runningSolution.YPreviousActuallySet[OldCenterIndex] = false;
		ParallelClustering.runningSolution.YPreviousActuallySet[NewCenterIndex] = false;
		Box<double[]> tempRef_Object = new Box<>(ParallelClustering.runningSolution.Sigma_k_i_[NewCenterIndex]);
		Program.CalculateSigma(ParallelClustering.runningSolution.Y_k_i_[NewCenterIndex], tempRef_Object);
		ParallelClustering.runningSolution.Sigma_k_i_[NewCenterIndex] = tempRef_Object.content;
		ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;

		ParallelClustering.runningSolution.Splittable_k_[OldCenterIndex] = -1;
		ParallelClustering.runningSolution.Splittable_k_[NewCenterIndex] = -1;
		ParallelClustering.runningSolution.SplitPriority_k_[OldCenterIndex] = 2;
		ParallelClustering.runningSolution.SplitPriority_k_[NewCenterIndex] = 2;
		ParallelClustering.runningSolution.LocalStatus[NewCenterIndex] = ParallelClustering.runningSolution.LocalStatus[OldCenterIndex];
		ParallelClustering.runningSolution.LocalSplitCreatedIndex[NewCenterIndex] = 0;
		final int CreatedIndex_child = ClusteringSolution.SetCreatedIndex(NewCenterIndex);
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			ParallelClustering.runningSolution.LocalSplitCreatedIndex[OldCenterIndex] = -1 - CreatedIndex_child;
			ParallelClustering.runningSolution.LocalSplitCreatedIndex[NewCenterIndex] = +1 + ParallelClustering.runningSolution.LocalCreatedIndex[OldCenterIndex];
		}
		ParallelClustering.runningSolution.P_k_[OldCenterIndex] = 0.5 * ParallelClustering.runningSolution.P_k_[OldCenterIndex];
		ParallelClustering.runningSolution.P_k_[NewCenterIndex] = ParallelClustering.runningSolution.P_k_[OldCenterIndex];
		ParallelClustering.runningSolution.C_k_[OldCenterIndex] = 0.5 * ParallelClustering.runningSolution.C_k_[OldCenterIndex];
		ParallelClustering.runningSolution.C_k_[NewCenterIndex] = ParallelClustering.runningSolution.C_k_[OldCenterIndex];

		// Increase number of Clusters
		ParallelClustering.runningSolution.Ncent_ThisNode++;

		//  Exit if Distributed
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			return;
		}

		ParallelClustering.runningSolution.SetActiveClusters();
		if (Program.UseTriangleInequality_DA > 0)
		{
			DATriangleInequality.SplitCenter(OldCenterIndex, ParallelClustering.runningSolution.LocalCreatedIndex[OldCenterIndex], NewCenterIndex, CreatedIndex_child);
		}

		final int OldActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[OldCenterIndex];

		//  Parallel Section Splitting Malpha_k_ for non distributed mode
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1,
                    (threadIndex) ->  //  End Parallel Section Splitting Cluster
                    {
                        int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                        int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                        for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                            int IndirectClusterIndex = ParallelClustering.runningSolution.MapClusterToIndirect(alpha,
                                    OldActiveClusterIndex);
                            if (IndirectClusterIndex < 0) {
                                continue;
                            }
                            double newvalueofM = ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] * 0.5;
                            ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = newvalueofM;
                            int NewIndirectClusterIndex = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                            ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][NewIndirectClusterIndex] = CreatedIndex_child;
                            ++ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                            ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][NewIndirectClusterIndex] = newvalueofM;
                        }

                    });
        });

    } // End dothesplit

    public static boolean CheckNumberofClustersTooBig() throws MPIException {
        boolean changed = false;
        GlobalReductions.FindIntSum NumberChanged = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);

        //  Check Number of Clusters for each point. This could be stricter
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    if (ParallelClustering.runningSolution.NumClusters_alpha_[alpha] <= ClusteringSolution.TargetClustersperPoint) {
                        continue;
                    }
                    ParallelClustering.runningSolution.SetClustersforaPoint(alpha);
                    NumberChanged.addapoint(threadIndex, 1);
                }

            }); //  End Section Setting Cluster numbers for points
        });
        NumberChanged.sumoverthreadsandmpi();

        if (NumberChanged.TotalInt != 0)
        {
            if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
            {
                DAVectorUtility.SALSAPrint(1, " Points Changed " + NumberChanged.TotalInt);
            }
            changed = true;
        }
        return changed;

    } // End CheckNumberofClustersTooBig()

	public static void CompleteCleanUp()
	{
		DAVectorUtility.SALSAPrint(0, "Complete Clean Up T " + ParallelClustering.runningSolution.Temperature);
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    ParallelClustering.runningSolution.SetClustersforaPoint(alpha);
                }

            }); //  End Section Setting Cluster numbers for points
        });

        VectorAnnealIterate.HammyNotSet = true;
		ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;

	} // End CompleteCleanUp()

	//	Find initial Temperature Cluster Center and other initial parameters
	public static void InitializeSolution(ClusteringSolution StartSolution) throws MPIException {

		//  Deterministic Annealing
		double[] Initial_Y = new double[Program.ParameterVectorDimension];
		int SpongePosition = StartSolution.SpongeCluster;
		final double InitialM;
		StartSolution.P_k_[0] = 1.0;
		final int FirstRealCluster;
		if (SpongePosition >= 0)
		{
			InitialM = 0.5;
			StartSolution.P_k_[0] = 0.5;
			StartSolution.P_k_[1] = 0.5;
			FirstRealCluster = 1;
		} else {
            FirstRealCluster = 0;
            InitialM = 1.0;
        }
		StartSolution.SplitPriority_k_[FirstRealCluster] = -1;
		StartSolution.Splittable_k_[FirstRealCluster] = 0;
		StartSolution.LocalSplitCreatedIndex[FirstRealCluster] = 0;
		StartSolution.LocalStatus[FirstRealCluster] = 1;
		int CreatedIndex = ClusteringSolution.SetCreatedIndex(FirstRealCluster);

		GlobalReductions.FindArrayMean SystemCoG = new GlobalReductions.FindArrayMean(DAVectorUtility.ThreadCount, Program.ParameterVectorDimension);
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			if (SpongePosition >= 0)
			{
				StartSolution.Y_k_i_[SpongePosition][VectorIndex] = 0.0; // Set "center" of sponge to Origin
			}
		}

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    SystemCoG.addapoint(threadIndex, Program.PointPosition[alpha]);
                    StartSolution.Map_alpha_PointertoCreatedIndex[alpha][FirstRealCluster] = CreatedIndex;
                    StartSolution.M_alpha_kpointer_[alpha][FirstRealCluster] = InitialM;
                    if (SpongePosition >= 0) {
                        StartSolution.M_alpha_kpointer_[alpha][SpongePosition] = InitialM;
                        StartSolution.Map_alpha_PointertoCreatedIndex[alpha][SpongePosition] = StartSolution.LocalCreatedIndex[SpongePosition];
                    }
                }

            }); // End Parallel Section calculating initial center positions
        });

        SystemCoG.sumoverthreadsandmpi();
		for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
		{
			Initial_Y[VectorIndex] = SystemCoG.Totalmean[VectorIndex];
			StartSolution.Y_k_i_[FirstRealCluster][VectorIndex] = Initial_Y[VectorIndex];
		}
		Box<double[]> tempRef_Object = new Box<>(StartSolution.Sigma_k_i_[FirstRealCluster]);
		Program.CalculateSigma(Initial_Y, tempRef_Object);
		StartSolution.Sigma_k_i_[FirstRealCluster] = tempRef_Object.content;

		GlobalReductions.FindDoubleMean InitialAverages = new GlobalReductions.FindDoubleMean(DAVectorUtility.ThreadCount);

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int LocalPointIndex = beginpoint; LocalPointIndex < indexlen + beginpoint; LocalPointIndex++) {
                    InitialAverages.addapoint(threadIndex,
                            DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(LocalPointIndex,
                                    FirstRealCluster, StartSolution)
                    );
                }

            }); // End Parallel Section calculating average distance
        });
        InitialAverages.sumoverthreadsandmpi();

		//  Estimate of Initial Temperature is average scaled squared distance
		//  Fudge factor of 1.5 over estimated critical temperature for first split
		StartSolution.Temperature = 1.5 * InitialAverages.Totalmean;
		if (Program.Tminimum > 0.0)
		{
			VectorAnnealIterate.Tmin = Program.Tminimum;
		}
		else
		{
			VectorAnnealIterate.Tmin = StartSolution.Temperature / Math.abs(Program.Tminimum);
		}
		Program.ActualStartTemperature = StartSolution.Temperature; // For Output
		Program.TargetEndTemperature = VectorAnnealIterate.Tmin; // For Output
		DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global + " Initial Temperature " + String.format("%1$5.4f", StartSolution.Temperature) + " Minimum " + String.format("%1$5.4f", VectorAnnealIterate.Tmin));

		StartSolution.ActualCoolingFactor = Program.InitialCoolingFactor;

		StartSolution.PairwiseHammy = 0.0;
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			StartSolution.PairwiseHammy += 0.5 * ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.runningSolution.C_k_[RealClusterIndex];
		}

    } // End InitializeSolution

	//  This assumes essentially that Sponge Used   and is cluster 0 in both files
	//  Assumes that second run processes sponge from first
	public static void Restart() throws MPIException {
		//  Reset Input File selection
		Program.SelectedInputLabel = Program.RestartSelectedInputLabel;
		Program.InputFileType = 1;

		// Read old Center Positions
		Program.Replicate = 1;
		ParallelClustering.runningSolution.Ncent_Global = 0;
		if (ParallelClustering.runningSolution.SpongeCluster >= 0)
		{
			ParallelClustering.runningSolution.Ncent_Global = 1;
		}
		--ParallelClustering.runningSolution.Ncent_ThisNode;
		DAVectorReadData.ReadDataFromFile(Program.RestartClusterFile, 2);
		int InitialClusterCount = ParallelClustering.runningSolution.Ncent_Global;
		int ExtraClusterCount = 0;
		if (Program.config.LabelFile.length() > 0)
		{
			DAVectorReadData.ReadDataFromFile(Program.config.LabelFile, 1);
			DAVectorReadData.ReadDataFromFile(Program.config.LabelFile, 2);
			ExtraClusterCount = ParallelClustering.runningSolution.Ncent_Global - InitialClusterCount;
		}
		ParallelClustering.runningSolution.Ncent_ThisNode = ParallelClustering.runningSolution.Ncent_Global;
		Program.InitialNcent = ParallelClustering.runningSolution.Ncent_Global;

		// Set up Clusters including CreatedIndex. This is in Global not Distributed mode
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			ParallelClustering.runningSolution.P_k_[RealClusterIndex] = 1.0 / ((double) ParallelClustering.runningSolution.Ncent_ThisNode);
			ParallelClustering.runningSolution.FreezingMeasure_k_[RealClusterIndex] = 0.0;
			ParallelClustering.runningSolution.SplitPriority_k_[RealClusterIndex] = -1;
			ParallelClustering.runningSolution.Splittable_k_[RealClusterIndex] = 0;
			ParallelClustering.runningSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
			ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] = 0;
			if (RealClusterIndex != ParallelClustering.runningSolution.SpongeCluster)
			{
				ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] = 1;
				Box<double[]> tempRef_Object = new Box<>(ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex]);
				Program.CalculateSigma(ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex], tempRef_Object);
				ParallelClustering.runningSolution.Sigma_k_i_[RealClusterIndex] = tempRef_Object.content;
			}
			int CreatedIndex = ClusteringSolution.SetCreatedIndex(RealClusterIndex); // Sets LocalCreatedIndex
			ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex] = 0;
		}

		ParallelClustering.runningSolution.Temperature = Program.RestartTemperature;
		if (Program.Tminimum > 0.0)
		{
			VectorAnnealIterate.Tmin = Program.Tminimum;
		}
		else
		{
			VectorAnnealIterate.Tmin = ParallelClustering.runningSolution.Temperature / Math.abs(Program.Tminimum);
		}
		Program.ActualStartTemperature = ParallelClustering.runningSolution.Temperature; // For Output
		Program.TargetEndTemperature = VectorAnnealIterate.Tmin; // For Output

		ParallelClustering.runningSolution.SetActiveClusters();

		//  Set initial clusters for every point!
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int NumClustersperPoint = 0;
                    int AssignedCluster = Program.PointLabel[alpha];
                    if (AssignedCluster >= ParallelClustering.runningSolution.Ncent_Global) {
                        DAVectorUtility.printAndThrowRuntimeException(
                                " Illegal Cluster Number " + AssignedCluster + " At Point " + (alpha + DAVectorUtility.PointStart_Process) + " Total " + ParallelClustering.runningSolution.Ncent_Global);

                    }
                    if (ParallelClustering.runningSolution.SpongeCluster >= 0) {
                        // Always have a sponge entry for each point
                        NumClustersperPoint = 1;
                        ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][0] = ParallelClustering.runningSolution.LocalCreatedIndex[ParallelClustering.runningSolution.SpongeCluster];
                        ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][0] = 1.0;
                    }
                    if ((AssignedCluster != 0) || (ParallelClustering.runningSolution.SpongeCluster < 0)) {
                        ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][NumClustersperPoint] = ParallelClustering.runningSolution.LocalCreatedIndex[AssignedCluster];
                        ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][NumClustersperPoint] = 1.0;
                        if (ParallelClustering.runningSolution.SpongeCluster >= 0) {
                            ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][0] = 0.0;
                        }
                        ++NumClustersperPoint;
                    }
                    ParallelClustering.runningSolution.NumClusters_alpha_[alpha] = NumClustersperPoint;
                }

            }); //  End Parallel Section
        });

        ParallelClustering.runningSolution.FindOccupationCounts();
		double wgt = 1.0 / DAVectorUtility.PointCount_Global;
		for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
		{
			if (ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex] == 0)
			{
				// Exception e = DAVectorUtility.printAndThrowRuntimeException(" No points for Cluster Number " + RealClusterIndex.ToString() + " Total " + ParallelClustering.RunningDAVectorSolt.Ncent_Global.ToString());
				//
				DAVectorUtility.SALSAPrint(0, " No points for Cluster Number " + RealClusterIndex + " Total " + ParallelClustering.runningSolution.Ncent_Global);
			}
			ParallelClustering.runningSolution.C_k_[RealClusterIndex] = (double)ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex];
			ParallelClustering.runningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.runningSolution.C_k_[RealClusterIndex];
		}
		String message = "Restart Set up Initial Clusters " + InitialClusterCount + " Extra Clusters " + ExtraClusterCount;
		if (ParallelClustering.runningSolution.SpongeCluster >= 0)
		{
			message += " Sponge Count " + ParallelClustering.runningSolution.OccupationCounts_k_[ParallelClustering.runningSolution.SpongeCluster];
		}
		else
		{
			message += " No Sponge";
		}
		DAVectorUtility.SALSAPrint(1, message);

		if ((ParallelClustering.runningSolution.SpongeCluster != -1) && (Program.SpongePoption == 1))
		{
			wgt = 1.0 - ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster];
			ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster] = Program.SpongePWeight / ParallelClustering.runningSolution.Ncent_Global;
			wgt = (1.0 - ParallelClustering.runningSolution.P_k_[ParallelClustering.runningSolution.SpongeCluster]) / wgt;
			for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; RealClusterIndex++)
			{
				if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
				{
					continue;
				}
				if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] > 0)
				{
					ParallelClustering.runningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.runningSolution.P_k_[RealClusterIndex];
				}
			}
		}
		ParallelClustering.runningSolution.SetClusterWidths();
		double AverageClusterwidth = ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth;
		double[] oldwidth = new double[Program.ParameterVectorDimension];
		if (Program.CalculateIndividualWidths)
		{
            System.arraycopy(ParallelClustering.runningSolution.AverageWidth, 0, oldwidth, 0,
                    Program.ParameterVectorDimension);
		}
		Program.ActualSpongeTemperature = ParallelClustering.runningSolution.Temperature;
		Program.ActualWidth = AverageClusterwidth;
		DAVectorUtility.InterimTiming();
		Program.TimeatSponge = DAVectorUtility.HPDuration;

		VectorAnnealIterate.DistributeNextTime = true;
		ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;

		VectorAnnealIterate.CompleteCleanUp();
		ParallelClustering.runningSolution.SetClusterWidths();
		String widthmessage = "";
		if (Program.CalculateIndividualWidths)
		{
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				widthmessage += VectorIndex + " Component Widths " + String.format("%1$5.4E", oldwidth[VectorIndex]) + " " + String.format("%1$5.4E", ParallelClustering.runningSolution.AverageWidth[VectorIndex]) + " ";
			}
		}
		else
		{
			widthmessage = "Total width " + String.format("%1$5.4E", ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth);
		}
		DAVectorUtility.SALSAPrint(1, widthmessage);

	} // End Restart


	public static boolean CheckCloseClusters(int LocalActiveClusterIndex1)
	{
		int RealClusterIndex1 = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex1];
		for (int LocalActiveClusterIndex2 = 0; LocalActiveClusterIndex2 < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex2++)
		{
			if (LocalActiveClusterIndex2 == LocalActiveClusterIndex1)
			{
				continue;
			}
			int RealClusterIndex2 = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex2];
			if (RealClusterIndex2 == ParallelClustering.runningSolution.SpongeCluster)
			{
				continue;
			}
			int test = ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex1] - ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex2];
			if (test > 0)
			{
				continue;
			}
			if ((test == 0) && (LocalActiveClusterIndex1 < LocalActiveClusterIndex2))
			{
				continue;
			}
			double TestDistce = DAVectorParallelism.getSquaredScaledDistanceTweenActiveClusters(LocalActiveClusterIndex1, LocalActiveClusterIndex2, ParallelClustering.runningSolution);
			if (TestDistce < Program.ScaledSquaredDistanceatClosenessTest)
			{
				return true;
			}
		}
		return false;

	} // End CheckCloseClusters

	// Return False if current final solution has a problem
	//  Return True if there is no problem or no way of fixing problem
	public static boolean CheckValidSolution(boolean UseBestSolution) throws MPIException {
		if (UseBestSolution)
		{
			if (!ParallelClustering.bestSolution.SolutionSet)
			{
				DAVectorUtility.SALSAPrint(1, "Solution Not changed as no best solution");
			}
			if ((ParallelClustering.bestSolution.Ncent_ThisNode != ParallelClustering.runningSolution.Ncent_ThisNode) || (ParallelClustering.bestSolution.IterationSetAt != ParallelClustering.runningSolution.IterationSetAt))
			{
				DAVectorUtility.SALSAPrint(1, " Best Solution Not Used to restart as Number of Centers " + ParallelClustering.bestSolution.Ncent_ThisNode + " Different from Running Solution with " + ParallelClustering.runningSolution.Ncent_ThisNode + " Best Iteration# " + ParallelClustering.bestSolution.IterationSetAt + " Current " + ParallelClustering.runningSolution.IterationSetAt);
			}
			else
			{
				boolean hammychange = ParallelClustering.runningSolution.PairwiseHammy > ParallelClustering.bestSolution.PairwiseHammy;
				hammychange = DAVectorUtility.SynchronizeMPIBoolean(hammychange);
				if (hammychange)
				{
					ClusteringSolution.CopySolution(ParallelClustering.bestSolution, ParallelClustering.runningSolution);
					DistributedClusteringSolution.ManageMajorSynchronization(true);
					DAVectorUtility.SALSAPrint(1, "Best Solution Used");
					ParallelClustering.runningSolution.OldHammy = 0.0;
					ParallelClustering.runningSolution.PairwiseHammy = 0.0;
					ParallelClustering.runningSolution.Eigenvectorset = false;
					ParallelClustering.runningSolution.ClustertoSplit = -1;
				}
			}
		}

		boolean CheckCloseness = false;
		if (ParallelClustering.runningSolution.Temperature < Program.TemperatureforClosenessTest)
		{
			if ((VectorAnnealIterate.TemperatureatLastCloseTest < 0.0) || (ParallelClustering.runningSolution.Temperature < VectorAnnealIterate.TemperatureatLastCloseTest - 0.25))
			{
				CheckCloseness = true;
			}
		}
		CheckCloseness = DAVectorUtility.SynchronizeMPIBoolean(CheckCloseness);
		if (CheckCloseness)
		{
			VectorAnnealIterate.TemperatureatLastCloseTest = ParallelClustering.runningSolution.Temperature;
		}

		boolean[] RemoveCluster = new boolean[ParallelClustering.runningSolution.Ncent_ThisNode];
		int Numbertoosmall01 = 0;
		int Numbertoosmall2 = 0;
		double ProbabilitySum01 = 0.0;
		double ProbabilitySum2 = 0.0;
		String documentit = "";
		ParallelClustering.runningSolution.FindOccupationCounts();
		double C_k_Test = Program.MinimumCountforCluster_C_k;
		if (ParallelClustering.runningSolution.SpongeCluster >= 0)
		{
			C_k_Test = Program.MinimumCountforCluster_C_kwithSponge;
		}
		int NumberSmall_C = 0;
		int NumberSmall_OccCount = 0;
		int NumberClose = 0;
		for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
			RemoveCluster[RealClusterIndex] = false;
			if (RealClusterIndex == ParallelClustering.runningSolution.SpongeCluster)
			{
				if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
				{
					ProbabilitySum2 += ParallelClustering.runningSolution.P_k_[RealClusterIndex];
				}
				else
				{
					ProbabilitySum01 += ParallelClustering.runningSolution.P_k_[RealClusterIndex];
				}
				continue;
			}
			if (!(ParallelClustering.runningSolution.DistributedExecutionMode && (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] < 2)))
			{
				if (ParallelClustering.runningSolution.C_k_[RealClusterIndex] < C_k_Test)
				{
					RemoveCluster[RealClusterIndex] = true;
					++NumberSmall_C;
				}
				else if (ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex] < Program.MinimumCountforCluster_Points)
				{
					RemoveCluster[RealClusterIndex] = true;
					++NumberSmall_OccCount;
				}
				if (CheckCloseness && (!RemoveCluster[RealClusterIndex]))
				{
					RemoveCluster[RealClusterIndex] = CheckCloseClusters(LocalActiveClusterIndex);
					if (RemoveCluster[RealClusterIndex])
					{
						++NumberClose;
					}
				}
			}
		}

		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - int - sum
//                NumberSmall_C = DAVectorUtility.MPI_communicator.<Integer>Allreduce(NumberSmall_C, Operation<Integer>.Add);
                NumberSmall_C = DAVectorUtility.mpiOps.allReduce(NumberSmall_C, MPI.SUM);
                // Note - MPI Call - Allreduce - int - sum
//                NumberSmall_OccCount = DAVectorUtility.MPI_communicator.<Integer>Allreduce(NumberSmall_OccCount, Operation<Integer>.Add);
                NumberSmall_OccCount = DAVectorUtility.mpiOps.allReduce(NumberSmall_OccCount, MPI.SUM);
                // Note - MPI Call - Allreduce - int - sum
//                NumberClose = DAVectorUtility.MPI_communicator.<Integer>Allreduce(NumberClose, Operation<Integer>.Add);
                NumberClose = DAVectorUtility.mpiOps.allReduce(NumberClose, MPI.SUM);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		}
		else
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - boolean []
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Boolean>Broadcast(tempRef_RemoveCluster, 0);
                DAVectorUtility.mpiOps.broadcast(RemoveCluster, 0);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
		}
		Program.TotalClustersDeleted_CSmall += NumberSmall_C;
		Program.TotalClustersDeleted_OccCount += NumberSmall_OccCount;
		Program.TotalClustersDeleted_Close += NumberClose;

		for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
			if (RemoveCluster[RealClusterIndex])
			{
				if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
				{
					++Numbertoosmall2;
				}
				else
				{
					++Numbertoosmall01;
				}
				documentit += RealClusterIndex + "(" + ParallelClustering.runningSolution.LocalCreatedIndex[RealClusterIndex] + ") C_k_ " + String.format("%1$2.1f", ParallelClustering.runningSolution.C_k_[RealClusterIndex]) + " Points " + ParallelClustering.runningSolution.OccupationCounts_k_[RealClusterIndex];
				if (Program.ParameterVectorDimension == 2)
				{
					documentit += " Y " + String.format("%1$4.3f", ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][0]) + " " + String.format("%1$4.3f", ParallelClustering.runningSolution.Y_k_i_[RealClusterIndex][1]) + " * ";
				}
				else
				{
					documentit += " * ";
				}
				continue;
			}
			if (ParallelClustering.runningSolution.LocalStatus[RealClusterIndex] == 2)
			{
				ProbabilitySum2 += ParallelClustering.runningSolution.P_k_[RealClusterIndex];
			}
			else
			{
				ProbabilitySum01 += ParallelClustering.runningSolution.P_k_[RealClusterIndex];
			}
        }

		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            if (DAVectorUtility.MPI_Size > 1){
                // Note - MPI Call - Allreduce - int - sum
//                Numbertoosmall2 = DAVectorUtility.MPI_communicator.<Integer>Allreduce(Numbertoosmall2, Operation<Integer>.Add);
                Numbertoosmall2 = DAVectorUtility.mpiOps.allReduce(Numbertoosmall2, MPI.SUM);
                // Note - MPI Call - Allreduce - int - sum
//                ProbabilitySum2 = DAVectorUtility.MPI_communicator.<Double>Allreduce(ProbabilitySum2, Operation<Double>.Add);
                ProbabilitySum2 = DAVectorUtility.mpiOps.allReduce(ProbabilitySum2, MPI.SUM);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
		}
		double ProbabilitySum = ProbabilitySum01 + ProbabilitySum2;
		int GlobalNumberTooSmall = Numbertoosmall01 + Numbertoosmall2;
		if (GlobalNumberTooSmall == 0)
		{
			return true;
		}

		//  Need to change solution as one or more clusters too small
		if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
		{
			if (ParallelClustering.runningSolution.DistributedExecutionMode)
			{
				DAVectorUtility.SALSASyncPrint(1, GlobalNumberTooSmall + " Clusters Too Small ", documentit);
			}
			else
			{
				DAVectorUtility.SALSAPrint(1, GlobalNumberTooSmall + " Clusters Too Small " + documentit);
			}
		}

		for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
		{
			int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
			if (!RemoveCluster[RealClusterIndex])
			{
				if (Program.ContinuousClustering)
				{
					ParallelClustering.runningSolution.P_k_[RealClusterIndex] = ParallelClustering.runningSolution.P_k_[RealClusterIndex] / ProbabilitySum;
				}
            }
		}
		if (!ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			int oldNcent = ParallelClustering.runningSolution.Ncent_ThisNode;
			for (int clusterindex = oldNcent - 1; clusterindex >= 0; clusterindex--)
			{
				if (RemoveCluster[clusterindex])
				{
					ParallelClustering.runningSolution.RemoveCluster(clusterindex); // This both changes cluster count and shifts up clusters
				}
			}
		}
		else
		{
			for (int clusterindex = 0; clusterindex < ParallelClustering.runningSolution.Ncent_ThisNode; clusterindex++)
			{
				if (RemoveCluster[clusterindex])
				{
					ParallelClustering.runningSolution.LocalStatus[clusterindex] = -1;
				}
			}
		}
		DistributedClusteringSolution.ManageMajorSynchronization(true);
		ClusteringSolution.CopySolution(ParallelClustering.runningSolution, ParallelClustering.bestSolution);
		return false;

	} // End CheckValidSolution

	public static boolean AddSpongeCluster() throws MPIException {
		if (ParallelClustering.runningSolution.SpongeCluster >= 0)
		{
			return false;
		}
		int CurrentMax = ParallelClustering.runningSolution.Ncent_ThisNode;
		if (CurrentMax >= Program.maxNcentperNode)
		{
			return false;
		}

		double AverageClusterwidth = ParallelClustering.runningSolution.TotaloverVectorIndicesAverageWidth;

		boolean Addasponge = (Program.CreateSpongeScaledSquaredWidth > 0.0) && (AverageClusterwidth <= Program.CreateSpongeScaledSquaredWidth);
		if (!Addasponge)
		{
			Addasponge = (Program.SpongeTemperature1 > 0.0) && (ParallelClustering.runningSolution.Temperature < Program.SpongeTemperature1);
		}
		Addasponge = DAVectorUtility.SynchronizeMPIBoolean(Addasponge);
		if (!Addasponge)
		{
			return false;
		}
		VectorAnnealIterate.OutputClusteringResults("StartSponge");

		//  Note Cluster numbers different in each node if distributed execution. Created Index identical across nodes
		ParallelClustering.runningSolution.SpongeCluster = CurrentMax;
		int SpongeCreatedIndex = 0;
		if (ParallelClustering.runningSolution.DistributedExecutionMode)
		{
			if (DAVectorUtility.MPI_Rank == 0)
			{
				SpongeCreatedIndex = ClusteringSolution.SetCreatedIndex(CurrentMax);
			}
			DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int
            if (DAVectorUtility.MPI_Size > 1){
//                DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_SpongeCreatedIndex, 0);
                SpongeCreatedIndex = DAVectorUtility.mpiOps.broadcast(SpongeCreatedIndex,0);
            }
			DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
			if (DAVectorUtility.MPI_Rank != 0)
			{
				ClusteringSolution.UniversalMapping[SpongeCreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, 1 + CurrentMax);
				ParallelClustering.runningSolution.LocalCreatedIndex[CurrentMax] = SpongeCreatedIndex;
			}
		}
		else
		{
			SpongeCreatedIndex = ClusteringSolution.SetCreatedIndex(CurrentMax);
		}

		Program.ActualSpongeTemperature = ParallelClustering.runningSolution.Temperature;
		Program.ActualWidth = AverageClusterwidth;
		DAVectorUtility.InterimTiming();
		Program.TimeatSponge = DAVectorUtility.HPDuration;
		DAVectorUtility.SALSAPrint(1, "Sponge added with Created Index " + SpongeCreatedIndex + " Time " + String.format("%1$5.4E", Program.TimeatSponge));

		//  Add Sponge to every point!
        // Note - parallel for
        final int SpongeCreatedIndexLoopVar = SpongeCreatedIndex;
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    int NewIndirectClusterIndex = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][NewIndirectClusterIndex] = SpongeCreatedIndexLoopVar;
                    ++ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][NewIndirectClusterIndex] = 0.0;
                }

            }); //  End Parallel Section
        });

        //  Set up Sponge Cluster
		ParallelClustering.runningSolution.Splittable_k_[CurrentMax] = 0;
		ParallelClustering.runningSolution.SplitPriority_k_[CurrentMax] = 0;
		ParallelClustering.runningSolution.C_k_[CurrentMax] = 0.0;
		ParallelClustering.runningSolution.ClusterScaledSquaredWidth_k_[CurrentMax] = 0.0;
		ParallelClustering.runningSolution.FreezingMeasure_k_[CurrentMax] = 0.0;
		ParallelClustering.runningSolution.LocalStatus[CurrentMax] = 0;
		ParallelClustering.runningSolution.LocalSplitCreatedIndex[CurrentMax] = 0;

		ParallelClustering.runningSolution.P_k_[CurrentMax] = Program.SpongePWeight / ParallelClustering.runningSolution.Ncent_Global;
		double wgt = (1.0 - ParallelClustering.runningSolution.P_k_[CurrentMax]);
		for (int ClusterIndex = 0; ClusterIndex < ParallelClustering.runningSolution.Ncent_ThisNode; ClusterIndex++)
		{
			if (ParallelClustering.runningSolution.LocalStatus[ClusterIndex] < 0)
			{
				continue;
			}
			ParallelClustering.runningSolution.P_k_[ClusterIndex] = wgt * ParallelClustering.runningSolution.P_k_[ClusterIndex];
		}

		ParallelClustering.runningSolution.DiffMalpha_k_Set = -1;
		ParallelClustering.runningSolution.Ncent_ThisNode = CurrentMax + 1;
		++ParallelClustering.runningSolution.Ncent_Global;
		DistributedClusteringSolution.ManageMajorSynchronization(true);
		return true;

	} // End AddSpongeCluster()

	public static void CalculateClusterStatus() throws MPIException {
		if (Program.DoLCMS)
		{
			LCMSAnalyze.LCMSCalculateClusterStatus();
		}
	}

	// Unlike Earlier versions, cluster labels are 0-based
	public static void OutputClusteringResults(String FileLabel) throws MPIException {
		boolean DoFileOutput = true;
		if (FileLabel.length() == 0)
		{
			DoFileOutput = false;
		}

		Program.ClusterNumberOutput = ParallelClustering.runningSolution.Ncent_Global;

		//  Set Global Counts and Positions
		double[][] GlobalCenters = new double[ParallelClustering.runningSolution.Ncent_Global][];
		for (int GlobalClusterIndex = 0; GlobalClusterIndex < ParallelClustering.runningSolution.Ncent_Global; GlobalClusterIndex++)
		{
			GlobalCenters[GlobalClusterIndex] = new double[Program.ParameterVectorDimension];
		}
		ClusteringSolution.SetGlobalClusterNumbers();

		// Generate Cluster Labels (cluster number) from current solution
		int[] labels = new int[DAVectorUtility.PointCount_Process];
		String[] Extralabels = new String[DAVectorUtility.PointCount_Process];

		//  Parallel Section setting cluster labels
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
            {
                int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++) {
                    double distmax = -1.0;
                    int NearestClusterCreatedIndex = 0;
                    int IndirectSize = ParallelClustering.runningSolution.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++) {
                        if (ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] > distmax) {
                            distmax = ParallelClustering.runningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            NearestClusterCreatedIndex = ParallelClustering.runningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        }
                    }
                    int position = ClusteringSolution.UniversalMapping[NearestClusterCreatedIndex].GlobalClusterNumber;
                    if (position < 0) {
                        DAVectorUtility.printAndThrowRuntimeException(
                                " CreatedIndex has no Global Position " + NearestClusterCreatedIndex);

                    }
                    labels[alpha] = position;
                    Extralabels[alpha] = "";
                    if (Program.CompareSolution > 0) {
                        // Add older clustering labels
                        int MedeaIndex = Program.MedeaClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int MclustIndex = Program.MclustClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int GoldenIndex = Program.GoldenPeaks.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int OurIndex = Program.OurClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        Extralabels[alpha] = OurIndex + " " + MedeaIndex + " " + MclustIndex + " " + GoldenIndex;
                    }
                }

            }); // End Parallel Section setting cluster label
        });


        String directory = Paths.get(Program.config.ClusterFile).getParent().toString();
		String file = Files.getNameWithoutExtension(Program.config.ClusterFile) + FileLabel + "-M" + Program.maxNcentperNode + "-C" + ParallelClustering.runningSolution.Ncent_Global + "." + Files.getFileExtension(
                Program.config.ClusterFile);
		String ClusternumberFileName = Paths.get(directory, file).toString();

		int MPItag = 100;
		int MPItag2 = 101;
		if (DAVectorUtility.MPI_Rank == 0)
		{
			if (DoFileOutput)
			{
				VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, labels, Program.PointPosition, Program.PointOriginalExprment, Extralabels, DAVectorUtility.PointCount_Process, DAVectorUtility.PointStart_Process, false);
			}

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, DAVectorUtility.ThreadCount - 1, (threadIndex) ->
                {
                    int indexlen = DAVectorUtility.PointsperThread[threadIndex];
                    int beginpoint = DAVectorUtility.StartPointperThread[threadIndex] - DAVectorUtility.PointStart_Process;
                    System.arraycopy(labels, beginpoint, Program.ClusterAssignments,
                            beginpoint + DAVectorUtility.PointStart_Process, indexlen + beginpoint - beginpoint);

                }); // End Parallel Section setting cluster assignments in process 0
            });

            // Note - MPI Call Block (with MPI_Size = 1 this loop will not happen)
			for (int MPISource = 1; MPISource < DAVectorUtility.MPI_Size; MPISource++)
			{
                // Note - MPI Call - Receive - MPIPacket<Integer>
				MPIPacket fromsource = DAVectorUtility.mpiOps.receive(MPISource, MPItag, MPIPacket.Type.Integer);
				int AwayArraySize = fromsource.getNumberOfPoints();
                int firstPoint = fromsource.getFirstPoint();
                for (int index = 0; index < AwayArraySize; index++)
				{
                    Program.ClusterAssignments[index + firstPoint] = fromsource.getMArrayIntAt(index);
				}

				double[][] AwayPointPositions = new double[AwayArraySize][];
				int[] AwayOriginalExprment = new int[AwayArraySize];
				for (int index = 0; index < AwayArraySize; index++)
				{
					AwayPointPositions[index] = new double[Program.ParameterVectorDimension];
				}
				MPIPacket fromsourcedouble;
				MPIPacket fromsouceint;
				fromsouceint = DAVectorUtility.mpiOps.receive(MPISource, MPItag2, MPIPacket.Type.Integer);

				for (int LocalPointIndex = 0; LocalPointIndex < AwayArraySize; LocalPointIndex++)
				{
					AwayOriginalExprment[LocalPointIndex] = fromsouceint.getMArrayIntAt(LocalPointIndex);
				}

				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
					fromsourcedouble = DAVectorUtility.mpiOps.receive(MPISource, MPItag, MPIPacket.Type.Double);
					for (int LocalPointIndex = 0; LocalPointIndex < AwayArraySize; LocalPointIndex++)
					{
						AwayPointPositions[LocalPointIndex][VectorIndex] = fromsourcedouble.getMArrayDoubleAt(LocalPointIndex);
					}
				}
				if (DoFileOutput)
				{
					VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, fromsource::getMArrayIntAt, AwayPointPositions,AwayOriginalExprment, Extralabels, AwayArraySize, fromsource.getFirstPoint(), true);
				}
			}

			int CenterLabelSize = ParallelClustering.runningSolution.Ncent_Global;
			if (ParallelClustering.runningSolution.SpongeCluster >= 0)
			{
				--CenterLabelSize;
			}
			int[] CenterLabel = new int[CenterLabelSize];
			double[][] CenterPositions = new double[CenterLabelSize][];
			String[] ExtraCenterlabels = new String[CenterLabelSize];

			int decrement = 0;
			for (int StrippedClusterIndex = 0; StrippedClusterIndex < CenterLabelSize; StrippedClusterIndex++)
			{
				CenterPositions[StrippedClusterIndex] = new double[Program.ParameterVectorDimension];
				if (ClusteringSolution.TotalClusterSummary.SpongeCluster == StrippedClusterIndex)
				{
					decrement = 1;
				}
				int ClusterIndex = StrippedClusterIndex + decrement;
                System.arraycopy(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], 0,
                        CenterPositions[StrippedClusterIndex], 0, Program.ParameterVectorDimension);
				CenterLabel[StrippedClusterIndex] = Math.max(100000000, CenterLabelSize);
				int OtherIndices = -1;
				ExtraCenterlabels[StrippedClusterIndex] = "";
				if (Program.CompareSolution > 0)
				{
                    ExtraCenterlabels[StrippedClusterIndex] = ClusterIndex + " " + OtherIndices + " " + OtherIndices + " " + OtherIndices;
				}

			}

			// In first column for centers, we output cluster number + Total Point Count
			if (DoFileOutput)
			{
				VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, CenterLabel, CenterPositions, new int[0], ExtraCenterlabels, CenterLabelSize, DAVectorUtility.PointCount_Global, true);
			}

		} // End Root Process that receives cluster assignments from afar
		else
		{
            // Note - MPI Call Block (with MPI_Size = 1 this else will not happen)
			MPIPacket tosend = MPIPacket.newIntegerPacket(DAVectorUtility.PointCount_Process);
			tosend.setFirstPoint(DAVectorUtility.PointStart_Process);
			tosend.setNumberOfPoints(DAVectorUtility.PointCount_Process);
            for (int i = 0; i < labels.length; ++i){
                tosend.setMArrayIntAt(i,labels[i]);
            }
            // Note - MPI Call - Send - MPIPacket<Integer>
			DAVectorUtility.mpiOps.send(tosend, 0, MPItag);
			MPIPacket tosenddouble = MPIPacket.newDoublePacket(DAVectorUtility.PointCount_Process);
			MPIPacket tosendint = MPIPacket.newIntegerPacket(DAVectorUtility.PointCount_Process);

			tosendint.setFirstPoint(DAVectorUtility.PointStart_Process);
			tosendint.setNumberOfPoints(DAVectorUtility.PointCount_Process);

			for (int LocalPointIndex = 0; LocalPointIndex < DAVectorUtility.PointCount_Process; LocalPointIndex++)
			{
				tosendint.setMArrayIntAt(LocalPointIndex, Program.PointOriginalExprment[LocalPointIndex]);
			}
			DAVectorUtility.mpiOps.send(tosendint, 0, MPItag2);

			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				tosenddouble.setFirstPoint(DAVectorUtility.PointStart_Process);
				tosenddouble.setNumberOfPoints(DAVectorUtility.PointCount_Process);
				for (int LocalPointIndex = 0; LocalPointIndex < DAVectorUtility.PointCount_Process; LocalPointIndex++)
				{
					tosenddouble.setMArrayDoubleAt(LocalPointIndex, Program.PointPosition[LocalPointIndex][VectorIndex]);
				}
                // Note - MPI Call - Send - MPIPacket<Double>
				DAVectorUtility.mpiOps.send(tosenddouble, 0, MPItag);
			}
		}
        // Note - MPI Call - Broadcast - int
//		DAVectorUtility.MPI_communicator.<Integer>Broadcast(tempRef_ClusterAssignments, 0);
        DAVectorUtility.mpiOps.broadcast(Program.ClusterAssignments, 0);
        // Note - MPI Call - Barrier
		DAVectorUtility.mpiOps.barrier();

        /* Restart Timer - requires reset and start */
		DAVectorUtility.PreciseTimer.reset();
		DAVectorUtility.PreciseTimer.start();

	} // End OutputClusteringResults

	public static void Output3DClusterLabels(String extraname)
	{
		if (DAVectorUtility.MPI_Rank != 0)
		{
			return;
		}
		if (Program.RW3DData <= 0)
		{
			return;
		}

        String directory = Paths.get(Program.config.ClusterFile).getParent().toString();
		String file = Files.getNameWithoutExtension(Program.config.ClusterFile) + "-Plot3D" + "-M" + Program.maxNcentperNode + "-C" + ParallelClustering.runningSolution.Ncent_Global + "." + Files.getFileExtension(
                Program.config.ClusterFile);
		String ClusternumberFileName = Paths.get(directory, extraname + "-" + file).toString();

		Write3DClusterFile(ClusternumberFileName, Program.ClusterAssignments, 0, false);
	} // End Output3DClusterLabels()

    public static void WriteClusterFile(String fname, int[] labels, double[][] PointPositions,int[] experimentNumbers, String[] ExtraLabels, int dataPoints, int startposition, boolean append){
        WriteClusterFile(fname, i -> labels[i], PointPositions,experimentNumbers, ExtraLabels, dataPoints, startposition, append);
    }
	public static void WriteClusterFile(String fname, IntArray labels, double[][] PointPositions,int[] experimentNumbers, String[] ExtraLabels, int dataPoints, int startposition, boolean append)
	{

        OpenOption mode = append ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;

        try (PrintWriter writer = new PrintWriter(
                java.nio.file.Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset(), mode), true)) {
			System.out.println("experimentNumbers ----------------------- " + experimentNumbers.length);
			System.out.println("data points ----------------------- " + dataPoints);
			System.out.println("PointPositions points ----------------------- " + PointPositions.length);
			System.out.println("PointOriginalIndex points ----------------------- " + Program.PointOriginalIndex.length);
			for (int i = 0; i < dataPoints; i++) {
                String line = String.valueOf(i + startposition);
				line += " " + experimentNumbers[i];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++) {
                    line += " " + String.format("%1$7.6f", PointPositions[i][VectorIndex]);
                }
                if (Program.ParameterVectorDimension == 2) {
                    line += " 0.0";
                }
                line += " " + labels.get(i) + " " + ExtraLabels[i];
                writer.println(line);
            }
            writer.close();
        } catch (IOException e) {
            System.err.format("Failed writing cluster file due to I/O exception: %s%n", e);
        }

		if (Program.SigmaMethod > 0)
		{
			fname = fname.replace(".txt", "Scaled.txt");
            try (PrintWriter writer = new PrintWriter(
                    java.nio.file.Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset(), mode), true)) {
                double tmp;
                for (int i = 0; i < dataPoints; i++)
                {
                    String line = String.valueOf(i + startposition);
					line += " " + Program.PointOriginalExprment[i];
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        tmp = PointPositions[i][VectorIndex];
                        if ((Program.SigmaMethod >= 2) && (VectorIndex == 0) && (tmp != 0.0))
                        {
                            tmp = Math.log(tmp);
                        }
                        tmp = tmp / Program.SigmaVectorParameters_i_[VectorIndex];
                        line += " " + String.format("%1$7.6f", tmp);
                    }
                    if (Program.ParameterVectorDimension == 2)
                    {
                        line += " 0.0";
                    }
                    line += " " + labels.get(i);
                    writer.println(line);
                }
                writer.close();
            } catch (IOException e) {
                System.err.format("Failed writing scaled cluster results due to I/O exception: %s%n", e);
            }
		}
	} // End WriteClusterFile

	public static void Write3DClusterFile(String fname, int[] labels, int startposition, boolean append)
	{
		OpenOption mode = append ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;
        try (PrintWriter writer = new PrintWriter(
                java.nio.file.Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset(), mode), true)) {
            for (int i = 0; i < DAVectorUtility.PointCount_Global; i++)
            {
                String line = String.valueOf(i + startposition);
                for (int VectorIndex = 0; VectorIndex < Program.RW3DData; VectorIndex++)
                {
                    line += " " + String.format("%1$7.6f", Program.FullPoint3DPosition[i][VectorIndex]);
                }
                if (Program.RW3DData == 2)
                {
                    line += " 0.0";
                }
                line += " " + labels[i];
                writer.println(line);
            }
        } catch (IOException e) {
            System.err.format("Failed writing 3D cluster file due to I/O exception: %s%n", e);
        }
	} // End Write3DClusterFile

	// Unlike Earlier versions, cluster labels are 0-based
	public static void SimpleOutputClusteringResults(String FileLabel, double[][] ClusterCenters) throws MPIException { // Output using precalculated Cluster Assignment Array
        String directory = Paths.get(Program.config.ClusterFile).getParent().toString();
		String file = Files.getNameWithoutExtension(Program.config.ClusterFile) + FileLabel + "-M" + Program.maxNcentperNode + "-C" + ParallelClustering.runningSolution.Ncent_Global+ "." + Files.getFileExtension(
                Program.config.ClusterFile);
		String ClusternumberFileName = Paths.get(directory, file).toString();

		int MPItag = 100;
		if (DAVectorUtility.MPI_Rank == 0)
		{
			Box<int[]> tempRef_ClusterAssignments = new Box<>(Program.ClusterAssignments);
			Box<double[][]> tempRef_PointPosition = new Box<>(Program.PointPosition);
			VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, tempRef_ClusterAssignments, tempRef_PointPosition, DAVectorUtility.PointCount_Process, DAVectorUtility.PointStart_Process, DAVectorUtility.PointStart_Process, false);
			Program.ClusterAssignments = tempRef_ClusterAssignments.content;
			Program.PointPosition = tempRef_PointPosition.content;

            // Note - MPI Call Block (with MPI_Size = 1 this loop will not happen)
			for (int MPISource = 1; MPISource < DAVectorUtility.MPI_Size; MPISource++)
			{
                // Note - MPI Call - Receive - MPIPacket<Integer>
				MPIPacket fromsource = DAVectorUtility.mpiOps.receive(MPISource, MPItag, MPIPacket.Type.Integer);
				int AwayArraySize = fromsource.getNumberOfPoints();

				double[][] AwayPointPositions = new double[AwayArraySize][];
				for (int index = 0; index < AwayArraySize; index++)
				{
					AwayPointPositions[index] = new double[Program.ParameterVectorDimension];
				}
				MPIPacket fromsourcedouble;
				for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
				{
                    // Note - MPI Call - Receive - MPIPacket<Double>
					fromsourcedouble = DAVectorUtility.mpiOps.receive(MPISource, MPItag, MPIPacket.Type.Double);
					for (int LocalPointIndex = 0; LocalPointIndex < AwayArraySize; LocalPointIndex++)
					{
						AwayPointPositions[LocalPointIndex][VectorIndex] = fromsourcedouble.getMArrayDoubleAt(LocalPointIndex);
					}
				}
				Box<int[]> tempRef_ClusterAssignments2 = new Box<>(Program.ClusterAssignments);
				Box<double[][]> tempRef_AwayPointPositions = new Box<>(AwayPointPositions);
				VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, tempRef_ClusterAssignments2, tempRef_AwayPointPositions, AwayArraySize, fromsource.getFirstPoint(), fromsource.getFirstPoint(), true);
				Program.ClusterAssignments = tempRef_ClusterAssignments2.content;
				AwayPointPositions = tempRef_AwayPointPositions.content;
			}

			int CenterLabelSize = ParallelClustering.runningSolution.Ncent_Global;
			if (ParallelClustering.runningSolution.SpongeCluster >= 0)
			{
				--CenterLabelSize;
			}
			int[] CenterLabel = new int[CenterLabelSize];

			int decrement = 0;
			for (int StrippedClusterIndex = 0; StrippedClusterIndex < CenterLabelSize; StrippedClusterIndex++)
			{
				if (ClusteringSolution.TotalClusterSummary.SpongeCluster == StrippedClusterIndex)
				{
					decrement = 1;
				}
				int ClusterIndex = StrippedClusterIndex + decrement;
				CenterLabel[StrippedClusterIndex] = Math.max(100000000, CenterLabelSize);
			}

			// In first column for centers, we output cluster number + DAVectorUtility.PointCount_Global
			Box<int[]> tempRef_CenterLabel = new Box<>(CenterLabel);
			Box<double[][]> tempRef_ClusterCenters = new Box<>(ClusterCenters);
			VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, tempRef_CenterLabel, tempRef_ClusterCenters, CenterLabelSize, DAVectorUtility.PointCount_Global, 0, true);
			CenterLabel = tempRef_CenterLabel.content;
			ClusterCenters = tempRef_ClusterCenters.content;

		} // End Root Process that receives cluster assignments from afar
		else
		{
            // Note - MPI Call Block (with MPI_Size = 1 this else will not happen)
			MPIPacket tosend = MPIPacket.newIntegerPacket(DAVectorUtility.PointCount_Process);
			tosend.setFirstPoint(DAVectorUtility.PointStart_Process);
			tosend.setNumberOfPoints(DAVectorUtility.PointCount_Process);
            int numberOfPoints = tosend.getNumberOfPoints();
            int firstPoint = tosend.getFirstPoint();
            for (int index = 0; index < numberOfPoints; index++)
			{
                tosend.setMArrayIntAt(index, Program.ClusterAssignments[index + firstPoint]);
			}
            // Note - MPI Call - Send - MPIPacket<Integer>
			DAVectorUtility.mpiOps.send(tosend, 0, MPItag);
			MPIPacket tosenddouble = MPIPacket.newDoublePacket(DAVectorUtility.PointCount_Process);
			for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
			{
				tosenddouble.setFirstPoint(DAVectorUtility.PointStart_Process);
				tosenddouble.setNumberOfPoints(DAVectorUtility.PointCount_Process);
				for (int LocalPointIndex = 0; LocalPointIndex < DAVectorUtility.PointCount_Process; LocalPointIndex++)
				{
					tosenddouble.setMArrayDoubleAt(LocalPointIndex, Program.PointPosition[LocalPointIndex][VectorIndex]);
				}
                // Note - MPI Call - Send - MPIPacket<Double>
				DAVectorUtility.mpiOps.send(tosenddouble, 0, MPItag);
			}
		}

        /* Restart Timer - requires reset and start */
		DAVectorUtility.PreciseTimer.reset();
		DAVectorUtility.PreciseTimer.start();

	} // End SimpleOutputClusterLabels(filename)


	public static void SimpleWriteClusterFile(String fname, Box<int[]> labels, Box<double[][]> PointPositions, int dataPoints, int startposition1, int startposition2, boolean append)
	{
        OpenOption mode = append ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;

        try (PrintWriter writer = new PrintWriter(
                java.nio.file.Files.newBufferedWriter(Paths.get(fname), Charset.defaultCharset(), mode), true)) {
            for (int i = 0; i < dataPoints; i++)
            {
                String line = String.valueOf(i + startposition1);
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    line += " " + String.format("%1$7.6f", PointPositions.content[i][VectorIndex]);
                }
                if (Program.ParameterVectorDimension == 2)
                {
                    line += " 0.0";
                }
                line += " " + labels.content[i + startposition2];
                writer.println(line);
            }
        } catch (IOException e) {
            System.err.format("Failed writing simple cluster file due to I/O exception: %s%n", e);
        }
	} // End SimpleWriteClusterFile

} // End class dist
 // End namespace cluster
