package edu.indiana.soic.spidal.davs;

import com.google.common.base.Strings;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.regex.Pattern;

public class DAVectorReadData {

    // Experiment Dependent code only applicable to Mani's 2016 dataset
    public static void ReadExperimentNumbers(String comparisonClusterFile) {
        int MinSplitSize = 8;
        int ExpermentNumberPosition = 8;
        int SplitPosition = 3;
        DAVectorUtility.SALSAPrint(0, "Reading Experiment Numbers");

        boolean success = false;
        int count = 0;
        if (Strings.isNullOrEmpty(comparisonClusterFile)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try(BufferedReader br = Files.newBufferedReader(Paths.get(comparisonClusterFile),Charset.defaultCharset())){
            String line;
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.printAndThrowRuntimeException("Count " + count + "Illegal data length on Point " +
                            "file " + splits.length + " " + MinSplitSize + " " + line);
                }

                Integer parsedInt = Ints.tryParse(splits[SplitPosition]);
                if (parsedInt == null || parsedInt != Program.SelectedInputLabel) continue;

                Program.ExperimentNumberAssigments[count] = Integer.valueOf(splits[ExpermentNumberPosition]);

                count++;

            }
            DAVectorUtility.SALSAPrint(0, "Experiment Numbers Read " + count);

            success = true;
            br.close();
        } catch (IOException e) {
            System.err.format("Failed reading Points data due to I/O exception: %s%n", e);
        }

        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File Analyze error " + comparisonClusterFile);
        }
    }
    public static void ReadLabelsFromFile(String fname) {
        int MinSplitSize = 8;
        int SplitPosition = 3;
        int MedeaPosition = 7;
        int GoldenIDPosition = 5;
        int MclustPosition = 6;
        int GoldenLabelPosition = 4;
        int ExpermentNumberPosition = 8;

        boolean success = false;
        int count = 0;

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try {
            BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset());
            String line;
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.printAndThrowRuntimeException("Count " + count + "Illegal data length on Point " +
                            "file " + splits.length + " " + MinSplitSize + " " + line);
                }

                Integer parsedInt = Ints.tryParse(splits[SplitPosition]);
                if (parsedInt == null || parsedInt != Program.SelectedInputLabel) continue;

                parsedInt = Ints.tryParse(splits[GoldenIDPosition]);
                parsedInt = parsedInt == null ? -1 : parsedInt;

                Program.GoldenPeaks.PointstoClusterIDs[count] = parsedInt;
                GoldenExamination.GoldenID[count] = parsedInt;
                GoldenExamination.GoldenLabel[count] = splits[GoldenLabelPosition];
                GoldenExamination.PeakPosition[count][0] = Double.parseDouble(splits[1]);
                GoldenExamination.PeakPosition[count][1] = Double.parseDouble(splits[2]);
                parsedInt = Ints.tryParse(splits[MclustPosition]);
                if (parsedInt == null) {
                    Double tryagain = Doubles.tryParse(splits[MclustPosition]);
                    if (tryagain == null) {
                        DAVectorUtility.printAndThrowRuntimeException("Count " + count + "Illegal Mclust value on " +
                                "Point file " + splits.length + " " + line);
                    } else {
                        parsedInt = (int) Math.floor(tryagain + 0.001);
                    }
                }
                Program.MclustClusters.PointstoClusterIDs[count] = parsedInt;
                Program.ExperimentNumberAssigments[count] = Integer.valueOf(splits[ExpermentNumberPosition]);
                parsedInt = Ints.tryParse(splits[MedeaPosition]);
                if (parsedInt == null) {
                    Double tryagain = Doubles.tryParse(splits[MedeaPosition]);
                    if (tryagain == null) {
                        DAVectorUtility.printAndThrowRuntimeException("Count " + count + "Illegal Medea value on " +
                                "Point file " + splits.length + " " + line);
                    } else {
                        parsedInt = (int) Math.floor(tryagain + 0.001);
                    }
                }
                Program.MedeaClusters.PointstoClusterIDs[count] = parsedInt;
                count++;
            }
            success = true;
            br.close();
        } catch (IOException e) {
            System.err.format("Failed reading Points data due to I/O exception: %s%n", e);
        }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File Analyze error " + fname);
        }
        Program.GoldenPeaks.setup();
        Program.MclustClusters.setup();
        Program.MedeaClusters.setup();

    } // End ReadLabelsFromFile

    // read point data from file to memory
    // OutputFileType = 0 Harvard Format Cosmic Index mz RT Charge Peptide Peptide.id Mclust Medea
    // OutputFileType = 1 Ingest output file:   Index mz RT  0.0 Cluster#
    public static void AnalyzeDataFromFile(String fname) {
        int Maxcounts = 200;
        int[] CountLabels = new int[Maxcounts];
        for (int ChargeIndex = 0; ChargeIndex < Maxcounts; ChargeIndex++) {
            CountLabels[ChargeIndex] = 0;
        }
        int MinSplitSize = Program.ParameterVectorDimension + 3;
        int SplitPosition = 1 + Program.ParameterVectorDimension;
        if (Program.InputFileType == 1) {
            MinSplitSize = 5;
            SplitPosition = 4;
        }
        boolean success = false;
        int count = 0;

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try {
            BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset());
            String line;
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.printAndThrowRuntimeException("Count " + count + "Illegal data length on Point " +
                            "file " + splits.length + " " + MinSplitSize + " " + line);
                }

                Integer parsedInt = Ints.tryParse(splits[SplitPosition]);
                if (parsedInt == null || ((Program.SelectedInputLabel < 0) && (parsedInt == -Program
                        .SelectedInputLabel))) {
                    continue;
                }
                int position = 1 + parsedInt;
                position = Math.max(0, position);
                position = Math.min(Maxcounts - 1, position);
                ++CountLabels[position];
                count++;
            }
            success = true;
            if (Program.SelectedInputLabel >= 0) {
                DAVectorUtility.PointCount_Global = CountLabels[1 + Program.SelectedInputLabel];
            } else {
                DAVectorUtility.PointCount_Global = 0;
                for (int selection = 0; selection < Maxcounts; selection++) {
                    DAVectorUtility.PointCount_Global += CountLabels[selection];
                }
            }
            DAVectorUtility.SALSAPrint(1, "File Analyzed " + fname + " Total " + count);
            String label = "Charge";
            if (Program.InputFileType == 1) {
                label = "Cluster";
            }
            if (CountLabels[0] > 0) {
                DAVectorUtility.SALSAPrint(1, "Negative " + label + "s " + CountLabels[0]);
            }
            for (int LabelIndex = 1; LabelIndex < (Maxcounts - 1); LabelIndex++) {
                if (CountLabels[LabelIndex] > 0) {
                    DAVectorUtility.SALSAPrint(1, label + " " + (LabelIndex - 1) + " Count " + CountLabels[LabelIndex]);
                }
            }
            if (CountLabels[Maxcounts - 1] > 0) {
                DAVectorUtility.SALSAPrint(1, "Overflow " + label + "s " + CountLabels[Maxcounts - 1]);
            }
            br.close();
        } catch (IOException e) {
            System.err.format("Failed reading Points data due to I/O exception: %s%n", e);
        }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File Analyze error " + fname);
        }
    } // End AnalyzeDataFromFile

    //  Readvectors = 0 Normal
    //  Readvectors = 1 Read points changing cluster labels
    //  Readvectors = 2 Read Cluster Centers
    public static void ReadDataFromFile(String fname, int ReadVectorsOption) {
        int FirstPointPosition = 0;
        int TotalNumberPointstoRead = 0;
        if (ReadVectorsOption <= 1) {
            FirstPointPosition = DAVectorUtility.PointStart_Process;
            TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;
        }
        int BeginLabel = 0;
        if (ReadVectorsOption == 1) {
            BeginLabel = ParallelClustering.runningSolution.Ncent_Global - 1;
        }

        int MinSplitSize = Program.ParameterVectorDimension + 3;
        int SplitPosition = 1 + Program.ParameterVectorDimension;
        int LabelPosition = 4 + Program.ParameterVectorDimension;
        int ExpermentNumberPosition = 8;

        if (Program.InputFileType == 1) {
            MinSplitSize = 5;
            SplitPosition = 4;
            LabelPosition = 4;
        }
        boolean success = false;
        String line = " Unset";
        int CountLinesinFile = 0;
        int CurrentLocalPointtoUse = 0;

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }
    //    DAVectorUtility.SALSAPrint(0,"File Name " + fname + "ReadVectorsOption " + ReadVectorsOption + " rank " + DAVectorUtility.MPI_Rank);
        try {
            BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset());
            Pattern pattern = Pattern.compile("[\t ]");

            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.printAndThrowRuntimeException("Count " + CountLinesinFile + "Illegal data length " +
                            "on Point file " + splits.length + " " + MinSplitSize + " " + line);
                }

                Integer parsedInt = Ints.tryParse(splits[SplitPosition]);
                if (parsedInt == null) {
                    continue;
                }

                if (Program.SelectedInputLabel >= 0) {
                    if (parsedInt != Program.SelectedInputLabel) {
                        continue;
                    }
                } else {
                    if (ReadVectorsOption >= 2) {
                        if (parsedInt != -Program.SelectedInputLabel) {
                            continue;
                        }
                    } else {
                        if (parsedInt == -Program.SelectedInputLabel) {
                            continue;
                        }
                    }
                }


                if ((ReadVectorsOption <= 0) && (CountLinesinFile < FirstPointPosition)) {
                    CountLinesinFile += Program.Replicate;
                    continue;
                }

                int ActualPointPosition = 0;
                if (ReadVectorsOption == 0) {
                    ActualPointPosition = CountLinesinFile - FirstPointPosition;
                } else if (ReadVectorsOption == 1) {
                    double CurrentLineY0 = Double.parseDouble(splits[1]);
                    boolean enditall = false;
                    boolean endinputline = false;
                    while (true) {
                        if (CurrentLocalPointtoUse >= TotalNumberPointstoRead) {
                            enditall = true;
                            break;
                        }
                        double CurrentPointY0 = Program.PointPosition[CurrentLocalPointtoUse][0];
                        if (CurrentLineY0 < (CurrentPointY0 - 0.000001)) {
                            endinputline = true;
                            break;
                        }
                        if (CurrentLineY0 > (CurrentPointY0 + 0.000001)) {
                            ++CurrentLocalPointtoUse;
                        } else {
                            if (Math.abs(Program.PointPosition[CurrentLocalPointtoUse][1] - Double.parseDouble
                                    (splits[2])) > 0.00001) {
                                ++CurrentLocalPointtoUse;
                                continue;
                            }
                            if (Program.PointLabel[CurrentLocalPointtoUse] != 0) {
                                DAVectorUtility.printAndThrowRuntimeException(fname + " Position " +
                                        (CurrentLocalPointtoUse + FirstPointPosition) + " Inconsistent Cluster Number" +
                                        " " + Program.PointLabel[CurrentLocalPointtoUse] + " Y0 " + String
                                        .format("%1$6.5f", Program.PointPosition[CurrentLocalPointtoUse][0]) + "" +
                                        " " + line);
                            }
                            ActualPointPosition = CurrentLocalPointtoUse;
                            ++CurrentLocalPointtoUse;
                            break;
                        }
                    }
                    if (endinputline) {
                        continue;
                    }
                    if (enditall) {
                        break;
                    }
                }
                for (int CountReplicas = 0; CountReplicas < Program.Replicate; CountReplicas++) {
                    if (ReadVectorsOption >= 2) {
                        ParallelClustering.runningSolution.Y_k_i_[ParallelClustering.runningSolution.Ncent_Global][0]
                                = Double.parseDouble(splits[1]);
                        ParallelClustering.runningSolution.Y_k_i_[ParallelClustering.runningSolution.Ncent_Global][1]
                                = Double.parseDouble(splits[2]);
                        if (Program.ParameterVectorDimension > 2) {
                            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension;
                                 VectorIndex++) {
                                ParallelClustering.runningSolution.Y_k_i_[ParallelClustering.runningSolution
                                        .Ncent_Global][VectorIndex] = Double.parseDouble(splits[VectorIndex + 1]);
                            }
                        }
                        ++ParallelClustering.runningSolution.Ncent_Global;
                    } else {
                        if (ReadVectorsOption == 0) {
                            Program.PointPosition[ActualPointPosition][0] = Double.parseDouble(splits[1]);
                            Program.PointPosition[ActualPointPosition][1] = Double.parseDouble(splits[2]);
                            Program.PointOriginalIndex[ActualPointPosition] = Integer.parseInt(splits[0]);
                            Program.PointOriginalExperimentNumber[ActualPointPosition] = Integer.parseInt(splits[ExpermentNumberPosition]);
                            if (Program.ParameterVectorDimension > 2) {
                                for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension;
                                     VectorIndex++) {
                                    Program.PointPosition[ActualPointPosition][VectorIndex] = Double
                                            .parseDouble(splits[VectorIndex + 1]);
                                }
                            }
                        }

                        Double parsedDouble = Doubles.tryParse(splits[LabelPosition]);
                        if (parsedDouble == null) {
                            parsedDouble = 0.0;
                        }
                        Program.PointLabel[ActualPointPosition] = (int) (parsedDouble + 0.001);
                        if (Program.PointLabel[ActualPointPosition] != 0) {
                            Program.PointLabel[ActualPointPosition] += BeginLabel;
                        }
                    }
                    ++ActualPointPosition;
                }
                CountLinesinFile += Program.Replicate;
                if ((ReadVectorsOption <= 0) && (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead))) {
                    break;
                }

            }
            if ((ReadVectorsOption <= 0) && (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead))) {
                DAVectorUtility.printAndThrowRuntimeException("Illegal count on Points file " + fname + " Rank " +
                        DAVectorUtility.MPI_Rank + " Lines in File " + CountLinesinFile + " Number to Read " +
                        TotalNumberPointstoRead);
            }
            success = true;
            br.close();
        } catch (Exception e) {
            DAVectorUtility.printAndThrowRuntimeException("Failed reading Points data " + DAVectorUtility.MPI_Rank +
                    " " + CountLinesinFile + " Start " + FirstPointPosition + " Number " + TotalNumberPointstoRead +
                    " " + line + e.getMessage());
            }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File read error " + fname + " rank " + DAVectorUtility.MPI_Rank);
        }
    } // End ReadDataFromFile

    public static void ReadDataFromFile(String fname, int ClusterPosition, int[] InitialPointAssignment,
                                        int StartPointPosition) {
        int FirstPointPosition;
        int TotalNumberPointstoRead;
        FirstPointPosition = DAVectorUtility.PointStart_Process;
        TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;

        int MinSplitSize = ClusterPosition + 1;
        if (StartPointPosition >= 0) {
            MinSplitSize = Math.max(MinSplitSize, StartPointPosition + Program.ParameterVectorDimension);
        }

        boolean success = false;
        String line = " Unset";
        int CountLinesinFile = 0;
        String stringtest = "";

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try {
            BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset());
            Pattern pattern = Pattern.compile("[\t ,]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.SALSAPrint(0, "Count " + CountLinesinFile + " Illegal data length on Point file "
                            + splits.length + " " + MinSplitSize + " " + line);
                    continue;
                } // Skip header lines

                Double junk = Doubles.tryParse(splits[StartPointPosition]);
                if (junk == null) continue; // Skip header lines

                if (CountLinesinFile < FirstPointPosition) {
                    CountLinesinFile += 1;
                    continue;
                }

                int ActualPointPosition = CountLinesinFile - FirstPointPosition;

                if (StartPointPosition >= 0) {
                    stringtest = "0 *" + splits[StartPointPosition];
                    Program.PointPosition[ActualPointPosition][0] = Double.parseDouble
                            (splits[StartPointPosition]);
                    stringtest = "1 *" + splits[StartPointPosition + 1];
                    Program.PointPosition[ActualPointPosition][1] = Double.parseDouble
                            (splits[StartPointPosition + 1]);
                    if (Program.ParameterVectorDimension > 2) {
                        for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension;
                             VectorIndex++) {
                            stringtest = VectorIndex + " *" + splits[StartPointPosition + VectorIndex];
                            Program.PointPosition[ActualPointPosition][VectorIndex] = Double.parseDouble
                                    (splits[VectorIndex + StartPointPosition]);
                        }
                    }
                }
                if (ClusterPosition >= 0) {
                    Integer label = Ints.tryParse(splits[ClusterPosition]);
                    if (label == null) {
                        label = 1;
                    }
                    InitialPointAssignment[ActualPointPosition] = label - 1;
                }

                ++ActualPointPosition;
                ++CountLinesinFile;
                if (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead)) {
                    break;
                }

            }
            if (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead)) {
                DAVectorUtility.printAndThrowRuntimeException("Illegal count on Points file " + fname + " Rank " +
                        DAVectorUtility.MPI_Rank + " Lines in File " + CountLinesinFile + " Number to Read " +
                        TotalNumberPointstoRead);
            }
            success = true;
            br.close();
        } catch (Exception e) {
            DAVectorUtility.printAndThrowRuntimeException(stringtest + "* Failed reading Points data " +
                    DAVectorUtility.MPI_Rank + " " + CountLinesinFile + " Start " + FirstPointPosition + " Number " +
                    TotalNumberPointstoRead + " " + line + e.getMessage());
        }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File read error " + fname);
        }
    } // End ReadDataFromFile

    public static void Read3DDataFromFile(String fname, int VectorSize, int StartPointPosition) {
        int FirstPointPosition = 0;
        int TotalNumberPointstoRead;
        TotalNumberPointstoRead = DAVectorUtility.PointCount_Global;

        int MinSplitSize = StartPointPosition + VectorSize;

        boolean success = false;
        String line = " Unset";
        int CountLinesinFile = 0;

        if (Strings.isNullOrEmpty(fname)) {
            DAVectorUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        try {
            BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset());
            Pattern pattern = Pattern.compile("[\t ,]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());
                if (splits.length < MinSplitSize) {
                    DAVectorUtility.SALSAPrint(0, "Count " + CountLinesinFile + " Illegal data length on Point file "
                            + splits.length + " " + MinSplitSize + " " + line);
                    continue;
                } // Skip header lines

                Double junk = Doubles.tryParse(splits[StartPointPosition]);
                if (junk == null) continue; // Skip header lines

                if (CountLinesinFile < FirstPointPosition) {
                    CountLinesinFile += 1;
                    continue;
                }

                int ActualPointPosition = CountLinesinFile - FirstPointPosition;
                if (StartPointPosition >= 0) {
                    Program.FullPoint3DPosition[ActualPointPosition][0] = Double.parseDouble
                            (splits[StartPointPosition]);
                    Program.FullPoint3DPosition[ActualPointPosition][1] = Double.parseDouble
                            (splits[StartPointPosition + 1]);
                    if (VectorSize > 2) {
                        for (int VectorIndex = 2; VectorIndex < VectorSize; VectorIndex++) {
                            Program.FullPoint3DPosition[ActualPointPosition][VectorIndex] = Double.parseDouble
                                    (splits[VectorIndex + StartPointPosition]);
                        }
                    }
                }

                ++ActualPointPosition;
                ++CountLinesinFile;
                if (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead)) {
                    break;
                }

            }
            if (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead)) {
                DAVectorUtility.printAndThrowRuntimeException("Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank + " Lines in File " + CountLinesinFile + " Number to Read " + TotalNumberPointstoRead);
            }
            success = true;
            br.close();
        } catch (Exception e) {
            DAVectorUtility.printAndThrowRuntimeException("Failed reading Points data " + DAVectorUtility.MPI_Rank + " " + CountLinesinFile + " Start " + FirstPointPosition + " Number " + TotalNumberPointstoRead + " " + line + e.getMessage());
        }
        if (!success) {
            DAVectorUtility.printAndThrowRuntimeException("DA Vector File read error " + fname);
        }
    } // End Read3DDataFromFile
} // End class DAVectorReadData