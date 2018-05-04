using System;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RES.CSVReader;
using System.IO;
using Clerc2011;
using System.Globalization;
using Microsoft.SolverFoundation.Common;
using Microsoft.SolverFoundation.Services;

namespace MastWeightOptimizer
{
    class QuasiMWO
    {
        private static readonly CultureInfo _culture = new CultureInfo("en-GB");

        static void Main(string[] args)
        {
/****** To Change the Stack Size of the .exe change the value of the stackSize Variable***********************************************************************************************************************************************************/
            const int stackSize = 419430400;
            Thread th = new Thread(() =>
            {
/**************************************************************************************************************************************************************************************************************************************************/
/************************Part 1: Declare and Initialize all Variables and Read Data from CSVs *************************************************************************************************************************************************/
                /******** Indicate the batch size by changing MaxTurbines***********************/
                int MaxTurbines = 10;
                /*******************************************************************************/

                CSVReader mastReader = new CSVReader(new FileStream("MastLocations.csv", FileMode.Open));
                var value = int.Parse(mastReader["MastID"][0]);

                CSVReader turbineReader = new CSVReader(new FileStream("TurbineLocationsAndMastWeights.csv", FileMode.Open));

                var masts = new Mast[mastReader.RowCount];
                var mastIDs = new int[mastReader.RowCount];

                var turbines = new Turbine[turbineReader.RowCount];
                var turbineIDs = new int[turbineReader.RowCount];
                var Weight = new double[turbineReader.RowCount, mastReader.RowCount];

                var _binTopoUncertaintyCalculator = new BinTopoUncertaintyCalculator(new CoefficientOfVariationCalculator(new MastToTurbineVariationCalculator(0.2, 3000), new SpeedupVariationCalculator(0.5)), new CorrelationCalculator(new AngleCorrelationCalculator(), new DistanceCorrelationCalculator(1000)));

                var binReader = new CSVReader(new FileStream("BinInformation.csv", FileMode.Open));

                int count = binReader.RowCount;
                int number = turbineReader.RowCount * mastReader.RowCount;

                BinInformation[] bins = new BinInformation[binReader.RowCount];

                double[] variance = new double[count];

                int[] T = new int[binReader.RowCount];
                int[] M = new int[binReader.RowCount];

                int[] IT = new int[binReader.RowCount];
                int[] IM = new int[binReader.RowCount];

                for (int i = 0; i < mastReader.RowCount; i++)
                {
                    masts[i] = new Mast(double.Parse(mastReader["MastX"][i], _culture), double.Parse(mastReader["MastY"][i], _culture));
                    mastIDs[i] = int.Parse(mastReader["MastID"][i]);
                }

                for (int i = 0; i < turbineReader.RowCount; i++)
                {
                    turbines[i] = new Turbine(double.Parse(turbineReader["TurbineX"][i], _culture), double.Parse(turbineReader["TurbineY"][i], _culture));
                    turbineIDs[i] = int.Parse(turbineReader["TurbineID"][i]);
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        Weight[i, j] = double.Parse(turbineReader["Mast" + (j + 1) + "Weight"][i]);

                    }
                }

                for (int i = 0; i < binReader.RowCount; i++)
                {
                    int mastID = int.Parse(binReader["MastID"][i]);
                    int turbineID = int.Parse(binReader["TurbineID"][i]);

                    int mastIndex = GetIndex(mastIDs, mastID);
                    int turbineIndex = GetIndex(turbineIDs, turbineID);

                    bins[i].BinEnergy = double.Parse(binReader["Energy"][i], _culture);
                    bins[i].Direction = double.Parse(binReader["Direction"][i], _culture);
                    bins[i].SpeedUp = double.Parse(binReader["SpeedUp"][i], _culture);
                    bins[i].SensitivityFactor = double.Parse(binReader["Sensitivity"][i], _culture);
                    bins[i].Turbine = turbines[turbineIndex];
                    bins[i].Mast = masts[mastIndex];
                    T[i] = int.Parse(binReader["TurbineID"][i], _culture);
                    M[i] = int.Parse(binReader["MastID"][i], _culture);
                }

                double[,] MastWeights = new double[turbineReader.RowCount, mastReader.RowCount];

                int step = 0;

                SolverContext context = new SolverContext();

/**********************************Part 2: Define the solving algorithm and its parameters*********************************************************************************************************/
                HybridLocalSearchDirective MastWeightOptimizer = new HybridLocalSearchDirective();

                /*******Update these solver parameters by changing these variables******************/
                MastWeightOptimizer.RunUntilTimeout = false;
                // MastWeightOptimizer.TimeLimit = 60000 * 60 * 8;
                MastWeightOptimizer.EqualityTolerance = 0.001;
                MastWeightOptimizer.WaitLimit = -1;
                /*************************************************************************************/

/*****************************************Part 3: Begin the batch optimisation Process *******************************************************************************/
                for (int x = 0; x < turbineReader.RowCount / MaxTurbines; x++)
                {
/**********************************Part 3.1: Create Batch Solver Model, Decisions and Constraints-Re-initialise each time **********************************************************************************/
                    Term total = 0;
                    Term Energy = 0;
                    Term Uncertainty = 0;

                    Model model = context.CreateModel();

                    Decision[,] MWeight = new Decision[MaxTurbines, mastReader.RowCount];

                    int y = x * MaxTurbines;
                    for (int row = 0; row < MaxTurbines; row++)
                    {
                        for (int col = 0; col < mastReader.RowCount; col++)
                        {
                            MWeight[row, col] = new Decision(Domain.RealNonnegative, "MastWeight" + x + (row+1) + (col+1));
                            model.AddDecision(MWeight[row, col]);
                            model.AddConstraint("limit" + row + col, 0 <= MWeight[row, col] <= 1);
                            MWeight[row, col].SetInitialValue(Weight[(y + row), col]);
                        }
                    }

                    /******* Sum of mast Weights must equal 1 constraint****************************/
                    for (int i = 0; i < MaxTurbines; i++)
                    {
                        Term balance = 0;
                        Term holder = 0;
                        for (int j = 0; j < mastReader.RowCount; j++)
                        {
                            holder = MWeight[i, j];
                            balance += holder;
                        }

                        model.AddConstraint("Balances" + x + i, balance == 1);
                    }
                    /********************************************************************************/

                    BinInformation[] Ibins = new BinInformation[MaxTurbines * 12 * mastReader.RowCount];

/************************Part 3.2: Create batch bin information as subset of full bin information CSV*******************************************************************************************/

                    int TurbID = 0;
                    int tracker = 0;

                    for (int i = 0; i < binReader.RowCount; i++)
                    {
                        TurbID = int.Parse(binReader["TurbineID"][i], _culture);

                        if ((x * MaxTurbines) < TurbID && ((x + 1) * MaxTurbines) >= TurbID)
                        {
                            int mastID = int.Parse(binReader["MastID"][i]);
                            int turbineID = int.Parse(binReader["TurbineID"][i]);

                            int mastIndex = GetIndex(mastIDs, mastID);
                            int turbineIndex = GetIndex(turbineIDs, turbineID);

                            Ibins[tracker].BinEnergy = double.Parse(binReader["Energy"][i], _culture);
                            Ibins[tracker].Direction = double.Parse(binReader["Direction"][i], _culture);
                            Ibins[tracker].SpeedUp = double.Parse(binReader["SpeedUp"][i], _culture);
                            Ibins[tracker].SensitivityFactor = double.Parse(binReader["Sensitivity"][i], _culture);
                            Ibins[tracker].Turbine = turbines[turbineIndex];
                            Ibins[tracker].Mast = masts[mastIndex];
                            IT[tracker] = int.Parse(binReader["TurbineID"][i], _culture);
                            IM[tracker] = int.Parse(binReader["MastID"][i], _culture);

                            tracker++;

                        }
                    }

/***************************Part 3.3: Define the Objective Function: I.e. Uncertainty************************************************************************************************************/

                    for (int i = 0; i < Ibins.Length; i++)
                    {
                        int mast1 = IM[i] - 1;
                        int turbine1 = (IT[i] - 1) - (x * 10);

                        if (Ibins[i].BinEnergy > 0 && Ibins[i].SensitivityFactor > 0)
                        {
                            variance[i] = _binTopoUncertaintyCalculator.CalculateVariance(Ibins[i]);

                            if (variance[i] > 0)
                            {
                                for (int j = 0; j <= i; j++)
                                {
                                    int mast2 = IM[j] - 1;
                                    int turbine2 = (IT[j] - 1) - (x * 10);

                                    if (Ibins[j].BinEnergy > 0 && Ibins[j].SensitivityFactor > 0 && variance[j] > 0)
                                    {

                                        var singleResult = _binTopoUncertaintyCalculator.Calculate(Ibins[i], Ibins[j], variance[i], variance[j]) * Ibins[i].BinEnergy * Ibins[j].BinEnergy * MWeight[turbine1, mast1] * MWeight[turbine2, mast2];

                                        if (i != j)
                                            singleResult *= 2;

                                        total += singleResult;
                                    }
                                }
                            }
                        }
                        Energy += Ibins[i].BinEnergy * MWeight[turbine1, mast1];

                    }

                    total = Model.Sqrt(total);

                    Uncertainty = total / Energy;

/**************************************Part 3.4: Solve the Model ********************************************************************************************************************/
                    model.AddGoal("Uncertainty", GoalKind.Minimize, Uncertainty);

                    Solution solution = context.Solve(MastWeightOptimizer);

                    Report report = solution.GetReport();

/**********************************Part 3.5: Print Batch Report to the Console *******************************************************************************************************/
                    Console.WriteLine("Working ");
                    Console.WriteLine("{0}", report);

/*************Part 3.6: Convert optimised mast-to-turbine weightings to double floating points and append to array**********************************************************************************/

                    for (int i = 0; i < MaxTurbines; i++)
                    {
                        for (int j = 0; j < mastReader.RowCount; j++)
                        {
                            MastWeights[(i + (y)), j] = MWeight[i, j].ToDouble();
                        }
                    }

/***************************par 3.7: Clear the model and begin the batch optimisation again ***************************************************************************/

                    Console.WriteLine("y: {0}", y);
                    context.ClearModel();
                    step++;
                }

/**********************************************************************************************************************************************************************/
/*********************************Part 4: Repeat the batch process for the remaining turbines less than Max Turbines **************************************************/

                int OtherTurbines = turbineReader.RowCount % MaxTurbines;
                Console.WriteLine("{0}", OtherTurbines);
                Console.WriteLine("{0}", step);

                Term Othertotal = 0;
                Term OtherEnergy = 0;
                Term OtherUncertainty = 0;

                SolverContext Othercontext = SolverContext.GetContext();
                Model Othermodel = Othercontext.CreateModel();

                Decision[,] OtherMWeight = new Decision[OtherTurbines, mastReader.RowCount];

                for (int row = 0; row < OtherTurbines; row++)
                {
                    for (int col = 0; col < mastReader.RowCount; col++)
                    {
                        OtherMWeight[row, col] = new Decision(Domain.RealNonnegative, "MastWeight" + step + (row+1) + (col+1));
                        Othermodel.AddDecision(OtherMWeight[row, col]);
                        Othermodel.AddConstraint("limit" + row + col, 0 <= OtherMWeight[row, col] <= 1);
                        OtherMWeight[row, col].SetInitialValue(Weight[((step * MaxTurbines) + row), col]);
                    }
                }

                for (int i = 0; i < OtherTurbines; i++)
                {
                    Term Otherbalance = 0;
                    Term Otherholder = 0;
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        Otherholder = OtherMWeight[i, j];
                        Otherbalance += Otherholder;
                    }

                    Othermodel.AddConstraint("Balances" + step + i, Otherbalance == 1);
                }

                BinInformation[] OtherIbins = new BinInformation[OtherTurbines * 12 * mastReader.RowCount];

                int OtherTurbID = 0;
                int Othertracker = 0;

                int[] OtherIT = new int[binReader.RowCount];
                int[] OtherIM = new int[binReader.RowCount];

                for (int i = 0; i < binReader.RowCount; i++)
                {
                    OtherTurbID = int.Parse(binReader["TurbineID"][i], _culture);
                    if ((step * MaxTurbines) < OtherTurbID && ((step * MaxTurbines) + OtherTurbines) >= OtherTurbID)
                    {
                        int mastID = int.Parse(binReader["MastID"][i]);
                        int turbineID = int.Parse(binReader["TurbineID"][i]);

                        int mastIndex = GetIndex(mastIDs, mastID);
                        int turbineIndex = GetIndex(turbineIDs, turbineID);

                        OtherIbins[Othertracker].BinEnergy = double.Parse(binReader["Energy"][i], _culture);
                        OtherIbins[Othertracker].Direction = double.Parse(binReader["Direction"][i], _culture);
                        OtherIbins[Othertracker].SpeedUp = double.Parse(binReader["SpeedUp"][i], _culture);
                        OtherIbins[Othertracker].SensitivityFactor = double.Parse(binReader["Sensitivity"][i], _culture);
                        OtherIbins[Othertracker].Turbine = turbines[turbineIndex];
                        OtherIbins[Othertracker].Mast = masts[mastIndex];
                        OtherIT[Othertracker] = int.Parse(binReader["TurbineID"][i], _culture);
                        OtherIM[Othertracker] = int.Parse(binReader["MastID"][i], _culture);

                        Othertracker++;
                    }
                }

                Console.WriteLine("{0}", Othertracker);

                for (int i = 0; i < OtherIbins.Length; i++)
                {
                    int mast1 = OtherIM[i] - 1;
                    int turbine1 = (OtherIT[i] - 1) - (step * MaxTurbines);

                    if (OtherIbins[i].BinEnergy > 0 && OtherIbins[i].SensitivityFactor > 0)
                    {
                        variance[i] = _binTopoUncertaintyCalculator.CalculateVariance(OtherIbins[i]);

                        if (variance[i] > 0)
                        {
                            for (int j = 0; j <= i; j++)
                            {
                                int mast2 = OtherIM[j] - 1;
                                int turbine2 = (OtherIT[j] - 1) - (step * MaxTurbines);

                                if (OtherIbins[j].BinEnergy > 0 && OtherIbins[j].SensitivityFactor > 0 && variance[j] > 0)
                                {

                                    var singleResult = _binTopoUncertaintyCalculator.Calculate(OtherIbins[i], OtherIbins[j], variance[i], variance[j]) * OtherIbins[i].BinEnergy * OtherIbins[j].BinEnergy * OtherMWeight[turbine1, mast1] * OtherMWeight[turbine2, mast2];

                                    if (i != j)
                                        singleResult *= 2;

                                    Othertotal += singleResult;
                                }
                            }
                        }
                    }
                    OtherEnergy += OtherIbins[i].BinEnergy * OtherMWeight[turbine1, mast1];

                }

                Othertotal = Model.Sqrt(Othertotal);
                OtherUncertainty = Othertotal / OtherEnergy;

                Othermodel.AddGoal("Uncertainty", GoalKind.Minimize, OtherUncertainty);

                Solution Othersolution = Othercontext.Solve(MastWeightOptimizer);

                Report Otherreport = Othersolution.GetReport();

                Console.WriteLine("Working ");
                Console.WriteLine("{0}", Otherreport);

                for (int i = 0; i < OtherTurbines; i++)
                {
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        MastWeights[(i + (step * MaxTurbines)), j] = OtherMWeight[i, j].ToDouble();
                    }
                }

/************Part 5: Using optimised mast-to-turbine weightings from batches calculate overall uncertainty, energy and sigma*************************************************************************************************************************************************************************************************************/

                double FinalEnergy = 0;
                double FinalSigma = 0;
                double FinalUncertainty = 0;

                var DubBinTopoUncertaintyCalculator = new BinTopoUncertaintyCalculator(new CoefficientOfVariationCalculator(new MastToTurbineVariationCalculator(0.2, 3000), new SpeedupVariationCalculator(0.5)), new CorrelationCalculator(new AngleCorrelationCalculator(), new DistanceCorrelationCalculator(1000)));

                for (int i = 0; i < count; i++)
                {
                    int mast1 = M[i] - 1;
                    int turbine1 = T[i] - 1;

                    if (bins[i].BinEnergy > 0 && bins[i].SensitivityFactor > 0)
                    {
                        variance[i] = DubBinTopoUncertaintyCalculator.CalculateVariance(bins[i]);

                        if (variance[i] > 0)
                        {
                            for (int j = 0; j <= i; j++)
                            {
                                int mast2 = M[j] - 1;
                                int turbine2 = T[j] - 1;

                                if (bins[j].BinEnergy > 0 && bins[j].SensitivityFactor > 0 && variance[j] > 0)
                                {

                                    double singleResult = DubBinTopoUncertaintyCalculator.Calculate(bins[i], bins[j], variance[i], variance[j]) * bins[i].BinEnergy * bins[j].BinEnergy * MastWeights[turbine1, mast1] * MastWeights[turbine2, mast2];

                                    if (i != j)
                                        singleResult *= 2;

                                    FinalSigma += singleResult;
                                }
                            }
                        }
                    }

                    FinalEnergy += bins[i].BinEnergy * MastWeights[turbine1, mast1];
                }

                FinalSigma = Math.Sqrt(FinalSigma);

                FinalUncertainty = FinalSigma / FinalEnergy;

/***************Part 6: Print the optimised mast-to-turbine weightings, uncertainty, energy and sigma to the console*************************************************************************************************************************/

                Console.WriteLine("************* Optimised Mast-To-Turbine Weihghtings**********************************************");
                Console.WriteLine("");
                Console.Write("Turbine#     ");
                for (int i = 0; i < mastReader.RowCount; i++)
                {
                    Console.Write("MastWeight{0}    ", (i + 1));
                }

                for (int i = 0; i < turbineReader.RowCount; i++)
                {
                    Console.WriteLine("");
                    Console.Write("{0}              ", (i + 1));
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        Console.Write("{0:0.00}             ", MastWeights[i, j]);
                    }
                }

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("************** Optimised Uncertainty, Energy and Sigma*********************************************");
                Console.WriteLine("");
                Console.WriteLine("Energy: {0}", FinalEnergy);
                Console.WriteLine("Sigma: {0}", FinalSigma);
                Console.WriteLine("");
                Console.WriteLine("Uncertainty: {0}", FinalUncertainty);
                Console.WriteLine("");
                Console.WriteLine("****************************************************************************************************");


                int GetIndex(int[] ids, int id)
                {
                    for (int i = 0; i < ids.Length; i++)
                    {
                        if (ids[i] == id) return i;
                    }
                    return -1;

                }
            }, stackSize);
            th.Start();
            th.Join();
        }

    }
}


